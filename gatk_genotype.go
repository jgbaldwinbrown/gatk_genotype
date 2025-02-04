package main

import (
	"bufio"
	"fmt"
	"io"
	"os/exec"
	"os"
	"flag"
	"log"
	"encoding/csv"
)

type Flags struct {
	RefPath string
	SeqPairsPath string
	Outpre string
	Threads int
	MemoryGb int
}

func BwaIndex(ref string) error {
	cmd := exec.Command("bwa", "index", ref)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func BwaMem(ref, forward, reverse, out string, threads int) (err error) {
	w, err := os.Create(out)
	if err != nil {
		return err
	}
	defer func() {
		e := w.Close()
		if err == nil {
			err = e
		}
	}()
	bw := bufio.NewWriter(w)
	defer func() {
		e := bw.Flush()
		if err == nil {
			err = e
		}
	}()

	bwaCmd := exec.Command("bwa", "mem", ref, forward, reverse, "-t", fmt.Sprint(threads))
	bwaCmd.Stderr = os.Stderr
	bwaOut, err := bwaCmd.StdoutPipe()
	if err != nil {
		return err
	}

	viewCmd := exec.Command("samtools", "view", "-bS")
	viewCmd.Stdin = bwaOut
	viewCmd.Stderr = os.Stderr
	viewOut, err := viewCmd.StdoutPipe()
	if err != nil {
		return err
	}

	sortCmd := exec.Command("samtools", "sort")
	sortCmd.Stderr = os.Stderr
	sortCmd.Stdin = viewOut
	sortCmd.Stdout = bw

	if e := bwaCmd.Start(); e != nil {
		return e
	}
	if e := viewCmd.Start(); e != nil {
		return e
	}
	if e := sortCmd.Start(); e != nil {
		return e
	}

	if e := bwaCmd.Wait(); e != nil {
		return e
	}
	if e := viewCmd.Wait(); e != nil {
		return e
	}
	if e := sortCmd.Wait(); e != nil {
		return e
	}

	return nil
}

func AddRG(inpath, outpath, name string) error {
	rgchange := fmt.Sprintf("@RG\tID:%v\tSM:%v", name, name)
	cmd := exec.Command("samtools", "addreplacerg", "-r", rgchange, inpath, "-o", outpath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func Faidx(fapath string) error {
	cmd := exec.Command("samtools", "faidx", fapath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func SamIndex(bampath string) error {
	cmd := exec.Command("samtools", "index", bampath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func CreateDict(fapath string) error {
	cmd := exec.Command("picard-tools", "CreateSequenceDictionary", "-R", fapath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func HaplotypeCall(fapath, bampath, gvcfpath, name string, memoryGb int) error {
	memstr := fmt.Sprintf("-Xmx%vg", memoryGb)
	cmd := exec.Command(
		"gatk",
		"--java-options", memstr,
		"HaplotypeCaller", 
		"-R", fapath,
		"-I", bampath,
		"-O", gvcfpath,
		"--sample-name", name,
		"-ERC", "GVCF",
	)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if e := cmd.Run(); e != nil {
		return fmt.Errorf("HaplotypeCall: %w", e)
	}
	return nil
}

func CombineGvcf(fapath, gvcfoutpath, vcfoutpath string, memoryGb int, gvcfpaths ...string) error {
	memstr := fmt.Sprintf("-Xmx%vg", memoryGb)
	combsl := []string{
		"gatk",
		"--java-options", memstr,
		"CombineGVCFs",
		"-R", fapath,
		"-O", gvcfoutpath,
	}
	for _, gvcfpath := range gvcfpaths {
		combsl = append(combsl, "--variant", gvcfpath)
	}
	comb := exec.Command(combsl[0], combsl[1:]...)
	comb.Stdout = os.Stdout
	comb.Stderr = os.Stderr
	if e := comb.Run(); e != nil {
		return e
	}
	geno := exec.Command(
		"gatk",
		"--java-options", memstr,
		"GenotypeGVCFs",
		"-R", fapath,
		"-V", gvcfoutpath,
		"-O", vcfoutpath,
	)
	if e := geno.Run(); e != nil {
		return fmt.Errorf("GenotypeGVCFs: %w", e)
	}
	return nil
}

type ReadSet struct {
	Name string
	ForwardPath string
	ReversePath string
}

func ParseReadSets(path string) (pairs []ReadSet, err error) {
	r, e := os.Open(path)
	if e != nil {
		return nil, e
	}
	defer func() {
		e := r.Close()
		if err == nil {
			err = e
		}
	}()
	cr := csv.NewReader(r)
	cr.Comma = '\t'
	for l, e := cr.Read(); e != io.EOF; l, e = cr.Read() {
		if len(l) < 3 {
			return nil, fmt.Errorf("ParseReadSets: len(l) %v < 3; l %v", len(l), l)
		}
		pairs = append(pairs, ReadSet{
			Name: l[0],
			ForwardPath: l[1],
			ReversePath: l[2],
		})
	}
	return pairs, nil
}

func FullFQFMimic(f Flags) error {
	if e := BwaIndex(f.RefPath); e != nil {
		return e
	}
	if e := Faidx(f.RefPath); e != nil {
		return e
	}
	if e := CreateDict(f.RefPath); e != nil {
		return e
	}
	sets, e := ParseReadSets(f.SeqPairsPath)
	if e != nil {
		return e
	}

	var gvcfpaths []string
	for _, set := range sets {
		bampath := f.Outpre + "_" + set.Name + ".bam"
		if e := BwaMem(f.RefPath, set.ForwardPath, set.ReversePath, bampath, f.Threads); e != nil {
			return e
		}

		bampathrg := f.Outpre + "_" + set.Name + "_rg.bam"
		if e := AddRG(bampath, bampathrg, set.Name); e != nil {
			return e
		}
		if e := SamIndex(bampathrg); e != nil {
			return e
		}

		gvcfpath := f.Outpre + "_" + set.Name + ".g.vcf.gz"
		if e := HaplotypeCall(f.RefPath, bampathrg, gvcfpath, set.Name, f.MemoryGb); e != nil {
			return e
		}
		gvcfpaths = append(gvcfpaths, gvcfpath)
	}

	goutpath := f.Outpre + ".g.vcf.gz"
	outpath := f.Outpre + ".vcf.gz"
	return CombineGvcf(f.RefPath, goutpath, outpath, f.MemoryGb, gvcfpaths...)
	
}

func main() {
	var f Flags
	flag.StringVar(&f.RefPath, "r", "", "Path to reference .fa file (required)")
	flag.StringVar(&f.SeqPairsPath, "s", "", "Path to tab-separated table containing pairs of forward and reverse read paths, one line per sample (required). Format: name (tab) forward.fq.gz (tab) reverse.fq.gz")
	flag.StringVar(&f.Outpre, "o", "out", "Output prefix")
	flag.IntVar(&f.Threads, "t", 1, "Threads to use")
	flag.IntVar(&f.MemoryGb, "m", 8, "Memory to use (integer, gigabytes)")
	flag.Parse()
	if (f.RefPath == "" || f.SeqPairsPath == "") {
		log.Fatal("missing -r or -s")
	}

	if e := FullFQFMimic(f); e != nil {
		log.Fatal(e)
	}
}
