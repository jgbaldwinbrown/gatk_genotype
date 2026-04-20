package main

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os/exec"
	"os"
	"flag"
	"log"
	"encoding/csv"

	"golang.org/x/sync/errgroup"
)

type Flags struct {
	NoAln bool
	RefPath string
	SeqPairsPath string
	Outpre string
	Threads int
	MemoryGb int
	Nproc int
	Trim bool
	BamPathsPath string
	DeleteTempFiles bool
}

const adaptersFa = `>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC>PrefixPE/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA>TruSeq2_SE
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_f
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_PE_r
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA`

func CreateTempAdapters() (string, error) {
	fp, e := os.CreateTemp("", "adapters*.fa")
	if e != nil {
		return "", e
	}
	path := fp.Name()
	if _, e := fmt.Fprintf(fp, "%v\n", adaptersFa); e != nil {
		os.Remove(path)
		return "", e
	}
	if e := fp.Close(); e != nil {
		os.Remove(path)
		return "", e
	}
	return path, nil
}

func TrimCore(forward, reverse, forwardTrimmed, forwardTrimmedUnpaired, reverseTrimmed, reverseTrimmedUnpaired string, threads int) (err error) {
	adapterPath, e := CreateTempAdapters()
	if e != nil {
		return e
	}
	defer func() {
		e := os.Remove(adapterPath)
		if err == nil {
			err = e
		}
	}()
	clip := fmt.Sprintf("ILLUMINACLIP:%v:2:30:10", adapterPath)

	cmd := exec.Command(
		"java", "-jar", "/usr/share/java/trimmomatic-0.39.jar",
		"PE",
		"-threads", fmt.Sprint(threads),
		"-phred33",
		forward, reverse,
		forwardTrimmed, forwardTrimmedUnpaired,
		reverseTrimmed, reverseTrimmedUnpaired,
		clip,
		"LEADING:20",
		"TRAILING:20",
		"MINLEN:30",
	)
	cmd.Stderr = os.Stderr
	cmd.Stdout = os.Stdout
	return cmd.Run()
}

func Trim(forward, reverse, outpre string, threads int) error {
	return TrimCore(
		forward, reverse,
		outpre + "_forward_trimmed.fq.gz", outpre + "_forward_trimmed_unpaired.fq.gz",
		outpre + "_reverse_trimmed.fq.gz", outpre + "_reverse_trimmed_unpaired.fq.gz",
		threads,
	)
}

func FileExists(path string) bool {
	_, err := os.Stat(path)
	return err == nil
}

func DeleteTrim(f Flags, s ReadSet) error {
	outpre := f.Outpre + "_" + s.Name
	trimouts := []string {
		outpre + "_forward_trimmed.fq.gz",
		outpre + "_forward_trimmed_unpaired.fq.gz",
		outpre + "_reverse_trimmed.fq.gz",
		outpre + "_reverse_trimmed_unpaired.fq.gz",
	}
	for _, trimout := range trimouts {
		if FileExists(trimout) {
			if e := os.Remove(trimout); e != nil {
				return e
			}
		}
	}
	return nil
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

func DeleteBam(f Flags, s ReadSet) error {
	bampath := f.Outpre + "_" + s.Name + ".bam"
	if FileExists(bampath) {
		if e := os.Remove(bampath); e != nil {
			return e
		}
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

func DeleteRGBam(f Flags, s ReadSet) error {
	bampathrg := f.Outpre + "_" + s.Name + "_rg.bam"
	if FileExists(bampathrg) {
		if e := os.Remove(bampathrg); e != nil {
			return e
		}
	}
	return nil
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

func DeleteGvcfs(f Flags, s ReadSet) error {
	gvcfpath := f.Outpre + "_" + s.Name + ".g.vcf.gz"
	if FileExists(gvcfpath) {
		if e := os.Remove(gvcfpath); e != nil {
			return e
		}
	}
	return nil
}

func DeletePath(path string) error {
	if FileExists(path) {
		if e := os.Remove(path); e != nil {
			return e
		}
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
	BamPath string
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

func ParseBamSets(path string) (pairs []ReadSet, err error) {
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
		if len(l) < 2 {
			return nil, fmt.Errorf("ParseReadSets: len(l) %v < 2; l %v", len(l), l)
		}
		pairs = append(pairs, ReadSet{
			Name: l[0],
			BamPath: l[1],
		})
	}
	return pairs, nil
}

func FullFQFMimic(f Flags) (err error) {
	var sets []ReadSet
	if e := BwaIndex(f.RefPath); e != nil {
		return e
	}
	if e := Faidx(f.RefPath); e != nil {
		return e
	}
	if e := CreateDict(f.RefPath); e != nil {
		return e
	}
	if !f.NoAln {
		var e error
		sets, e = ParseReadSets(f.SeqPairsPath)
		if e != nil {
			return e
		}
	} else {
		var e error
		sets, e = ParseBamSets(f.BamPathsPath)
		if e != nil {
			return e
		}
	}

	var group errgroup.Group
	group.SetLimit(f.Nproc)
	var gvcfpaths []string
	for _, set := range sets {
		set := set

		bampath := f.Outpre + "_" + set.Name + ".bam"
		if f.NoAln {
			bampath = set.BamPath
		}

		bampathrg := f.Outpre + "_" + set.Name + "_rg.bam"
		gvcfpath := f.Outpre + "_" + set.Name + ".g.vcf.gz"
		gvcfpaths = append(gvcfpaths, gvcfpath)
		fwd := set.ForwardPath
		rev := set.ReversePath

		group.Go(func() error {
			if !f.NoAln {
				if f.Trim {
					fwd = f.Outpre + "_" + set.Name + "_forward_trimmed.fq.gz"
					rev = f.Outpre + "_" + set.Name + "_reverse_trimmed.fq.gz"
					if e := Trim(set.ForwardPath, set.ReversePath, f.Outpre + "_" + set.Name, f.Threads); e != nil {
						return e
					}
					if f.DeleteTempFiles {
						defer func() {
							errors.Join(err, DeleteTrim(f, set))
						}()
					}
				}
				if e := BwaMem(f.RefPath, fwd, rev, bampath, f.Threads); e != nil {
					return e
				}
				if f.DeleteTempFiles {
					defer func() {
						errors.Join(err, DeleteBam(f, set))
					}()
				}
			}
			if e := AddRG(bampath, bampathrg, set.Name); e != nil {
				return e
			}
			if f.DeleteTempFiles {
				defer func() {
					errors.Join(err, DeleteRGBam(f, set))
				}()
			}
			if e := SamIndex(bampathrg); e != nil {
				return e
			}
			if e := HaplotypeCall(f.RefPath, bampathrg, gvcfpath, set.Name, f.MemoryGb); e != nil {
				return e
			}
			return nil
		})
		if f.DeleteTempFiles {
			defer func() {
				errors.Join(err, DeleteGvcfs(f, set))
			}()
		}
	}
	if e := group.Wait(); e != nil {
		return e
	}

	goutpath := f.Outpre + ".g.vcf.gz"
	if f.DeleteTempFiles {
		defer func() {
			errors.Join(err, DeletePath(goutpath))
		}()
	}
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
	flag.IntVar(&f.Nproc, "n", 1, "Number of simultaneous runs of BWA / picard to run")
	flag.BoolVar(&f.Trim, "T", false, "Also trim input files with trimmomatic")
	flag.BoolVar(&f.NoAln, "noaln", false, "Skip aligning (it is already done). Not compatible with -s.")
	flag.BoolVar(&f.DeleteTempFiles, "d", false, "Delete all intermediate files except for index files and the final .vcf.gz file.")
	flag.StringVar(&f.BamPathsPath, "bams", "", "Path to file containing paths to bams, one line per sample. Format: name (tab) path.bam. Not compatible with -s and requires -noaln.")
	flag.Parse()

	if f.RefPath == "" {
		log.Fatal("missing -r")
	}

	if !f.NoAln {
		if (f.SeqPairsPath == "") {
			log.Fatal("missing -s")
		}
	} else {
		if (f.BamPathsPath == "") {
			log.Fatal("missing -bams")
		}
	}

	if e := FullFQFMimic(f); e != nil {
		log.Fatal(e)
	}
}
