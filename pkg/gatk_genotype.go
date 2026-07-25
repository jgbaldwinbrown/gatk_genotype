package ggtype

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os/exec"
	"os"
	"encoding/csv"
	"strings"
	"regexp"

	"golang.org/x/sync/errgroup"
	"github.com/jgbaldwinbrown/zfile"
)

type Span struct {
	Chr string
	Start int64
	End int64
	FullChr bool
}

var extRe = regexp.MustCompile(`\.[^/]*$`)

func StripExtension(path string) string {
	return extRe.ReplaceAllString(path, "")
}

type Flags struct {
	NoAln bool
	NoHaploCall bool
	NoCombine bool
	RefPath string
	SeqPairsPath string
	Outpre string
	Threads int
	MemoryGb int
	SerialMemoryGb int
	Nproc int
	Trim bool
	BamPathsPath string
	DeleteTempFiles bool
	PicardCmd string
	Gogogo bool
	ChrsPath string
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

func Trim(forward, reverse, outpre string, threads int, gogogo bool) error {
	fwd := outpre + "_forward_trimmed.fq.gz"
	if !gogogo && IsDone(fwd) {
		return nil
	}

	err := TrimCore(
		forward, reverse,
		outpre + "_forward_trimmed.fq.gz", outpre + "_forward_trimmed_unpaired.fq.gz",
		outpre + "_reverse_trimmed.fq.gz", outpre + "_reverse_trimmed_unpaired.fq.gz",
		threads,
	)
	if err != nil {
		return err
	}

	return CreateDone(fwd)
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

	donepath := outpre + "_forward_trimmed.fq.gz" + ".done"
	if FileExists(donepath) {
		if e := os.Remove(donepath); e != nil {
			return e
		}
	}
	return nil
}

func BwaIndex(ref string, gogogo bool) error {
	if !gogogo && IsDone(ref + ".bwt") {
		return nil
	}

	cmd := exec.Command("bwa", "index", ref)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	err := cmd.Run()
	if err != nil {
		return err
	}

	return CreateDone(ref + ".bwt")
}

func BwaMem(ref, forward, reverse, out string, threads int, gogogo bool) (err error) {
	if !gogogo && IsDone(out) {
		return nil
	}

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

	return CreateDone(out)
}

func DeleteBam(f Flags, s ReadSet) error {
	bampath := f.Outpre + "_" + s.Name + ".bam"
	if FileExists(bampath) {
		if e := os.Remove(bampath); e != nil {
			return e
		}
	}
	donepath := bampath + ".done"
	if FileExists(donepath) {
		if e := os.Remove(donepath); e != nil {
			return e
		}
	}
	return nil
}

func AddRG(inpath, outpath, name string, gogogo bool) error {
	if !gogogo && IsDone(outpath) {
		return nil
	}

	rgchange := fmt.Sprintf("@RG\tID:%v\tSM:%v", name, name)
	cmd := exec.Command("samtools", "addreplacerg", "-r", rgchange, inpath, "-o", outpath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	err := cmd.Run()
	if err != nil {
		return err
	}

	return CreateDone(outpath)
}

func DeleteRGBam(f Flags, s ReadSet) error {
	bampathrg := f.Outpre + "_" + s.Name + "_rg.bam"
	if FileExists(bampathrg) {
		if e := os.Remove(bampathrg); e != nil {
			return e
		}
	}
	donepath := bampathrg + ".done"
	if FileExists(donepath) {
		if e := os.Remove(donepath); e != nil {
			return e
		}
	}
	return nil
}

func Faidx(fapath string, gogogo bool) error {
	if !gogogo && IsDone(fapath + ".fai") {
		return nil
	}

	cmd := exec.Command("samtools", "faidx", fapath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		return err
	}

	return CreateDone(fapath + ".fai")
}

func SamIndex(bampath string, gogogo bool) error {
	if !gogogo && IsDone(bampath + ".bai") {
		return nil
	}

	cmd := exec.Command("samtools", "index", bampath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if e := cmd.Run(); e != nil {
		return e
	}

	return CreateDone(bampath + ".bai")
}

func CreateDict(picardcmd, fapath string, gogogo bool) error {
	dictpath := StripExtension(fapath) + ".dict"
	if !gogogo && IsDone(dictpath) {
		return nil
	}

	args := strings.Fields(picardcmd)
	args = append(args, "CreateSequenceDictionary", "-R", fapath)
	cmd := exec.Command(args[0], args[1:]...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if e := cmd.Run(); e != nil {
		return e
	}

	return CreateDone(dictpath)
}

func HaplotypeCall(fapath, bampath, gvcfpath, name string, memoryGb int, gogogo bool) error {
	if !gogogo && IsDone(gvcfpath) {
		return nil
	}

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

	return CreateDone(gvcfpath)
}

func HaplotypeCallSpans(fapath, bampath, gvcfpath, name string, memoryGb int, gogogo bool, spans ...Span) error {
	if !gogogo && IsDone(gvcfpath) {
		return nil
	}

	memstr := fmt.Sprintf("-Xmx%vg", memoryGb)
	args := []string{
		"gatk",
		"--java-options", memstr,
		"HaplotypeCaller", 
		"-R", fapath,
		"-I", bampath,
		"-O", gvcfpath,
		"--sample-name", name,
		"-ERC", "GVCF",
	}
	for _, span := range spans {
		if span.FullChr {
			args = append(args, "-L", span.Chr)
		} else {
			args = append(args, "-L", fmt.Sprintf("%v:%v-%v", span.Chr, span.Start, span.End))
		}
	}
	cmd := exec.Command(args[0], args[1:]...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if e := cmd.Run(); e != nil {
		return fmt.Errorf("HaplotypeCall: %w", e)
	}

	return CreateDone(gvcfpath)
}

var commentRe = regexp.MustCompile(`^#`)

func CatGvcf(w io.Writer, path string, writeHeader bool) (err error) {
	r, e := zfile.Open(path)
	if e != nil {
		return e
	}
	defer func() {
		err = errors.Join(err, r.Close())
	}()

	if writeHeader {
		_, e := io.Copy(w, r)
		return e
	}

	s := bufio.NewScanner(r)
	s.Buffer([]byte{}, 1e15)
	for s.Scan() {
		if commentRe.MatchString(s.Text()) {
			continue
		}
		if _, e := fmt.Fprintln(w, s.Text()); e != nil {
			return e
		}
	}
	return s.Err()
}

func Tabix(vcfpath string) error {
	if IsDone(vcfpath) {
		return nil
	}
	cmd := exec.Command("tabix", "-p", "vcf", vcfpath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func JoinGvcfs(gvcfpre string, gogogo bool, threads int, chrs ...Span) (err error) {
	outpath := gvcfpre + ".g.vcf.gz"
	if !gogogo && IsDone(outpath) {
		return nil
	}
	out, e := zfile.CreateBgz(outpath, threads)
	if e != nil {
		return e
	}
	defer func() {
		err = errors.Join(err, out.Close())
	}()

	writeHeader := true
	for _, chr := range chrs {
		gvcfpath := gvcfpre + "_chr_" + chr.Chr + ".g.vcf.gz"
		if e := CatGvcf(out, gvcfpath, writeHeader); e != nil {
			return e
		}
		writeHeader = false
	}
	return nil
}

func ParseChrs(path string) (spans []Span, err error) {
	r, e := os.Open(path)
	if e != nil {
		return nil, e
	}
	defer func() {
		err = errors.Join(r.Close())
	}()

	s := bufio.NewScanner(r)
	for s.Scan() {
		spans = append(spans, Span{Chr: s.Text(), FullChr: true})
	}
	return spans, s.Err()
}

func HaplotypeCallSplit(fapath, bampath, gvcfpre, name, chrspath string, memoryGb int, gogogo bool, threads int, rq *RunQueue) error {
	chrs, e := ParseChrs(chrspath)
	if e != nil {
		return e
	}
	var g errgroup.Group
	for _, chr := range chrs {
		g.Go(func() error {
			return rq.RunErr(func() error {
				return HaplotypeCallSpans(fapath, bampath, gvcfpre + "_chr_" + chr.Chr + ".g.vcf.gz", name, memoryGb, gogogo, chr)
			})
		})
	}
	if e := g.Wait(); e != nil {
		return e
	}
	e = rq.RunErr(func() error {
		if e := JoinGvcfs(gvcfpre, gogogo, threads, chrs...); e != nil {
			return e
		}
		if e := Tabix(gvcfpre + ".g.vcf.gz"); e != nil {
			return e
		}
		return nil
	})
	if e != nil {
		return e
	}
	return CreateDone(gvcfpre + ".g.vcf.gz")
}

func DeleteGvcfs(f Flags, s ReadSet) error {
	gvcfpath := f.Outpre + "_" + s.Name + ".g.vcf.gz"
	if FileExists(gvcfpath) {
		if e := os.Remove(gvcfpath); e != nil {
			return e
		}
	}
	donepath := gvcfpath + ".done"
	if FileExists(donepath) {
		if e := os.Remove(donepath); e != nil {
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

func IsDone(path string) bool {
	return FileExists(path + ".done")
}

func CreateDone(path string) error {
	_, e := os.Create(path + ".done")
	return e
}

func CombineGvcf(fapath, gvcfoutpath, vcfoutpath string, memoryGb int, gogogo bool, gvcfpaths ...string) error {
	memstr := fmt.Sprintf("-Xmx%vg", memoryGb)

	if gogogo || !IsDone(gvcfoutpath) {
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
		if e := CreateDone(gvcfoutpath); e != nil {
			return e
		}
	}

	if gogogo || !IsDone(vcfoutpath) {
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
		if e := CreateDone(vcfoutpath); e != nil {
			return e
		}
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
	if e := BwaIndex(f.RefPath, f.Gogogo); e != nil {
		return e
	}
	if e := Faidx(f.RefPath, f.Gogogo); e != nil {
		return e
	}
	if e := CreateDict(f.PicardCmd, f.RefPath, f.Gogogo); e != nil {
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

	rq := NewRunQueue(f.Nproc)

	var group errgroup.Group
	var gvcfpaths []string
	for _, set := range sets {
		set := set

		bampath := f.Outpre + "_" + set.Name + ".bam"
		if f.NoAln {
			bampath = set.BamPath
		}

		bampathrg := f.Outpre + "_" + set.Name + "_rg.bam"
		gvcfpath := f.Outpre + "_" + set.Name + ".g.vcf.gz"
		gvcfpre := f.Outpre + "_" + set.Name
		gvcfpaths = append(gvcfpaths, gvcfpath)
		fwd := set.ForwardPath
		rev := set.ReversePath

		group.Go(func() error {
			if !f.NoAln {
				if f.Trim {
					fwd = f.Outpre + "_" + set.Name + "_forward_trimmed.fq.gz"
					rev = f.Outpre + "_" + set.Name + "_reverse_trimmed.fq.gz"

					e := rq.RunErr(func() error {
						return Trim(set.ForwardPath, set.ReversePath, f.Outpre + "_" + set.Name, f.Threads, f.Gogogo)
					})
					if e != nil {
						return e
					}
					if f.DeleteTempFiles {
						defer func() {
							errors.Join(err, DeleteTrim(f, set))
						}()
					}
					if e != nil {
						return e
					}
				}
				e := rq.RunErr(func() error {
					return BwaMem(f.RefPath, fwd, rev, bampath, f.Threads, f.Gogogo)
				})
				if e != nil {
					return e
				}
				if f.DeleteTempFiles {
					defer func() {
						errors.Join(err, DeleteBam(f, set))
					}()
				}
			}

			if !f.NoHaploCall {
				e := rq.RunErr(func() error {
					return AddRG(bampath, bampathrg, set.Name, f.Gogogo)
				})
				if e != nil {
					return e
				}
				if f.DeleteTempFiles {
					defer func() {
						errors.Join(err, DeleteRGBam(f, set))
					}()
				}
				e = rq.RunErr(func() error {
					return SamIndex(bampathrg, f.Gogogo)
				})
				if e != nil {
					return e
				}
				if f.ChrsPath != "" {
					if e := HaplotypeCallSplit(f.RefPath, bampathrg, gvcfpre, set.Name, f.ChrsPath, f.MemoryGb, f.Gogogo, f.Threads, rq); e != nil {
						return e
					}
				} else {
					e := rq.RunErr(func() error {
						return HaplotypeCall(f.RefPath, bampathrg, gvcfpath, set.Name, f.MemoryGb, f.Gogogo)
					})
					if e != nil {
						return e
					}
				}
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

	if !f.NoCombine {
		goutpath := f.Outpre + ".g.vcf.gz"
		if f.DeleteTempFiles {
			defer func() {
				errors.Join(err, DeletePath(goutpath))
			}()
		}
		outpath := f.Outpre + ".vcf.gz"
		return CombineGvcf(f.RefPath, goutpath, outpath, f.SerialMemoryGb, f.Gogogo, gvcfpaths...)
	}

	return nil
}
