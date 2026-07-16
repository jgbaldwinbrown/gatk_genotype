package main

import (
	"flag"
	"log"

	"github.com/jgbaldwinbrown/gatk_genotype/pkg"
)

func main() {
	var f ggtype.Flags
	flag.StringVar(&f.RefPath, "r", "", "Path to reference .fa file (required)")
	flag.StringVar(&f.SeqPairsPath, "s", "", "Path to tab-separated table containing pairs of forward and reverse read paths, one line per sample (required). Format: name (tab) forward.fq.gz (tab) reverse.fq.gz")
	flag.StringVar(&f.Outpre, "o", "out", "Output prefix")
	flag.IntVar(&f.Threads, "t", 1, "Threads to use")
	flag.IntVar(&f.MemoryGb, "m", 8, "Memory to use for parallel processes (integer, gigabytes).")
	flag.IntVar(&f.SerialMemoryGb, "M", 8, "Memory to use for serial processes (integer, gigabytes).")
	flag.IntVar(&f.Nproc, "n", 1, "Number of simultaneous runs of BWA / picard to run")
	flag.BoolVar(&f.Trim, "T", false, "Also trim input files with trimmomatic")
	flag.BoolVar(&f.NoAln, "noaln", false, "Skip aligning (it is already done). Not compatible with -s.")
	flag.BoolVar(&f.NoCombine, "nocombine", false, "Skip combining gvcf files (will do this later).")
	flag.BoolVar(&f.DeleteTempFiles, "d", false, "Delete all intermediate files except for index files and the final .vcf.gz file.")
	flag.StringVar(&f.BamPathsPath, "bams", "", "Path to file containing paths to bams, one line per sample. Format: name (tab) path.bam. Not compatible with -s and requires -noaln.")
	flag.StringVar(&f.PicardCmd, "picard", "picard-tools", "Command to invoke picard tools.")
	flag.StringVar(&f.ChrsPath, "chrs", "", "Path to a newline-separated list of chrs to call haplotypes in; using this allows haplotypes to be called in pieces for cleaner restarts, one piece per chromosome.")
	flag.BoolVar(&f.Gogogo, "g", false, "Go go go: rebuild everything from scratch, ignoring \".done\" files.")
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

	if e := ggtype.FullFQFMimic(f); e != nil {
		log.Fatal(e)
	}
}
