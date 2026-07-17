package ggtype

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

const rawtemplate = `#!/bin/bash
#SBATCH -t 72:00:00    #max:    72 hours (24 on ash)
#SBATCH -N 1          #format: count or min-max
#SBATCH -A owner-guest    #values: yandell, yandell-em (ember), ucgd-kp (kingspeak)
#SBATCH -p redwood-guest    #kingspeak, ucgd-kp, kingspeak-freecycle, kingspeak-guest
#SBATCH -J %v        #Job name
#SBATCH --ntasks=28
#SBATCH --mem=128G

NUM_CORES="${SLURM_CPUS_ON_NODE}"
NUM_NODE_MEM="${SLURM_MEM_PER_NODE}"

echo $NUM_NODE_MEM
NUM_NODE_MEM_PER_CPU=%vecho "${NUM_NODE_MEM} / ${NUM_CORES}" | bc%v
echo $NUM_NODE_MEM_PER_CPU
NUM_NODE_MEM_GB=%vecho "${NUM_NODE_MEM} / 1000" | bc%v
echo $NUM_NODE_MEM_GB
NUM_NODE_MEM_PER_CPU_GB=%vecho "${NUM_NODE_MEM_PER_CPU} / 1000" | bc%v
echo $NUM_NODE_MEM_PER_CPU_GB

module load bwa
module load samtools
module load picard
module load gatk
module load golang
module load htslib

mkdir -p /scratch/general/pe-nfs1/u6012238/yang_run5
# rsync -avP ./ref/ /scratch/general/pe-nfs1/u6012238/yang_run5/
# touch /scratch/general/pe-nfs1/u6012238/yang_run5/human_g1k_v37_decoy.fasta.bwt.done
# touch /scratch/general/pe-nfs1/u6012238/yang_run5/human_g1k_v37_decoy.dict.done

gatk_genotype \
	-m "${NUM_NODE_MEM_PER_CPU_GB}" \
	-M "${NUM_NODE_MEM_GB}" \
	-t "${NUM_CORES}" \
	-o /scratch/general/pe-nfs1/u6012238/yang_run5/out \
	-r /scratch/general/pe-nfs1/u6012238/yang_run5/human_g1k_v37_decoy.fasta \
	-noaln \
	-nocombine \
	-bams %v \
	-n "${NUM_CORES}" \
	-chrs chrs_reforder.txt \
	-picard "java -jar /uufs/chpc.utah.edu/sys/installdir/r8/picard/3.3.0/picard.jar"
`

var template = fmt.Sprintf(rawtemplate, "%v", "`", "`", "`", "`", "`", "`", "%v")

func ReadPathLines(path string) (lines []string, err error) {
	r, e := os.Open(path)
	if e != nil {
		return nil, e
	}
	defer func() { err = errors.Join(err, r.Close()) }()

	s := bufio.NewScanner(r)
	s.Buffer([]byte{}, 1e15)
	for s.Scan() {
		lines = append(lines, s.Text())
	}
	return lines, s.Err()
}

func WriteToPath(path, content string) (err error) {
	w, e := os.Create(path)
	if e != nil {
		return e
	}
	defer func() { err = errors.Join(err, w.Close()) }()

	_, e = fmt.Fprint(w, content)
	return e
}

func WriteScript(scriptpath, line string) (err error) {
	name, _, ok := strings.Cut(line, "\t")
	if !ok {
		return fmt.Errorf("WriteScript: no tab in line\"%v\"", line)
	}
	bampathout := filepath.Join("scripts", name + "_dir", name + "_bampath.txt")
	script := fmt.Sprintf(template, name, bampathout)
	return WriteToPath(scriptpath, script)
}

func MakeProcessDir(line string) error {
	name, _, ok := strings.Cut(line, "\t")
	if !ok {
		return fmt.Errorf("MakeProcessDir: no tab in line\"%v\"", line)
	}
	dirpath := filepath.Join("scripts", name + "_dir")
	if e := os.MkdirAll(dirpath, 0755); e != nil {
		return e
	}
	bampathout := filepath.Join(dirpath, name + "_bampath.txt")
	if e := WriteToPath(bampathout, line); e != nil {
		return e
	}
	scriptpath := filepath.Join(dirpath, name + "_run.sh")
	if e := WriteScript(scriptpath, line); e != nil {
		return e
	}

	return nil
}

func MakeRunScript(bamlines []string) error {
	var b strings.Builder
	fmt.Fprintf(&b, `#!/bin/sh
set -e

`)

	for _, line := range bamlines {
		name, _, ok := strings.Cut(line, "\t")
		if !ok {
			return fmt.Errorf("MakeRunScript: no tab in line\"%v\"", line)
		}
		scriptpath := filepath.Join("scripts", name + "_dir", name + "_run.sh")
		fmt.Fprintf(&b, "sbatch %v\n", scriptpath)
	}
	return WriteToPath("runall.sh", b.String())
}
