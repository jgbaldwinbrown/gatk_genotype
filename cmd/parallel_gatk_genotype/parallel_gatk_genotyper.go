package main

import (
	"log"

	"github.com/jgbaldwinbrown/gatk_genotype/pkg"
)

func main() {
	bamlines, e := ggtype.ReadPathLines("bampaths.txt")
	if e != nil {
		log.Fatal(e)
	}

	for _, line := range bamlines {
		if e := ggtype.MakeProcessDir(line); e != nil {
			log.Fatal(e)
		}
	}

	if e := ggtype.MakeRunScript(bamlines); e != nil {
		log.Fatal(e)
	}
}
