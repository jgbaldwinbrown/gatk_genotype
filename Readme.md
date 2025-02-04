# gatk\_genotype

A thin Golang wrapper around GATK for genotyping

## Introduction

It can be challenging to set up a pipeline for running BWA and GATK for multiple-sample genotyping. This short Go program does exactly that.

## Prerequisites

- GATK 3 or 4
- picard tools
- samtools
- BWA
- The Golang build suite

## Installation

- In this directory, just run `go build gatk_genotype.go`, then copy `gatk_genotype` into a directory listed in your PATH.

## Use

```
Usage of gatk_genotype:
  -m int
    	Memory to use (integer, gigabytes) (default 8)
  -o string
    	Output prefix (default "out")
  -r string
    	Path to reference .fa file (required)
  -s string
    	Path to tab-separated table containing pairs of forward and reverse read paths, one line per sample (required). Format: name (tab) forward.fq.gz (tab) reverse.fq.gz
  -t int
    	Threads to use (default 1)
```
