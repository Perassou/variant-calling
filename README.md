# variant-calling

## Purpose

This workflow takes short read data (fastq files) as input and outputs variant calls (a vcf file).
The pipeline is as follow: bwa mem → samtools view/sort→ freebayes.


## How to build the image
The user will precise the reference file when building the docker image that has a `bwa` aligner index.

```
docker build -t bwa_alignment --build-arg genomes_dir=alignment/MN908947.3.fasta -f variant_calling_pipeline.Dockerfile .
```

## How to run 

```
docker run -v $(pwd):/data -t bwa_alignment -1 sample.R1.paired.fq -2 sample.R2.paired.fq
```

**Typical runtime for the example provided: 60 seconds**