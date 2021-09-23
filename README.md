## About
This realignment pipeline is reimplantation the pipeline found at https://github.com/BrenKenna/data_processing. The original pipeline is tailored to run on the grid infrastructure.

To increase reproducibility, portability and maintainability the code is reimplemented in Snakemake with identical use of commands and version of programs (BWA,GATK,samtools). To check for discrepanty between the grid version and this version, a line for line comparison is done.

for the header:

`diff <(samtools view -H sample1_new.final-gatk.cram) <(samtools view -H sample1_old.final-gatk.cram )`

for the reads:

`comm -23 <(samtools view -x RG sample1_new.final-gatk.cram) <(samtools view -x RG sample1_old.final-gatk.cram )`

This is done with an Hiseq 2500 sample and a Hiseq X sample. There were no changes found, except for some paths in the header of the file. This has no consequences on further analysis.


# Optimisations:
To make use of the current software version of tools and best practices there is also an updated workflow made. The pipeline is tailored to a system with 8GB of memory per core and has a fast "scratch" space mounted at /tmp.

- replaced bwa with bwa-mem2 
- updated samtools to 1.13 (from 1.9)
- Using bam format between realignment and BPQS recalibration (saves converting from bam to cram to bam again)
- creating an index is done when creating the cram/bam file with samtools --write-index option
- schedule 1 core for deduplication instead of 8 (about 700% more efficient)
- print reads uses one core instead of 8 cores (7 times more efficient in core hours) ([Recommendations for performance optimizations when using GATK3.8 and GATK4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6842142/) see fig 1 a)
- removed sorting a bam file
- HaplotypeCaller runs now in parallel instead of serial per chromosome (reduces wall time)
- merging gvcf files is done with bcftools 

To ensure the same results, the final cram produced uses an old version of samtools which has a minor bug when creating MD tags.



# computational costs

To give a ballpark figure, timings for 1 sample is show in this table.

|                                 | New pipeline     |            | original pipeline |            |
|---------------------------------|------------------|------------|-------------------|------------|
|                                 | Wall clock hours | core hours | Wall clock hours  | core hours |
| Realign                         | 8.5              | 68.2       | 21.9              | 174.9      |
| recalibrate                     | 14               | 14         | 23.5              | 187.8      |
| convert2cram_with_oldsamtools   | 1.1              | 1.1        | NA                | NA         |
| HaplotypeWGS  | 4.5              | 53         | 22.4              | 44.8       |
| merge_gvcf                      | 0.3              | 0.3        | NA                | NA         |
| Total                           | 28.4             | **136.6**      | 67.8              | **407.5**      |


# How to run

## Prerequisites
This pipeline needs access to CVMFS repository called "softdrive" to have access to the reference files.
conda and mamba installed with Bioconda repository configured.

## Step by Step
- Clone this repo to your  system

- Install snakemake and cookiecutter
```mamba install snakemake>6.8 cookiecutter```

- install slurm profile for snakemake. see: https://github.com/Snakemake-Profiles/slurm#quickstart

- Create a tab separate list  of sample id and paths of sample:

```
sampleid1	/here/is/my/sample1.cram
sampleid4	/another/is/here/sample4.cram
```
The input files should be cram files.

- edit `config/config.yaml` and correct the "list_of_files" line to the just created tab separate list with samples and paths.
Here is also the option to use the ooriginalor optimised pipeline: to use the optimised pipeline set the  "original" option to False.

- run the pipeline with:

```
snakemake --use-conda --profile slurm -j 1000  --shadow-prefix /tmp
```
 