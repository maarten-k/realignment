from collections import defaultdict

nested_dict = lambda: defaultdict(nested_dict)
nest = nested_dict()
z = [l.strip() for l in open("fastqlist.txt").readlines()]
for x in z:
    foundg = re.search("1534_(FR.+)/.+_(L\d{3})_(R[12])", x).groups()
    # print(f"{foundg[0]} {foundg[1]} {foundg[2]}\t {x}")
    nest[foundg[0]][foundg[1]][foundg[2]] = x


ruleorder: realign > realign_from_fastq


rule test_fastq:
    input:
        "FR26962644_TR-ALS2074_R1_merged.fastq.gz",
    message:
        "Hi!"
    shell:
        "echo hi"


rule fastp:
    input:
        R1=lambda wildcards: nest[wildcards.sample][wildcards.lane]["R1"],
        R2=lambda wildcards: nest[wildcards.sample][wildcards.lane]["R2"],
    output:
        R1="{sample}_{lane}_R1.fastq.gz",
        R2="{sample}_{lane}_R2.fastq.gz",
    conda:
        "fastp.yaml"
    resources:
        mem_mb=8000,
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}
        """


rule realign_from_fastq:
    input:
        fastq1=lambda wildcards: expand(
            "{{sample}}_{lane}_R1.fastq.gz", lane=nest[wildcards.sample].keys()
        ),
        fastq2=lambda wildcards: expand(
            "{{sample}}_{lane}_R2.fastq.gz", lane=nest[wildcards.sample].keys()
        ),
    output:
        bam=temp("{sample}.sorted.bam"),
        index=temp("{sample}.sorted.bam.bai"),
    threads: 8
    resources:
        mem_mb=24000,
    conda:
        "../envs/optimised.yaml"
    log:
        "{sample}_bwa.log",
    shadow:
        "copy-minimal"
    params:
        ref=config["fasta"],
        refold=config["fasta_input"],
    benchmark:
        "benchmarks/{sample}.bwa2.benchmark.txt"
    shell:
        """
        bwa-mem2 mem -K 100000000 -t {threads} -Y {params.ref} -R "@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"  <(cat {input.fastq1}) <(cat {input.fastq2}) 2>> {log} | samblaster -a --addMateTags| samtools sort -T /tmp/{wildcards.sample} --write-index  -m 7G -l 1 -@ {threads}  -o {output.bam}##idx##{output.index}
        """
