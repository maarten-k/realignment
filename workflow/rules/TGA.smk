configfile: "config/config.yaml"

from collections import defaultdict

nested_dict = lambda: defaultdict(nested_dict)
nest = nested_dict()
z = [l.strip() for l in open("TGA_fastq.txt").readlines()]
for x in z:
    foundg = re.search("([^\/]+ONCOR-\d+)_(\w+_L\d{3})_(R[12])", x).groups()
    print(f"{foundg[0]} {foundg[1]} {foundg[2]}\t {x}")
    nest[foundg[0]][foundg[1]][foundg[2]] = x


print(nest)

rule test_fastq:
    input:
        "FR26962644_TR-ALS2074_R1_merged.fastq.gz",
    message:
        "Hi!"
    shell:
        "echo hi"


rule realign_pipe_all:
    input:
        expand(
            "results/{SM}/{SM}.sorted_pipe.bam", SM=nest.keys()
        ),	


rule fastp:
    input:
        R1=lambda wildcards: nest[wildcards.sample][wildcards.lane]["R1"],
        R2=lambda wildcards: nest[wildcards.sample][wildcards.lane]["R2"],
    output:
        pipe("{sample}__{lane}.fastq"),
    threads:0.1
    log:
        html="results/{sample}/fastp/{sample}_{lane}.html",
        json="results/{sample}/fastp/{sample}_{lane}.json",
    conda:
        "../envs/fastp.yaml"
    resources:
        mem_mb=2000,
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -j {log.json} -h {log.html} --stdout > {output}
        """


rule realign_pipe:
    input:
        lambda wildcards: expand(
            "{{SM}}__{lane}.fastq", lane=nest[wildcards.SM].keys()
        ),
    output:
        bam="results/{SM}/{SM}.sorted_pipe.bam",
        index="results/{SM}/{SM}.sorted_pipe.bam.bai",
    threads: 28
    resources:
        mem_mb=46000,
        disk_mb=1000,
        runtime=1200
    conda:
        "../envs/realign.yaml"
    log:
        "results/{SM}/log/{SM}_bwa.log",
    params:
        ref=config["fasta"],
        refold=config["fasta_input"],
    benchmark:
        "benchmarks/{SM}.bwa2_pipe.benchmark.txt"
    shell:
        """
       mkdir -p tmp/{wildcards.SM}
       cat {input} |\
       bwa-mem2 mem -K 100000000 -pt {threads} -Y {params.ref} -R "@RG\\tID:{wildcards.SM}\\tLB:{wildcards.SM}\\tSM:{wildcards.SM}\\tPL:ILLUMINA" - 2>> {log} | \
       samblaster -a --addMateTags |\
       samtools sort -T tmp/{wildcards.SM}/{wildcards.SM} --write-index -m 4G -l 1 -@ 4  -o {output.bam}##idx##{output.index}
       rm -rf tmp/{wildcards.SM}
        """

