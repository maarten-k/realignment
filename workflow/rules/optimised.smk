chromosomes = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
]


wildcard_constraints:
    chr="chr[0-9XYM]{1,2}",
    SM="[A-Za-z0-9_-]+",


localrules:
    fetch_from_gridstorage,


ruleorder: index_cram > convert2cram_with_oldsamtools

rule realign_pipe:
    input:
        cram="ingress/{SM}.cram",
    output:
        bam=temp("results/{SM}/{SM}.sorted_pipe.bam"),
        index=temp("results/{SM}/{SM}.sorted_pipe.bam.bai"),
    threads: 8
    resources:
        mem_mb=54000,
    conda:
        "../envs/realign.yaml"
    log:
        "{SM}_bwa.log",
    shadow:
        "copy-minimal"
    params:
        ref=config["fasta"],
        refold=config["fasta_input"],
    benchmark:
        "benchmarks/{SM}.bwa2_pipe.benchmark.txt"
    shell:
        """

       samtools collate -@{threads} --output-fmt BAM -uOn 128 {input.cram} /tmp/{wildcards.SM}.tmp |\
       samtools bam2fq -@8 -t -s /dev/null |\
       bwa-mem2 mem -K 100000000 -pt {threads} -Y {params.ref} -R "@RG\\tID:{wildcards.SM}\\tLB:{wildcards.SM}\\tSM:{wildcards.SM}\\tPL:ILLUMINA" - 2>> {log} | \
       samblaster -a --addMateTags |\
       samtools sort -T /tmp/{wildcards.SM} --write-index  -m 2G -l 1 -@ {threads}  -o {output.bam}##idx##{output.index}

        """

# BQSR_Loci=$(echo -L chr{1..22} | sed 's/ / -L /g' | sed 's/-L -L/-L/g')


rule recalibrate_new:
    input:
        bam="results/{SM}/{SM}.sorted_pipe.bam",
        index="results/{SM}/{SM}.sorted_pipe.bam.bai",
    output:
        bam=temp("results/{SM}/{SM}.bqsr.bam"),
        metrics_dedup="results/{SM}/{SM}.dedupMetrics.txt",
    log:
        dedup="results/{SM}/log/{SM}.dedup.log",
        bqsr="results/{SM}/log/{SM}.bqsr.log",
    resources:
        mem_mb=8000,
        partition="fat"
    threads: 1
    priority: 50
    conda:
        "../envs/optimised.yaml"
    shadow:
        "shallow"
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
        bqsr_loci="-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22",
        Mills=config["Mills"],
        KnownIndels=config["KnownIndels"],
    benchmark:
        "benchmarks/recalibrate_new-{SM}.txt"
    shell:
        """
        echo "start $(date)"
        picard -Djava.io.tmpdir={resources.tmpdir} MarkDuplicates I={input.bam} AS=true O={wildcards.SM}.dedup.bam METRICS_FILE={output.metrics_dedup} QUIET=true COMPRESSION_LEVEL=1 2>> {log.dedup}
        echo "MarkDuplicated $(date)"
        samtools index -@ {threads} {wildcards.SM}.dedup.bam
        echo "indexed dedup $(date)"
        gatk3 -Xmx8G -Djava.io.tmpdir={resources.tmpdir} -T BaseRecalibrator -I {wildcards.SM}.dedup.bam -R {params.ref} -o {wildcards.SM}.recal -nct {threads} --downsample_to_fraction .1 {params.bqsr_loci} -knownSites {params.dbSNP} -knownSites {params.Mills} -knownSites {params.KnownIndels} 2>> {log.bqsr}
        echo "BaseRecalibrator done $(date)"

        gatk3 -Xmx8G -Djava.io.tmpdir={resources.tmpdir} -T PrintReads -I {wildcards.SM}.dedup.bam -R {params.ref} -nct {threads} --BQSR {wildcards.SM}.recal -o {output.bam} --globalQScorePrior -1.0 --preserve_qscores_less_than 6 --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30 --disable_indel_quals 2>> {log.bqsr}
        echo "printedreads $(date)"
        """


rule convert2cram_with_oldsamtools:
    input:
        bam="results/{SM}/{SM}.bqsr.bam",
    output:
        cram="results/{SM}/{SM}.final-gatk.cram",
        index="results/{SM}/{SM}.final-gatk.cram.crai",
    resources:
        mem_mb=2000,
    threads: 1
    priority: 60
    conda:
        "../envs/samtools1.9.yaml"
    params:
        ref=config["fasta"],
    benchmark:
        "benchmarks/convert2cram_with_oldsamtools-{SM}.txt"
    shell:
        """
        samtools view -C -O CRAM -T {params.ref} -@ {threads} {input.bam} -o {output.cram}
        samtools index -@ {threads} {output.cram}

        """


rule haplotype_per_chr:
    input:
        cram="results/{SM}/{SM}.final-gatk.cram",
        index="results/{SM}/{SM}.final-gatk.cram.crai",
    output:
        temp("results/{SM}/gvf/{SM}-{chr}.g.vcf.gz"),
    conda:
        "../envs/GATK_haplotyper.yaml"
    threads: 1
    resources:
        mem_mb=8000,
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
    benchmark:
        "benchmarks/haployper{SM}-{chr}.txt"
    shell:
        """
        gatk --java-options -Djava.io.tmpdir={resources.tmpdir} HaplotypeCaller -R {params.ref}  --dbsnp {params.dbSNP}  -I {input.cram} -O {output} -L {wildcards.chr} --seconds-between-progress-updates 100 -ERC GVCF --native-pair-hmm-threads {threads} 
        """


rule merge_gvcf:
    input:
        expand("{{SM}}/gvf/{{SM}}-{chr}.g.vcf.gz", chr=chromosomes),
    output:
        gvcf="{SM}.g.vcf.gz",
        index="{SM}.g.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=8000,
    benchmark:
        "benchmarks/merge_gvcf-{SM}.txt"
    shadow:
        "copy-minimal"
    threads: 1
    priority: 100
    shell:
        """
        bcftools concat -O z8  -o {output.gvcf} {input}
        tabix -p vcf  {output.gvcf}
        """


rule HaplotyperExome:
    input:
        cram="{SM}.final-gatk.cram",
        index="{SM}.final-gatk.cram.crai",
    output:
        gvcf="{SM}.WXS.g.vcf.gz",
        index="{SM}.WXS.g.vcf.gz.tbi",
    log:
        "{SM}.WXS.vcf.log",
    conda:
        "../envs/optimised.yaml"
    threads: 2
    resources:
        mem_mb=8000,
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
        tgt=config["tgt"],
    shadow:
        "copy-minimal"
    shell:
        """
        touch {input.index}
        gatk --java-options -Djava.io.tmpdir={resources.tmpdir} HaplotypeCaller -R {params.ref}  --dbsnp {params.dbSNP}  -I {input.cram} -O {output.gvcf} -ERC GVCF -L {params.tgt} --seconds-between-progress-updates 100 --native-pair-hmm-threads {threads} &>> {log}
        mv {output.gvcf} {output.gvcf}.tmp
        zcat {output.gvcf}.tmp |bgzip -@ {threads} -l 6 -c> {output.gvcf}

        tabix -fp vcf {output.gvcf}
        """


rule index_cram:
    input:
        cram="{SM}.final-gatk.cram",
    output:
        index="{SM}.final-gatk.cram.crai",
    resources:
        mem_mb=1000,
    conda:
        "../envs/optimised.yaml"
    threads: 1
    shell:
        """
        samtools index {input.cram}
        """
