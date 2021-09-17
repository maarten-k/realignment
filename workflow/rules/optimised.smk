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
    SM="[A-Za-z0-9]+",


rule realign:
    input:
        cram=lambda wildcards: samplesmap[wildcards.SM],
    output:
        bam="{SM}.sorted.bam",
        index="{SM}.sorted.bam.bai",
    threads: 8
    conda:
        "../envs/optimised.yaml"
    shadow:
        "copy-minimal"
    params:
        ref=config["fasta"],
        refold=config["fasta_input"],
    benchmark:
        "benchmarks/{SM}.bwa2.benchmark.txt"
    shell:
        """
        echo "start $(date)"
        samtools index -@ {threads} {input.cram}
        echo "indexed $(date)"
        samtools collate -@ {threads} --reference {params.refold} --output-fmt BAM -uOn 128 {input.cram} {wildcards.SM}.tmp |samtools bam2fq -@ {threads} -t -s /dev/null -1 {wildcards.SM}.R1.fq.gz -2 {wildcards.SM}.R2.fq.gz - > /dev/null
         echo "collated $(date)"
        bwa-mem2 mem -K 100000000 -t {threads} -Y {params.ref} -R "@RG\\tID:{wildcards.SM}\\tLB:{wildcards.SM}\\tSM:{wildcards.SM}\\tPL:ILLUMINA" {wildcards.SM}.R1.fq.gz {wildcards.SM}.R2.fq.gz 2>> {wildcards.SM}.log | samblaster -a --addMateTags | samtools view -h1 --threads {threads} -bS > {wildcards.SM}.aln.bam
         echo "bwa $(date)"
        rm -f {wildcards.SM}.R1.fq.gz {wildcards.SM}.R2.fq.gz
        samtools sort  -O bam -l 1 --write-index -@ {threads} {wildcards.SM}.aln.bam -o {output.bam}##idx##{output.index}
         echo "sorted $(date)"
        """


# BQSR_Loci=$(echo -L chr{1..22} | sed 's/ / -L /g' | sed 's/-L -L/-L/g')


rule recalibrate_new:
    input:
        bam="{SM}.sorted.bam",
        bamindex="{SM}.sorted.bam.bai",
    output:
        bam="{SM}.bqsr.bam",
        metrics_dedup="{SM}.dedupMetrics.txt",
    log:
        dedup="{SM}.dedup.log",
        bqsr="{SM}.bqsr.log",
    threads: 1
    conda:
        "../envs/optimised.yaml"
    shadow:
        "copy-minimal"
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
        picard -Djava.io.tmpdir={resources.tmpdir} MarkDuplicates I={input.bam} AS=true O={wildcards.SM}.dedup.bam METRICS_FILE={output.metrics_dedup} QUIET=true COMPRESSION_LEVEL=0 2>> {log.dedup}
        #why sort this stuff
        echo "MarkDuplicated $(date)"
        samtools sort -m 6G -O bam -l 1 --write-index  -@ {threads} {wildcards.SM}.dedup.bam -o {wildcards.SM}.dedup-sorted.bam##idx##{wildcards.SM}.dedup-sorted.bam.bai 
        echo "sorted dedup.bam $(date)"
        rm {wildcards.SM}.dedup.bam


        gatk3 -Djava.io.tmpdir={resources.tmpdir} -T BaseRecalibrator -I {wildcards.SM}.dedup-sorted.bam -R {params.ref} -o {wildcards.SM}.recal -nct {threads} --downsample_to_fraction .1 {params.bqsr_loci} -knownSites {params.dbSNP} -knownSites {params.Mills} -knownSites {params.KnownIndels} 2>> {log.bqsr}
        echo "BaseRecalibrator done $(date)"

        gatk3 -Djava.io.tmpdir={resources.tmpdir} -T PrintReads -I {wildcards.SM}.dedup-sorted.bam -R {params.ref} -nct {threads} --BQSR {wildcards.SM}.recal -o {output.bam} --globalQScorePrior -1.0 --preserve_qscores_less_than 6 --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30 --disable_indel_quals 2>> 2{log.bqsr}
        rm {wildcards.SM}.dedup-sorted.bam
        echo "printedreads $(date)"
        """


rule convert2cram_with_oldsamtools:
    input:
        bam="{SM}.bqsr.bam",
    output:
        cram="{SM}.final-gatk.cram",
        index="{SM}.final-gatk.cram.crai",
    shadow:
        "copy-minimal"
    threads: 1
    conda:
        "../envs/samtools1.9.yaml"
    params:
        ref=config["fasta"],
    benchmark:
        "benchmarks/convert2cram_with_oldsamtools-{SM}.txt"
    shell:
        """
        samtools view -C -O CRAM -T {params.ref} -@ {threads} {input.bam} -o {output.cram}
        samtools index {output.cram}

        """


rule haplotype_per_chr:
    input:
        cram="{SM}.final-gatk.cram",
        index="{SM}.final-gatk.cram.crai",
    output:
        "{SM}/gvf/{SM}-{chr}.g.vcf.gz",
    shadow:
        "copy-minimal"
    conda:
        "../envs/optimised.yaml"
    threads: 1
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
        "../envs/optimised.yaml"
    benchmark:
        "benchmarks/merge_gvcf-{SM}.txt"
    shadow:
        "copy-minimal"
    threads: 1
    shell:
        """
        bcftools concat -O z  -o {output.gvcf} {input}
        tabix -p vcf  {output.gvcf}
        """
