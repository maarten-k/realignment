
rule realign:
    input:
        cram=lambda wildcards: samplesmap[wildcards.SM],
    output:
        "{SM}.sorted.cram",
    threads: 8
    conda:
        "../envs/original.yaml"
    shadow:
        "copy-minimal"
    params:
        ref=config["fasta"],
        refold=config["fasta_input"],
    shell:
        """
         echo "start $(date)"
        samtools index -@ 8 {input.cram}
        echo "index $(date)"
        samtools bamshuf -@ 8 --reference {params.refold} --output-fmt BAM -uOn 128 {input.cram} {wildcards.SM}.tmp |samtools bam2fq -@ 8 -t -s /dev/null -1 {wildcards.SM}.R1.fq.gz -2 {wildcards.SM}.R2.fq.gz - > /dev/null
        echo "bamshuf $(date)"
        bwa mem -K 100000000 -t 8 -Y {params.ref} -R "@RG\\tID:{wildcards.SM}\\tLB:{wildcards.SM}\\tSM:{wildcards.SM}\\tPL:ILLUMINA" {wildcards.SM}.R1.fq.gz {wildcards.SM}.R2.fq.gz 2>> {wildcards.SM}.log | samblaster -a --addMateTags | samtools view -h --threads {threads} -bS > {wildcards.SM}.aln.bam
        echo "bwa $(date)"
        rm -f {wildcards.SM}.R1.fq.gz {wildcards.SM}.R2.fq.gz
        samtools sort -@ 8 {wildcards.SM}.aln.bam > {wildcards.SM}.sorted.bam
        echo "sorted $(date)"
        rm -f {wildcards.SM}.aln.bam
        samtools view -@ 8 -h -T {params.ref} -C {wildcards.SM}.sorted.bam > {output}
        rm -f {wildcards.SM}.sorted.bam
        echo "cramed $(date)"
        """


# BQSR_Loci=$(echo -L chr{1..22} | sed 's/ / -L /g' | sed 's/-L -L/-L/g')


rule recalibrate:
    input:
        cram="{SM}.sorted.cram",
    output:
        cram="{SM}.final-gatk.cram",
        index="{SM}.final-gatk.cram.crai",
        metrics_dedup="{SM}.dedupMetrics.txt",
    log:
        dedup="{SM}.dedup.log",
        bqsr="{SM}.bqsr.log",
    threads: 8
    conda:
        "../envs/original.yaml"
    shadow:
        "copy-minimal"
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
        bqsr_loci="-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22",
        Mills=config["Mills"],
        KnownIndels=config["KnownIndels"],
    shell:
        """
        wrk="."
        echo "start $(date)"
        samtools view -h -@ 8 -T {params.ref} {input.cram} -b > {wildcards.SM}.bam
         echo "cram2bam input $(date)"
        samtools index -@ 8 {wildcards.SM}.bam
        echo "indexed input $(date)"
        picard -Djava.io.tmpdir=${{wrk}} MarkDuplicates I={wildcards.SM}.bam AS=true O={wildcards.SM}.dedup.bam METRICS_FILE={output.metrics_dedup} QUIET=true COMPRESSION_LEVEL=0 2>> {log.dedup}
        rm {wildcards.SM}.bam
        echo "MarkDuplicated $(date)"
        samtools sort -@ 8 -o {wildcards.SM}.dedup-sorted.bam {wildcards.SM}.dedup.bam
        echo "sorted dedup.bam $(date)"
        rm {wildcards.SM}.dedup.bam
        echo "removed dedup.bam $(date)"

        samtools index -@ 8 {wildcards.SM}.dedup-sorted.bam
        echo "indexed dedup-sorted.bam $(date)"

        gatk3 -Xmx8g -Djava.io.tmpdir=${{wrk}} -T BaseRecalibrator -I {wildcards.SM}.dedup-sorted.bam -R {params.ref} -o {wildcards.SM}.recal -nct {threads} --downsample_to_fraction .1 {params.bqsr_loci} -knownSites {params.dbSNP} -knownSites {params.Mills} -knownSites {params.KnownIndels} 2>> {log.bqsr}
        echo "BaseRecalibrator done $(date)"
        gatk3 -Xmx8g -Djava.io.tmpdir=${{wrk}}  -T PrintReads -I {wildcards.SM}.dedup-sorted.bam -R {params.ref} -nct {threads} --BQSR {wildcards.SM}.recal -o {wildcards.SM}.bqsr.bam --globalQScorePrior -1.0 --preserve_qscores_less_than 6 --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30 --disable_indel_quals 2>> {log.bqsr}
        echo "printedreads $(date)"
        rm {wildcards.SM}.dedup-sorted.bam
         echo "removed dedup-sorted $(date)"
        samtools sort -@ 8 {wildcards.SM}.bqsr.bam -o {wildcards.SM}.final.bam
        echo "sorted $(date)"
        rm {wildcards.SM}.bqsr.bam
        samtools view -h -@ 8 -T {params.ref} -C {wildcards.SM}.final.bam > {output.cram}
        echo "bam2cram output$(date)"
        samtools index -@ 8 {output.cram}
        echo "cramindex output$(date)"
        """


rule HaplotyperWGS:
    input:
        "{SM}.final-gatk.cram",
    output:
        gvcf="{SM}.g.vcf.gz",
        index="{SM}.g.vcf.gz.tbi",
    log:
        "{SM}.vcf.log",
    conda:
        "../envs/original.yaml"
    threads: 2
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
    shadow:
        "copy-minimal"
    shell:
        """
        wrk=$(pwd)
        # Call variants chr1
        echo -e "\\n\\nCalling chr1\\n" >> {log}
        mkdir -p ${{wrk}}/chr1
        cd ${{wrk}}/chr1
        set -x
        gatk --java-options -Djava.io.tmpdir=. HaplotypeCaller -R {params.ref}  --dbsnp {params.dbSNP}  -I ../{input} -O ${{wrk}}/chr1/{wildcards.SM}.chr1.g.vcf.gz -L chr1 -ERC GVCF --native-pair-hmm-threads {threads} &>> {log}


        # Initialize gVCF and remove temp chrom dir
        if [ -f ${{wrk}}/chr1/{wildcards.SM}.chr1.g.vcf.gz ]
        then
            zcat ${{wrk}}/chr1/{wildcards.SM}.chr1.g.vcf.gz > ${{wrk}}/{wildcards.SM}.g.vcf
            cd ${{wrk}}
            rm -fr chr1
        else
            echo -e "\\nError during variant, exiting\\n"
            cd ..
            rm -fr  ${{wrk}}/{wildcards.SM}.g.vcf
            exit 1
        fi


        # Iteratively haplotype caller and append to gVCF
        for chrom in chr{{2..22}} chr{{X..Y}}
            do

            # Call variants
            echo -e "\\n\\nCalling ${{chrom}}\\n" >> {log}
            mkdir -p ${{wrk}}/${{chrom}}
            cd ${{wrk}}/${{chrom}}
            gatk --java-options -Djava.io.tmpdir=.  HaplotypeCaller -R {params.ref} --dbsnp {params.dbSNP} -I ../{input} -O ${{wrk}}/${{chrom}}/{wildcards.SM}.${{chrom}}.g.vcf.gz -L ${{chrom}} -ERC GVCF --native-pair-hmm-threads 2 &>> ${{wrk}}/{wildcards.SM}.vcf.log

            # Append to gVCF
            zcat ${{wrk}}/${{chrom}}/{wildcards.SM}.${{chrom}}.g.vcf.gz | grep -v "#" >> ${{wrk}}/{wildcards.SM}.g.vcf
            cd ${{wrk}}
            rm -fr ${{chrom}}
        done


        bgzip -c ${{wrk}}/{wildcards.SM}.g.vcf > {output.gvcf}
        tabix -p vcf -f {output.gvcf}


        """
