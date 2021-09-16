
rule realign:
    input:
        cram=lambda wildcards: samplesmap[wildcards.SM],
    output:
        bam="{SM}.sorted.bam",
        index="{SM}.sorted.bam.bai",
    threads: 8
    conda:
        "envs/optimised.yaml"
    shadow:
        "shallow-copy"
    params:
        ref=config["fasta"],
        refold=config["fasta_input"],
    benchmark:
        "benchmarks/{SM}.bwa2.benchmark.txt"
    shell:
        """
        samtools index -@ {threads} {input.cram}
        samtools collate -@ {threads} --reference {params.refold} --output-fmt BAM -uOn 128 {input.cram} {wildcards.SM}.tmp |samtools bam2fq -@ {threads} -t -s /dev/null -1 {wildcards.SM}.R1.fq.gz -2 {wildcards.SM}.R2.fq.gz - > /dev/null
        bwa-mem2 mem -K 100000000 -t {threads} -Y {params.ref} -R "@RG\\tID:{wildcards.SM}\\tLB:{wildcards.SM}\\tSM:{wildcards.SM}\\tPL:ILLUMINA" {wildcards.SM}.R1.fq.gz {wildcards.SM}.R2.fq.gz 2>> {wildcards.SM}.log | samblaster -a --addMateTags | samtools view -h1 --threads {threads} -bS > {wildcards.SM}.aln.bam
        rm -f {wildcards.SM}.R1.fq.gz {wildcards.SM}.R2.fq.gz
        samtools sort  -O bam -l 1 --write-index -@ {threads} {wildcards.SM}.aln.bam -o {output.bam}##idx##{output.index}
        """


# BQSR_Loci=$(echo -L chr{1..22} | sed 's/ / -L /g' | sed 's/-L -L/-L/g')


rule recalibrate:
    input:
        bam="{SM}.sorted.bam",
        bamindex="{SM}.sorted.bam.bai",
    output:
        cram="{SM}.final-gatk.cram",
        index="{SM}.final-gatk.cram.crai",
        metrics_dedup="{SM}.dedupMetrics.txt",
    log:
        dedup="{SM}.dedup.log",
        bqsr="{SM}.bqsr.log",
    threads: 1
    conda:
        "envs/optimised.yaml"
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
        bqsr_loci="-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22",
        Mills=config["Mills"],
        KnownIndels=config["KnownIndels"],
    benchmark:
        "benchmarks/{SM}.recalibrate_opt.benchmark.txt"
    shell:
        """
        wrk="."

        picard -Djava.io.tmpdir=${{wrk}} MarkDuplicates I={input.bam} AS=true O={wildcards.SM}.dedup.bam METRICS_FILE={output.metrics_dedup} QUIET=true COMPRESSION_LEVEL=0 2>> {log.dedup}
        #why sort this stuff
        samtools sort -m 6G -O bam -l 1 --write-index  -@ {threads} {wildcards.SM}.dedup.bam -o {wildcards.SM}.dedup-sorted.bam##bai##{wildcards.SM}.dedup-sorted.bam.bai 
        rm {wildcards.SM}.dedup.bam
        echo "removed dedup.bam"

        gatk3 -Djava.io.tmpdir=${{wrk}} -T BaseRecalibrator -I {wildcards.SM}.dedup-sorted.bam -R {params.ref} -o {wildcards.SM}.recal -nct {threads} --downsample_to_fraction .1 {params.bqsr_loci} -knownSites {params.dbSNP} -knownSites {params.Mills} -knownSites {params.KnownIndels} 2>> {log.bqsr}
        echo "BRC done"

        gatk3 -Djava.io.tmpdir=${{wrk}}  -T PrintReads -I {wildcards.SM}.dedup-sorted.bam -R {params.ref} -nct {threads} --BQSR {wildcards.SM}.recal -o {wildcards.SM}.bqsr.bam --globalQScorePrior -1.0 --preserve_qscores_less_than 6 --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30 --disable_indel_quals 2>> 2{log.bqsr}
        rm {wildcards.SM}.dedup-sorted.bam
        #next line just might be a view (why sort an sorted file?)
        samtools sort -m 6G -O CRAM --reference {params.ref}  --write-index -@ {threads} {wildcards.SM}.bqsr.bam -o {output.cram}

        """


rule HaplotyperWGS:
    input:
        "{SM}.final-gatk.cram",
    output:
        "",
    log:
        "{SM}.vcf.log",
    conda:
        "envs/optimised.yaml"
    threads: 2
    params:
        ref=config["fasta"],
        dbSNP=config["dbSNP"],
    shell:
        """
        wrk=".""
        # Call variants chr1
        echo -e "\\n\\nCalling chr1\\n" >> {log}
        mkdir -p ${wrk}/chr1
        cd ${wrk}/chr1
        /usr/bin/time gatk4 -Djava.io.tmpdir=${wrk}/chr1 HaplotypeCaller -R {params.ref}  --dbsnp {params.dbSMP}  -I {input} -O ${wrk}/chr1/${SM}.chr1.g.vcf.gz -L chr1 -ERC GVCF --native-pair-hmm-threads {threads} &>> {log}


        # Initialize gVCF and remove temp chrom dir
        if [ -f ${wrk}/chr1/${SM}.chr1.g.vcf.gz ]
        then
            zcat ${wrk}/chr1/${SM}.chr1.g.vcf.gz > ${wrk}/${SM}.g.vcf
            cd $wrk
            rm -fr chr1
        else
            echo -e "\\nError during variant, exiting\\n"
            cd ..
            rm -fr ${SM}
            exit 1
        fi


        # Iteratively haplotype caller and append to gVCF
        for chrom in chr{2..22} chr{X..Y}
            do

            # Call variants
            echo -e "\\n\\nCalling ${chrom}\\n" >> {log}
            mkdir -p ${wrk}/${chrom}
            cd ${wrk}/${chrom}
            /usr/bin/time gatk4 -Djava.io.tmpdir=${wrk}/${chrom} -jar  HaplotypeCaller -R {params.ref} --dbsnp {params.dbSNP} -I {input} -O ${wrk}/${chrom}/${SM}.${chrom}.g.vcf.gz -L ${chrom} -ERC GVCF --native-pair-hmm-threads 2 &>> ${wrk}/${SM}.vcf.log

            # Append to gVCF
            zcat ${wrk}/${chrom}/${SM}.${chrom}.g.vcf.gz | grep -v "#" >> {${wrk}/${SM}.g.vcf}
            cd $wrk
            rm -fr ${chrom}
        done


        # Compress and index
        ${BGZIP} -c ${wrk}/${SM}.g.vcf > {output.gvcf}
        ${TABIX} -p vcf -f {output.gvcf}


        # Calculate md5sum and check gVCF

        acc=$(basename {output.gvcf} | cut -d \. -f 1)
        base=$(basename {output.gvcf})


        # Sanity check gVCF
        iid=$(zcat {output.gvcf} | head -n 10000 | grep "#CHROM" | cut -f 10)
        size=$(du -sh {output.gvcf} | awk '{print $1}')
        ${TABIX} -R ${tgt} {output.gvcf}| bgzip -c > ${base}-exome-query.vcf.gz

        length=$(zcat ${base}-exome-query.vcf.gz | wc -l)
        width=$(zcat {output.gvcf} | grep -v "\#" | awk '{print NF}' | sort | uniq -c | awk '{print $2}' | sed 's/\n/,/g')
        NVar=$(zgrep -c "MQ" ${base}-exome-query.vcf.gz)

        # Summarise variants
        GQ20=$(zgrep "MQ" ${base}-exome-query.vcf.gz | cut -f 10 | cut -d \: -f 4 | sort -n | awk '$1 > 20 {print}' | wc -l)
        GQ60=$(zgrep "MQ" ${base}-exome-query.vcf.gz | cut -f 10 | cut -d \: -f 4 | sort -n | awk '$1 > 60 {print}' | wc -l)
        GQ90=$(zgrep "MQ" ${base}-exome-query.vcf.gz | cut -f 10 | cut -d \: -f 4 | sort -n | awk '$1 > 90 {print}' | wc -l)
        variantSummary=${GQ20},${GQ60},${GQ90}
        rm -f ${base}-exome-query.vcf.gz

        # Create table
        echo -e "IID\\tAccession\\tgVCF\\tDisk_Usage\\tWidth\\tLength\\tN_Variants\\tN_dbSNP_Calls\\tGenome_GQ_Summary(GT_20,GT_60,GT_90)\\tVariant_GQ_Summary(GT_20,GT_60,GT_90)" > ${base}_checks.tsv
        echo -e "${iid}\\t${acc}\\t${base}\\t${size}\\t${width}\\t${length}\\t${NVar}\\t${variantSummary}" >> ${base}_checks.tsv
        """
