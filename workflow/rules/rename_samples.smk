# remame.txt should be a tab seperated file with in the first column the old ID and the second column the new id
mapx = {x.split()[0]: x.split()[1] for x in open("rename_all.txt").readlines()}
mapxreverse = {v: k for k, v in mapx.items()}
print(mapx)
print(mapxreverse)


localrules:
    UploadtoGridd,


rule rename_cram:
    input:
        cram=lambda wildcards: f"{mapxreverse[wildcards.newname]}.final-gatk.cram",
    output:
        cram=temp("renamed/{newname}.final-gatk.cram"),
        index=temp("renamed/{newname}.final-gatk.cram.crai"),
    resources:
        mem_mb=1000,
    params:
        old=lambda wildcards: mapxreverse[wildcards.newname],
    conda:
        "../envs/samtools1.15.yaml"
    threads: 1
    shell:
        """
        samtools view -H {input.cram} |sed -e "/^@RG/s/{params.old}/{wildcards.newname}/g" >temp.{params.old}.head
        samtools reheader temp.{params.old}.head {input.cram} >{output.cram}
        rm temp.{params.old}.head
        samtools index {output.cram}
        """


rule rename_vcf:
    input:
        vcf=lambda wildcards: mapxreverse[wildcards.newname] + ".g.vcf.gz",
    output:
        vcf="renamed/{newname}.g.vcf.gz",
        index="renamed/{newname}.g.vcf.gz.tbi",
    resources:
        mem_mb=1000,
    conda:
        "../envs/optimised.yaml"
    threads: 1
    shell:
        """
        bcftools reheader -s rename.txt  {input.vcf} -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule all:
    input:
        expand("renamed/{newname}.g.vcf.gz", newname=mapx.values()),
        expand("renamed/{newname}.final-gatk.cram", newname=mapx.values()),
    output:
        "rename.done",
    shell:
        "touch {output}"


rule UploadtoGridd:
    input:
        "renamed/{newname}.final-gatk.cram",
        "renamed/{newname}.final-gatk.cram.crai",
    output:
        "renamed/{newname}.uploaded",
    shell:
        """
        set +u
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate /cvmfs/softdrive.nl/projectmine_sw/software
            set -u
            cram=$(basename {input[0]})
            pm_copy  {input[0]} /projectmine-nfs/Realignment/hg38/Netherlands/Tape/NovaSeqBroad/{wildcards.newname}/${{cram}}
            index=$(basename {input[1]})
            pm_copy  {input[1]} /projectmine-nfs/Realignment/hg38/Netherlands/Tape/NovaSeqBroad/{wildcards.newname}/${{index}}
            set +u
            conda deactivate
            set -u
            touch {output}
        """


rule upload_all:
    input:
        expand("renamed/{newname}.uploaded", newname=mapx.values()),
    output:
        "upload_all.done",
    shell:
        "touch {output}"
