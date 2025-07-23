# Author: Dhatri Badri

rule quality_aftertrim:
    input:
        r1 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R1.fastq.gz",
        r2 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R2.fastq.gz",
    output:
        aftertrim_fastqc_report_fwd = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_human_dups_removed_fastqc.html",
        aftertrim_fastqc_report_rev = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Reverse/{sample}_R2_human_dups_removed_fastqc.html",
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/quality_aftertrim/{sample}/{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """