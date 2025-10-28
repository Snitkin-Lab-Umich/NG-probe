# This rule performs quality control on sequencing reads after trimming using FastQC.
rule quality_aftertrim:
    input:
        r1 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R1.fastq.gz",
        r2 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R2.fastq.gz",
    output:
        aftertrim_fastqc_report_fwd = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_clumpify_R1_fastqc.html",
        aftertrim_fastqc_report_rev = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Reverse/{sample}_clumpify_R2_fastqc.html",
        fastqc = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_clumpify_R1_fastqc/fastqc_data.txt",
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/quality_aftertrim/{sample}/{sample}",
        sample = "{sample}"
    benchmark:
        "benchmarks/{prefix}/quality_aftertrim/{sample}.benchmark.txt"
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
        cd {params.outdir}_Forward/
        unzip {params.sample}_clumpify_R1_fastqc.zip
        """