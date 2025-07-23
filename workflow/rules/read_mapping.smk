# How to redirect the outputs of these rules to paired end fastq files?
rule bowtie2:
    input:
       r1 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R1.fastq.gz",
        r2 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R2.fastq.gz",
    output:
        bam = "results/{prefix}/bowtie2/{sample}.bam"
    params:
        index = "/path/to/ng_reference"
    threads: 
        config["ncores"]
    singularity:
        "docker://staphb/bowtie2:2.5.4" 
    shell:
        """
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} \
        | samtools view -Sb - \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """


rule bwa_mem:
    input:
        r1 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R1.fastq.gz",
        r2 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R2.fastq.gz"
    output:
        bam = "results/{prefix}/bwa/{sample}.bam"
    params:
        index = "/path/to/ng_reference" 
    threads: 
        config["ncores"]
    singularity:
        "docker://staphb/bwa:0.7.19" 
    shell:
        """
        bwa mem -t {threads} {params.index} {input.r1} {input.r2} \
        | samtools view -Sb - \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """