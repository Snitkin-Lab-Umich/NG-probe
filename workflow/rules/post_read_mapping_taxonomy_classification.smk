# 
# These rules
rule kraken2:
    input:
        r1_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_binned_reads_R1.fastq.gz",
        r2_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_binned_reads_R2.fastq.gz",
        r1_bwa = "results/{prefix}/bwa/{sample}/{sample}_binned_reads_R1.fastq.gz",
        r2_bwa = "results/{prefix}/bwa/{sample}/{sample}_binned_reads_R2.fastq.gz",
    output:
        report_bowtie2 = "results/{prefix}/kraken2/{sample}_bowtie2.report",
        output_bowtie2 = "results/{prefix}/kraken2/{sample}_bowtie2.out",
        report_bwa = "results/{prefix}/kraken2/{sample}_bwa.report",
        output_bwa = "results/{prefix}/kraken2/{sample}_bwa.out"
    params:
        db = config["kraken_db"]
    shell:
        """
        kraken2 \
          --db {params.db} \
          --paired {input.r1_bowtie2} {input.r2_bowtie2} \
          --report {output.report_bowtie2} \
          --output {output.output_bowtie2}
        
        kraken2 \
          --db {params.db} \
          --paired {input.r1_bwa} {input.r2_bwa} \
          --report {output.report_bwa} \
          --output {output.output_bwa}
        """

rule centrifuge:
    input:
        r1_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_binned_reads_R1.fastq.gz",
        r2_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_binned_reads_R2.fastq.gz",
        r1_bwa = "results/{prefix}/bwa/{sample}/{sample}_binned_reads_R1.fastq.gz",
        r2_bwa = "results/{prefix}/bwa/{sample}/{sample}_binned_reads_R2.fastq.gz",
    output:
        report_bowtie2 = "results/{prefix}/centrifuge/{sample}_bowtie2.report",
        output_bowtie2 = "results/{prefix}/centrifuge/{sample}_bowtie2.out",
        report_bwa = "results/{prefix}/centrifuge/{sample}_bwa.report",
        output_bwa = "results/{prefix}/centrifuge/{sample}_bwa.out"
    params:
        db = "/path/to/centrifuge_db"
    shell:
        """
        centrifuge \
          -x {params.db} \
          -1 {input.r1_bowtie2} -2 {input.r2_bowtie2} \
          --report-file {output.report_bowtie2} \
          -S {output.output_bowtie2}
       
        centrifuge \
          -x {params.db} \
          -1 {input.r1_bwa} -2 {input.r2_bwa} \
          --report-file {output.report_bwa} \
          -S {output.output_bwa}
        """