# Initial taxonomic classification of reads using Kraken2 and Centrifuge pre-binning step.

rule kraken2:
    input:
        r1 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R1.fastq.gz",
        r2 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R2.fastq.gz",
    output:
        report = "results/{prefix}/kraken2/{sample}.report",
        output = "results/{prefix}/kraken2/{sample}.out"
    params:
        db = "/path/to/kraken2_db"
    shell:
        """
        kraken2 \
          --db {params.db} \
          --paired {input.r1} {input.r2} \
          --report {output.report} \
          --output {output.output}
        """

rule centrifuge:
    input:
        r1 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R1.fastq.gz",
        r2 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R2.fastq.gz",
    output:
        report = "results/{prefix}/centrifuge/{sample}.report",
        output = "results/{prefix}/centrifuge/{sample}.out"
    params:
        db = "/path/to/centrifuge_db"
    shell:
        """
        centrifuge \
          -x {params.db} \
          -1 {input.r1} -2 {input.r2} \
          --report-file {output.report} \
          -S {output.output}
        """