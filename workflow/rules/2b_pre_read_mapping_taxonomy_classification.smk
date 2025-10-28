# Initial taxonomic classification of reads using Kraken2 pre-binning step.
rule kraken2_pre_mapping:
    input:
        r1 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R1.fastq.gz"),
        r2 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R2.fastq.gz"),
    output:
        report = "results/{prefix}/kraken2/pre-binning/{sample}/{sample}.tsv",
        output = "results/{prefix}/kraken2/pre-binning/{sample}/{sample}.out"
    params:
        db = config["kraken_db"]
    benchmark:
        "benchmarks/{prefix}/kraken2/pre-binning/{sample}.benchmark.txt"
    singularity:
        "docker://staphb/kraken2:2.1.6"
    shell:
        """
        kraken2 \
          --db {params.db} \
          --paired {input.r1} {input.r2} \
          --report {output.report} \
          --output {output.output}
        """
