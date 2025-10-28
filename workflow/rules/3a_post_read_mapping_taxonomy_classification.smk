# Taxonomic classification of reads using Kraken2 post-binning step.
rule kraken2_post_mapping:
  input:
    r1_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_mapped_R1.fastq.gz",
    r2_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_mapped_R2.fastq.gz",
  output:
    report_bowtie2 = "results/{prefix}/kraken2/post-binning/{sample}/{sample}_bowtie2.tsv",
    output_bowtie2 = "results/{prefix}/kraken2/post-binning/{sample}/{sample}_bowtie2.out",
  params:
    db = config["kraken_db"]
  benchmark:
    "benchmarks/{prefix}/kraken2/post-binning/{sample}.benchmark.txt"
  singularity:
    "docker://staphb/kraken2:2.1.6"
  shell:
    """
    kraken2 \
        --db {params.db} \
        --paired {input.r1_bowtie2} {input.r2_bowtie2} \
        --report {output.report_bowtie2} \
        --output {output.output_bowtie2}
    """
