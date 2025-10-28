# Use clumpify from bbtools to remove duplicates from fastq files
rule clumpify:
    input:
        r1 = "results/{prefix}/kneaddata/{sample}/{sample}_R1_kneaddata_paired_1.fastq.gz",
        r2 = "results/{prefix}/kneaddata/{sample}/{sample}_R2_kneaddata_paired_2.fastq.gz"
    output:
        r1 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R1.fastq.gz",
        r2 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R2.fastq.gz"
    benchmark:
        "benchmarks/{prefix}/clumpify/{sample}.benchmark.txt"
    singularity:
         "docker://staphb/bbtools:39.26" 
    shell:
        """
        clumpify.sh \
          in1={input.r1} \
          in2={input.r2} \
          out1={output.r1} \
          out2={output.r2} \
        """
