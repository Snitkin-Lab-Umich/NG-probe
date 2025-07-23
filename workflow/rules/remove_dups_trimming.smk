# Remove PCR duplicates and trim reads using fastp
# Picard markdups requires bam file

rule clumpify:
    input:
        r1 = "results/{prefix}/kneaddata/{sample}_kneaddata.trimmed.1.fastq.gz",
        r2 = "results/{prefix}/kneaddata/{sample}_kneaddata.trimmed.2.fastq.gz"
    output:
        r1 = "results/{prefix}/clumpify/{sample}_clumpify_R1.fastq.gz",
        r2 = "results/{prefix}/clumpify/{sample}_clumpify_R2.fastq.gz"
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

rule fastp:
    input:
        r1 = "results/{prefix}/clumpify/{sample}_clumpify_R1.fastq.gz",
        r2 = "results/{prefix}/clumpify/{sample}_clumpify_R2.fastq.gz"
    output:
        r1 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R1.fastq.gz",
        r2 = "results/{prefix}/fastp/{sample}/{sample}_human_dups_removed_R2.fastq.gz",
        html = "results/{prefix}/fastp/{sample}/fastp_report.html",
        json = "results/{prefix}/fastp/{sample}/fastp_report.json"
    singularity:
        "docker://staphb/fastp:1.0.1" 
    shell:
        """
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --html {output.html} --json {output.json} 
        """
