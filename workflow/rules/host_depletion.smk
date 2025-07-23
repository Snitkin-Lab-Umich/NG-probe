# Include one line description of each rule

# This rule performs host depletion using kneaddata to clean the input reads before further processing.
rule kneaddata: # Download database first: kneaddata_database --download human_genome bowtie2 /nfs/turbo/umms-esnitkin/database
    input:
        r1 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R1.fastq.gz"),
        r2 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R2.fastq.gz")
    output:
        r1 = "results/{prefix}/kneaddata/{sample}_kneaddata.trimmed.1.fastq.gz",
        r2 = "results/{prefix}/kneaddata/{sample}_kneaddata.trimmed.2.fastq.gz"
    singularity:
        "docker://biobakery/kneaddata:0.10.0"
    params:
        outdir = "results/{prefix}/kneaddata/{sample}",
        db = "/path/to/human_db"
    shell:
        """
        kneaddata \
          --input1 {input.r1} --input2 {input.r2} \
          --output {params.outdir} \
          --reference-db {params.db} \
        
        gzip {params.outdir}/*.fastq
        """
