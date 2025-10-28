
# This rule performs host depletion using kneaddata to clean the input reads before further processing.
# Download database first: kneaddata_database --download human_genome bowtie2 /nfs/turbo/umms-esnitkin/database/human_db # Downloaded July 29,2025
rule kneaddata: 
    input:
        r1 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R1.fastq.gz"),
        r2 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R2.fastq.gz")
    output:
        r1 = "results/{prefix}/kneaddata/{sample}/{sample}_R1_kneaddata_paired_1.fastq.gz",
        r2 = "results/{prefix}/kneaddata/{sample}/{sample}_R2_kneaddata_paired_2.fastq.gz",
    conda:
        "envs/kneaddata.yaml"
    # singularity:
    #     "docker://biobakery/kneaddata:0.10.0"
    params:
        outdir = "results/{prefix}/kneaddata/{sample}",
        db = config["human_db"],
        sample = "{sample}",
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
    benchmark:
        "benchmarks/{prefix}/kneaddata/{sample}.benchmark.txt"
    shell:
        """
        kneaddata --remove-intermediate-output --input1 {input.r1} --input2 {input.r2} --output {params.outdir} --reference-db {params.db} --trimmomatic-options "ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length}"

        gzip {params.outdir}/*.fastq

        mv {params.outdir}/{params.sample}_R1_kneaddata_paired_2.fastq.gz {params.outdir}/{params.sample}_R2_kneaddata_paired_2.fastq.gz
        """
