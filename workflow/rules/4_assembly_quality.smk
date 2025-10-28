def downsample_reads(file, file2, out1, out2, genome_size):
    file = file.pop()
    file2 = file2.pop()
    out1 = out1.pop()
    out2 = out2.pop()

    # Extract basic fastq reads stats with seqtk

    gsize = genome_size.pop()

    print("Using Genome Size: %s to calculate coverage" % gsize)
    
    seqtk_check = "/nfs/esnitkin/bin_group/seqtk/seqtk fqchk -q3 %s > %s_fastqchk.txt" % (file, file)

    print(seqtk_check)

    try:
        os.system(seqtk_check)
    except sp.CalledProcessError:
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % file, 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    print('Average Read Length: %s' % avg_len)

    print('Total number of bases in fastq: %s' % total_bases)

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    print('Original Covarage Depth: %s x' % ori_coverage_depth)

    # proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    # (nproc, err) = proc.communicate()
    # nproc = nproc.strip()

    if ori_coverage_depth >= 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        # r1_sub = "/tmp/%s" % os.path.basename(file)
        r1_sub = out1

        # Downsample using seqtk
        try:
            #print("Generating seqtk Downsampling command")
            print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file, factor, r1_sub))

            seqtk_downsample = "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (
                file, factor, r1_sub)
            os.system(seqtk_downsample)
            #call(seqtk_downsample, logger)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if file2:
            # r2_sub = "/tmp/%s" % os.path.basename(file2)
            r2_sub = out2
            try:
                print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                os.system("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                #call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (file2, factor, nproc, os.path.basename(file2)), logger)
            except sp.CalledProcessError:
                print('Error running seqtk for downsampling raw fastq reads.')
                sys.exit(1)
        else:
            r2_sub = "None"

    elif ori_coverage_depth < 100:
        r1_sub = file
        r2_sub = file2
        os.system("cp %s %s" % (file, out1))
        os.system("cp %s %s" % (file2, out2))

# Rule for downsampling reads to 100x before assembly
rule downsample:
    input:
        fq1_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_mapped_R1.fastq.gz",
        fq2_bowtie2 = "results/{prefix}/bowtie2/{sample}/{sample}_mapped_R2.fastq.gz",
    output:
        outr1_bowtie2 = "results/{prefix}/downsample/bowtie2/{sample}/{sample}_R1_bowtie2_mapped_downsampled.fastq.gz",
        outr2_bowtie2 = "results/{prefix}/downsample/bowtie2/{sample}/{sample}_R2_bowtie2_mapped_downsampled.fastq.gz",
    params:
        gsize = config["genome_size"],
    benchmark:
        "benchmarks/{prefix}/downsample/{sample}.benchmark.txt"
    run:
        downsample_reads({input.fq1_bowtie2}, {input.fq2_bowtie2}, {output.outr1_bowtie2}, {output.outr2_bowtie2}, {params.gsize})

# Assemble paired-end reads using SPAdes
rule assembly:
    input:
        fq1_bowtie2 = "results/{prefix}/downsample/bowtie2/{sample}/{sample}_R1_bowtie2_mapped_downsampled.fastq.gz",
        fq2_bowtie2 = "results/{prefix}/downsample/bowtie2/{sample}/{sample}_R2_bowtie2_mapped_downsampled.fastq.gz",
    output:
        spades_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/contigs.fasta", 
    params:
        bowtie2_out_dir = "results/{prefix}/spades/bowtie2/{sample}/",
        threads= 10,
        prefix_dir = "results/{prefix}/",
        sample = "{sample}"
    #conda:
    #    "envs/spades.yaml"
    singularity:
        "docker://staphb/spades:4.0.0"
    #envmodules:
    #    "Bioinformatics",
    #    "spades/4.0.0"
    benchmark:
        "benchmarks/{prefix}/assembly/{sample}.benchmark.txt"
    shell:
        """
        spades.py -1 {input.fq1_bowtie2} -2 {input.fq2_bowtie2} -o {params.bowtie2_out_dir} --threads {params.threads}
        """

# Filter assembled contigs to retain those longer than 1000 bp and rename headers
rule bioawk:
    input:
        spades_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/contigs.fasta", 
    output:
        spades_l1000_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/{sample}_contigs_l1000.fasta",
    params:
        out_dir_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}",
        prefix = "{sample}",
    benchmark:
        "benchmarks/{prefix}/bioawk/{sample}.benchmark.txt"
    conda:
        "envs/bioawk.yaml"
    shell:
        """
        bioawk -c fastx '{{if(length($seq) > 1000) {{print ">"$name; print $seq }}}}' {input.spades_assembly_bowtie2} > {output.spades_l1000_assembly_bowtie2} && \
        grep '>' {output.spades_l1000_assembly_bowtie2} > {params.out_dir_bowtie2}/spades_assembly_header_info.txt && \
        sed -i 's/>NODE_/>{params.prefix}_/g' {output.spades_l1000_assembly_bowtie2} && \
        sed -i 's/_length_.*_cov_.*//g' {output.spades_l1000_assembly_bowtie2}
        """

# AMR gene detection using AMRFinderPlus
rule amrfinder:
    input:
       spades_l1000_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/{sample}_contigs_l1000.fasta",
    output:
       amrfinder_bowtie2 = "results/{prefix}/amrfinder/bowtie2/{sample}/{sample}_amrfinder.tsv",
    params: 
       outdir = "results/{prefix}/amrfinder/",
       prefix = "{sample}",
    benchmark:
        "benchmarks/{prefix}/amrfinder/{sample}.benchmark.txt"
    singularity:
       "docker://staphb/ncbi-amrfinderplus:4.0.23-2025-07-16.1"
    shell:
       """
       amrfinder --plus --output {output.amrfinder_bowtie2} --debug --log {params.outdir}/{params.prefix}.log --nucleotide_output {params.outdir}/bowtie2/{params.prefix}/{params.prefix}_reported_nucl.fna -n {input.spades_l1000_assembly_bowtie2} -O Neisseria_gonorrhoeae
       """

# Species identification using Skani
rule skani:
    input:
        spades_l1000_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/{sample}_contigs_l1000.fasta",
    output:
        skani_output_bowtie2 = "results/{prefix}/skani/bowtie2/{sample}/{sample}_skani_output.txt",
    params:
        skani_ani_db = config["skani_db"],
        threads = 6
    benchmark:
        "benchmarks/{prefix}/skani/{sample}.benchmark.txt"
    #conda:
    #    "envs/skani.yaml"
    singularity:
        "docker://staphb/skani:0.2.1"
    shell:
        """
        skani search {input.spades_l1000_assembly_bowtie2} -d {params.skani_ani_db} -o {output.skani_output_bowtie2} -t {params.threads}
        """

# Assembly quality assessment using QUAST
rule quast:
    input:
        spades_l1000_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/{sample}_contigs_l1000.fasta",
    output:
        quast_report_bowtie2 = "results/{prefix}/quast/bowtie2/{sample}/report.tsv",
    params: 
        outdir_bowtie2 = "results/{prefix}/quast/bowtie2/{sample}/",
        prefix = "{sample}",
    benchmark:
        "benchmarks/{prefix}/quast/{sample}.benchmark.txt"
    #conda:
    #    "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.3.0"
    #envmodules:
    #    "Bioinformatics",
    #    "quast"
    shell:
        """
        quast.py {input.spades_l1000_assembly_bowtie2} -o {params.outdir_bowtie2} --contig-thresholds 0,1000,5000,10000,25000,50000
        quast.py {input.spades_l1000_assembly_bwa} -o {params.outdir_bwa} --contig-thresholds 0,1000,5000,10000,25000,50000
        """

# MLST typing using MLST
rule mlst:
    input:
        spades_l1000_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/{sample}_contigs_l1000.fasta",
    output:
        mlst_report_bowtie2 = "results/{prefix}/mlst/bowtie2/{sample}/report.tsv",
    params: 
        outdir = "results/{prefix}/mlst/{sample}/",
        prefix = "{sample}",
    benchmark:
        "benchmarks/{prefix}/mlst/{sample}.benchmark.txt"
    #conda:
    #    "envs/mlst.yaml"
    singularity:
        "docker://staphb/mlst:2.23.0-2024-03"
    #envmodules:
    #    "Bioinformatics",
    #    "mlst"
    shell:
        """
        mlst {input.spades_l1000_assembly_bowtie2} > {output.mlst_report_bowtie2}
        """

# Assembly completeness assessment using BUSCO
rule busco:
    input:
        spades_l1000_assembly_bowtie2 = "results/{prefix}/spades/bowtie2/{sample}/{sample}_contigs_l1000.fasta",
    output:
        busco_out_bowtie2 = "results/{prefix}/busco/bowtie2/{sample}/busco.txt",
    params: 
        outdir_bowtie2 = "results/{prefix}/busco/bowtie2/{sample}/",
        prefix = "{sample}",
        busco_out = "short_summary.specific.bacteria_odb12.{sample}.txt",
    benchmark:
        "benchmarks/{prefix}/busco/{sample}.benchmark.txt"
    #conda:
    #    "envs/busco.yaml"
    singularity:
        "docker://staphb/busco:5.8.2-prok-bacteria_odb12_2024-11-14"
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        """
        busco -f -i {input.spades_l1000_assembly_bowtie2} -m genome -l bacteria_odb12 -o {params.outdir_bowtie2}

        cp {params.outdir_bowtie2}/{params.busco_out} {params.outdir_bowtie2}/busco.txt
        """
    