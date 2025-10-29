def post_map_kraken_report(outdir, prefix):
    """
    Create kraken report files for bowtie2 which includes 
    top 3 species names and their relative abundances. 
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    kraken_dir = os.path.join(outdir, 'kraken2', 'post-binning')

    bowtie2_output_data = []

    def get_top_species(report_file):
        # If file missing or empty -> return 3x (NA, 0.0)
        if not os.path.exists(report_file) or os.path.getsize(report_file) == 0:
            return [('NA', 0.0), ('NA', 0.0), ('NA', 0.0)]

        species_list = []
        with open(report_file, 'r') as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    continue
                rank_code = cols[3].strip()
                name = cols[5].strip()
                if rank_code == 'S':  # Species level
                    try:
                        percentage = float(cols[0].strip())
                        species_list.append((name, percentage))
                    except ValueError:
                        continue

        species_list.sort(key=lambda x: x[1], reverse=True)

        # Pad with NA until 3
        while len(species_list) < 3:
            species_list.append(('NA', 0.0))
        return species_list[:3]

    for sample in os.listdir(kraken_dir):
        bowtie2_file = os.path.join(kraken_dir, sample, f"{sample}_bowtie2.tsv")

        bowtie2_species = get_top_species(bowtie2_file)

        bowtie2_row = [sample]
        for name, perc in bowtie2_species:
            bowtie2_row.extend([name, perc])
        bowtie2_output_data.append(bowtie2_row)

    # Write Bowtie2 output
    bowtie2_output_csv = os.path.join(report_dir, f'{prefix}_post_map_Final_Kraken_Report.csv')
    with open(bowtie2_output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'sample',
            'post_map_top_species_1', 'post_map_relative_abundance_1',
            'post_map_top_species_2', 'post_map_relative_abundance_2',
            'post_map_top_species_3', 'post_map_relative_abundance_3'
        ])
        writer.writerows(bowtie2_output_data)

def pre_map_kraken_report(outdir, prefix):
    """
    Create kraken report file for pre-binned data which includes
    top 3 species names and their relative abundances.
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    kraken_dir = os.path.join(outdir, 'kraken2', 'pre-binning')

    output_data = []

    def get_top_species(report_file):
        species_list = []
        with open(report_file, 'r') as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    continue
                rank_code = cols[3].strip()
                name = cols[5].strip()
                if rank_code == 'S':  # Species
                    try:
                        percentage = float(cols[0].strip())
                        species_list.append((name, percentage))
                    except ValueError:
                        continue
        species_list.sort(key=lambda x: x[1], reverse=True)
        while len(species_list) < 3:
            species_list.append(('NA', 0.0))
        return species_list[:3]

    for sample in os.listdir(kraken_dir):
        file = os.path.join(kraken_dir, sample, f"{sample}.tsv")
       
        species = get_top_species(file)
        
        row = [sample]
        for name, perc in species:
            row.extend([name, perc])
        output_data.append(row)


    output_csv = os.path.join(report_dir, f'{prefix}_pre_map_Final_Kraken_Report.csv')
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'sample',
            'pre_map_top_species_1', 'pre_map_relative_abundance_1',
            'pre_map_top_species_2', 'pre_map_relative_abundance_2',
            'pre_map_top_species_3', 'pre_map_relative_abundance_3'
        ])
        writer.writerows(output_data)

def mapped_reads(outdir, prefix):
    """
    Read in bowtie2 mapping stats to extract mapped read counts.
    """
    prefix = prefix.pop() 
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")

    bowtie2_dir = os.path.join(outdir, "bowtie2")

    bowtie2_rows = []

    def extract_reads(report_path): 
        # total = None
        mapped = None
        if os.path.exists(report_path):
            with open(report_path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    label = parts[2].lower()
                    # if label.startswith("total"):
                    #     try:
                    #         total = int(parts[0])
                    #     except ValueError:
                    #         pass
                    if label.startswith("mapped"):
                        try:
                            mapped = int(parts[0])
                        except ValueError:
                            pass
        return mapped

    # Process bowtie2 samples
    for sample in os.listdir(bowtie2_dir):
        bowtie2_report = os.path.join(bowtie2_dir, sample, f"{sample}_bowtie2_full_alignment_info.tsv")
        mapped = extract_reads(bowtie2_report)
        if mapped is not None:
            bowtie2_rows.append([sample, mapped])

    # Write Final Bowtie2 file
    final_bowtie2_csv = os.path.join(report_dir, f"{prefix}_Final_bowtie2_mapped_reads.csv")
    with open(final_bowtie2_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample", "bowtie2_mapped_reads"])
        writer.writerows(bowtie2_rows)


def reads_from_kneaddata(outdir, prefix):
    """
    Read in kneaddata files to extract read counts at different stages:
    initial, after trimming, and after human removal.
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    kneaddata_dir = os.path.join(outdir, "kneaddata")

    rows = []
    for sample in os.listdir(kneaddata_dir):
        log_file = os.path.join(kneaddata_dir, sample, f"{sample}_R1_kneaddata.log")
        if not os.path.exists(log_file):
            continue

        initial = after_trimming = after_human_removal = 0
        with open(log_file) as f:
            for line in f:
                if "READ COUNT: raw pair1" in line:
                    initial = int(float(line.strip().split(":")[-1]))
                elif "READ COUNT: trimmed pair1" in line:
                    after_trimming = int(float(line.strip().split(":")[-1]))
                elif "READ COUNT: final pair1" in line:
                    after_human_removal = int(float(line.strip().split(":")[-1]))

        rows.append([sample, initial, after_trimming, after_human_removal])

    # Write one combined CSV
    out_csv = os.path.join(report_dir, f"{prefix}_Final_kneaddata_read_counts.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "initial_read_pair", "after_trimming_read_pair", "after_human_removal_read_pair"])
        writer.writerows(rows)

def aggregate_amrfinder_results(prefix):
    """
    Aggregate AMRFinder results from Bowtie2 assemblies into a single TSV.
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "amrfinder", "bowtie2")
    all_results = []

    for tsv in glob.glob(os.path.join(bowtie2_dir, "*", "*_amrfinder.tsv")):
        sample = os.path.basename(tsv).replace("_amrfinder.tsv", "")
        try:
            df = pd.read_csv(tsv, sep="\t")
            df["sample"] = sample
            all_results.append(df)
        except Exception as e:
            print(f"Skipping {tsv}: {e}")

    if not all_results:
        print("No Bowtie2 AMRFinder results found.")
        return

    merged = pd.concat(all_results, ignore_index=True)

    # Write out combined file
    out_tsv = os.path.join(report_dir, f"{prefix}_Final_aggregated_amrfinder_results.tsv")
    merged.to_csv(out_tsv, sep="\t", index=False)

    print(f"Aggregated AMRFinder results written to: {out_tsv}")

def parse_busco_output(prefix):
    """ 
    This function parses BUSCO summary reports specifically from the Bowtie2
    directory and aggregates key metrics into a single CSV report.
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "busco", "bowtie2")
    all_results = []

    for summary_file in glob.glob(os.path.join(bowtie2_dir, "*", "busco.txt")):
        stats = {}
        sample_name = os.path.basename(os.path.dirname(summary_file))

        with open(summary_file, "r") as f:
            for line in f:
                line = line.strip()
                if re.match(r"^\d+\s+Complete BUSCOs", line):
                    stats["Complete BUSCOs"] = int(line.split()[0])
                elif re.match(r"^\d+\s+Complete and single-copy BUSCOs", line):
                    stats["Single-copy BUSCOs"] = int(line.split()[0])
                elif re.match(r"^\d+\s+Complete and duplicated BUSCOs", line):
                    stats["Duplicated BUSCOs"] = int(line.split()[0])
                elif re.match(r"^\d+\s+Fragmented BUSCOs", line):
                    stats["Fragmented BUSCOs"] = int(line.split()[0])
                elif re.match(r"^\d+\s+Missing BUSCOs", line):
                    stats["Missing BUSCOs"] = int(line.split()[0])
                elif re.match(r"^\d+\s+Total BUSCO groups searched", line):
                    stats["Total BUSCO groups searched"] = int(line.split()[0])

        if stats:  # Only add if stats were found
            all_results.append({
                "sample": sample_name,
                **stats
            })

    if not all_results:
        print("No Bowtie2 BUSCO reports found.")
        return

    # Merge into DataFrame and write to CSV
    df = pd.DataFrame(all_results)
    out_csv = os.path.join(report_dir, f"{prefix}_Final_busco_stats.csv")
    df.to_csv(out_csv, index=False)

    print(f"Aggregated Bowtie2 BUSCO stats written to: {out_csv}")

def aggregate_quast_results(prefix):
    """ 
    Aggregates key assembly statistics from QUAST results
    and combines them into a single CSV file.
    """
    prefix = prefix.pop() 
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "quast", "bowtie2")

    # Columns to extract from QUAST report
    keep_cols = [
        "# contigs (>= 1000 bp)",
        "# contigs",
        "Largest contig",
        "Total length",
        "GC (%)",
        "N50"
    ]

    all_results = []

    for tr_file in glob.glob(os.path.join(bowtie2_dir, "*", "transposed_report.tsv")):
        try:
            df = pd.read_csv(tr_file, sep="\t", engine="python")
            df.columns = df.columns.str.strip()  # Remove whitespace from column names

            # Check that all expected columns exist
            if not all(col in df.columns for col in keep_cols):
                raise ValueError(f"Missing expected columns in {tr_file}")

            # Extract relevant columns and add metadata
            df_subset = df[keep_cols].copy()
            df_subset["sample"] = os.path.basename(os.path.dirname(tr_file))

            # Reorder columns for clarity
            df_subset = df_subset[["sample"] + keep_cols]
            all_results.append(df_subset)

        except Exception as e:
            print(f"Skipping {tr_file} due to error: {e}")

    if all_results:
        merged_df = pd.concat(all_results, ignore_index=True)
        outfile = os.path.join(report_dir, f"{prefix}_Final_quast_results.csv")
        merged_df.to_csv(outfile, index=False)
        print(f"Saved Bowtie2 QUAST summary to {outfile}")
    else:
        print("No quast report file was found.")

def Summary(outdir, prefix):
    prefix = prefix.pop() 

    outdir = f"results/{prefix}/{prefix}_Report"

    # File paths
    post_map_bowtie2_kraken_file = os.path.join(outdir, f"{prefix}_post_map_Final_Kraken_Report.csv")
    pre_map_kraken_file = os.path.join(outdir, f"{prefix}_pre_map_Final_Kraken_Report.csv")
    bowtie2_mapped_file = os.path.join(outdir, f"{prefix}_Final_bowtie2_mapped_reads.csv")
    kneaddata_results = os.path.join(outdir, f"{prefix}_Final_kneaddata_read_counts.csv")
    amrfinder_aggregated_file = os.path.join(outdir, f"{prefix}_Final_aggregated_amrfinder_results.tsv")
    busco_aggregated_file = os.path.join(outdir, f"{prefix}_Final_busco_stats.csv")
    quast_aggregated_file = os.path.join(outdir, f"{prefix}_Final_quast_results.csv")
    
    # Load CSVs
    df_bowtie2_kraken = pd.read_csv(post_map_bowtie2_kraken_file)
    df_bowtie2_mapped = pd.read_csv(bowtie2_mapped_file)
    df_pre_map_kraken = pd.read_csv(pre_map_kraken_file)
    df_kneaddata_results = pd.read_csv(kneaddata_results)
    df_amrfinder = pd.read_csv(amrfinder_aggregated_file, sep="\t")
    df_busco = pd.read_csv(busco_aggregated_file)
    df_quast = pd.read_csv(quast_aggregated_file)

    # Merge on sample
    df = df_bowtie2_kraken
    df = df.merge(df_bowtie2_mapped, on='sample', how='left')
    df = df.merge(df_bowtie2_kraken, on='sample', how='left')
    df = df.merge(df_pre_map_kraken, on='sample', how='left')
    df = df.merge(df_kneaddata_results, on='sample', how='left')
    df = df.merge(df_busco, on='sample', how='left')
    df = df.merge(df_quast, on='sample', how='left')

    # Reorder columns
    final_columns = [
        'sample',
        'initial_read_pair', 'after_trimming_read_pair', 'after_human_removal_read_pair',
        'bowtie2_mapped_reads',
        'pre_map_top_species_1','pre_map_relative_abundance_1',
        'pre_map_top_species_2','pre_map_relative_abundance_2',
        'pre_map_top_species_3','pre_map_relative_abundance_3',
        'post_map_top_species_1','post_map_relative_abundance_1',
        'post_map_top_species_2','post_map_relative_abundance_2',
        'post_map_top_species_3','post_map_relative_abundance_3',
        '# contigs (>= 1000 bp)', '# contigs', 'Largest contig',
        'Total length', 'GC (%)', 'N50',
        'Complete BUSCOs', 'Single-copy BUSCOs', 'Duplicated BUSCOs',
        'Fragmented BUSCOs', 'Missing BUSCOs', 'Total BUSCO groups searched'
    ]

    # Save final summary
    qc_summary_file = os.path.join(outdir, f"{prefix}_QC_summary.csv")
    df = df[final_columns]
    df.to_csv(qc_summary_file, index=False)
    print(f"Saved QC summary to {qc_summary_file}")

rule generate_reports:
    input:
        kraken_post_bin = lambda wc: expand(f"results/{prefix}/kraken2/post-binning/{sample}/{sample}_bowtie2.tsv",prefix=PREFIX, sample=SAMPLES),
        bowtie2_report = lambda wc: expand(f"results/{prefix}/bowtie2/{sample}/{sample}_bowtie2_full_alignment_info.tsv",prefix=PREFIX, sample=SAMPLES),
    output:
        post_map_bowtie2_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_post_map_Final_Kraken_Report.csv",
        pre_map_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_pre_map_Final_Kraken_Report.csv",
        final_bowtie2_mapped_reads = "results/{prefix}/{prefix}_Report/{prefix}_Final_bowtie2_mapped_reads.csv",
        final_kneaddata_results = "results/{prefix}/{prefix}_Report/{prefix}_Final_kneaddata_read_counts.csv",  
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    benchmark:
        "benchmarks/{prefix}/generate_reports/generate_reports.benchmark.txt"
    run:
        post_map_kraken_report({params.outdir}, {params.prefix})
        mapped_reads({params.outdir}, {params.prefix})
        pre_map_kraken_report({params.outdir}, {params.prefix})
        reads_from_kneaddata({params.outdir}, {params.prefix})
        aggregate_amrfinder_results({params.prefix})
        parse_busco_output({params.prefix})
        aggregate_quast_results({params.prefix})

rule summary_report:
    input:
        post_map_bowtie2_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_post_map_Final_Kraken_Report.csv",
        final_bowtie2_mapped_reads = "results/{prefix}/{prefix}_Report/{prefix}_Final_bowtie2_mapped_reads.csv",   
        final_kneaddata_results = "results/{prefix}/{prefix}_Report/{prefix}_Final_kneaddata_read_counts.csv",   
    output:
        QC_summary = "results/{prefix}/{prefix}_Report/{prefix}_QC_summary.csv"
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    benchmark:
        "benchmarks/{prefix}/summary_report/summary_report.benchmark.txt"
    run:
        Summary({params.outdir}, {params.prefix})
    
