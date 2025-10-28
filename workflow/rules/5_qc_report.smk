# Get top 3 species

def post_map_kraken_report(outdir, prefix):
    """
    Create Kraken report files for bowtie2 and bwa,
    including top 3 species names and their relative abundances.
    If Kraken report is empty, fill with NA species and 0 abundance.
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    kraken_dir = os.path.join(outdir, 'kraken2', 'post-binning')

    bowtie2_output_data = []
    bwa_output_data = []

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
        bwa_file     = os.path.join(kraken_dir, sample, f"{sample}_bwa.tsv")

        bowtie2_species = get_top_species(bowtie2_file)
        bwa_species     = get_top_species(bwa_file)

        bowtie2_row = [sample]
        for name, perc in bowtie2_species:
            bowtie2_row.extend([name, perc])
        bowtie2_output_data.append(bowtie2_row)

        bwa_row = [sample]
        for name, perc in bwa_species:
            bwa_row.extend([name, perc])
        bwa_output_data.append(bwa_row)

    # Write Bowtie2 output
    bowtie2_output_csv = os.path.join(report_dir, f'{prefix}_bowtie2_post_map_Final_Kraken_Report.csv')
    with open(bowtie2_output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'sample',
            'post_map_bowtie2_top_species_1', 'post_map_bowtie2_relative_abundance_1',
            'post_map_bowtie2_top_species_2', 'post_map_bowtie2_relative_abundance_2',
            'post_map_bowtie2_top_species_3', 'post_map_bowtie2_relative_abundance_3'
        ])
        writer.writerows(bowtie2_output_data)

    # Write BWA output
    bwa_output_csv = os.path.join(report_dir, f'{prefix}_bwa_post_map_Final_Kraken_Report.csv')
    with open(bwa_output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'sample',
            'post_map_bwa_top_species_1', 'post_map_bwa_relative_abundance_1',
            'post_map_bwa_top_species_2', 'post_map_bwa_relative_abundance_2',
            'post_map_bwa_top_species_3', 'post_map_bwa_relative_abundance_3'
        ])
        writer.writerows(bwa_output_data)


def pre_map_kraken_report(outdir, prefix):
    """
    Create Kraken report file pre-binned data,
    including top 3 species names and their relative abundances.
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


# Get total reads
# def total_reads(outdir, prefix):
#     prefix = prefix.pop()
#     outdir = "results/%s" % prefix
#     report_dir = str(outdir) + "/%s_Report" % prefix

#     results = []

#     quality_aftertrim_dir = os.path.join(outdir, 'quality_aftertrim') 

#     for sample in os.listdir(quality_aftertrim_dir):
#         fastqc_path = f"results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_clumpify_R1_fastqc/fastqc_data.txt"

#         if not os.path.exists(fastqc_path):
#             print(f"Skipping {sample}: fastqc file not found.")
#             continue

#         total_sequences = None

#         with open(fastqc_path, 'r') as f:
#             for line in f:
#                 if line.startswith("Total Sequences"):
#                     total_sequences = int(line.strip().split('\t')[-1])

#         if total_sequences is None:
#             print(f"Skipping {sample} could not extract Total Sequences.")
#             continue

#         results.append({'sample': sample, 'total_reads': total_sequences})

#     df_results = pd.DataFrame(results)
#     results_file_path = os.path.join(report_dir, f'{prefix}_Total_reads.csv') 
    
#     df_results.to_csv(results_file_path, index=False)
#     print(f"Saved reads report to {results_file_path}")

# Get mapped reads
def mapped_reads(outdir, prefix):
    prefix = prefix.pop() 
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "bowtie2")
    bwa_dir = os.path.join(outdir, "bwa")

    bowtie2_rows = []
    bwa_rows = []

    def extract_reads(report_path): # 
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

    # Process bwa samples
    for sample in os.listdir(bwa_dir):
        bwa_report = os.path.join(bwa_dir, sample, f"{sample}_bwa_full_mapping_stats.tsv")
        mapped = extract_reads(bwa_report)
        if mapped is not None:
            bwa_rows.append([sample, mapped])

    # Write Final Bowtie2 file
    final_bowtie2_csv = os.path.join(report_dir, f"{prefix}_Final_bowtie2_mapped_reads.csv")
    with open(final_bowtie2_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample", "bowtie2_mapped_reads"])
        writer.writerows(bowtie2_rows)

    # Write Final BWA file
    final_bwa_csv = os.path.join(report_dir, f"{prefix}_Final_bwa_mapped_reads.csv")
    with open(final_bwa_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample", "bwa_mapped_reads"])
        writer.writerows(bwa_rows)

# This function will read in kneaddata log file and get initial_read_pair, after_trimming, after_human_removal
# initial_reads = initial
# after_trimming = trimmed_pairs
# after_human_removal = final_pair1

def reads_from_kneaddata(outdir, prefix):
    """
    Parse kneaddata log files and write CSV with:
    sample, initial_read_pair, after_trimming_read_pair, after_human_removal_read_pair
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
    out_csv = os.path.join(report_dir, f"{prefix}_kneaddata_reads.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "initial_read_pair", "after_trimming_read_pair", "after_human_removal_read_pair"])
        writer.writerows(rows)

# def aggregate_resfinder_results(prefix):
#     """
#     Aggregate ResFinder + PointFinder results across bowtie2 and bwa into one combined CSV.
#     Columns are suffixed with _bowtie2 or _bwa.
#     """
#     prefix = prefix.pop()
#     outdir = f"results/{prefix}"
#     report_dir = os.path.join(outdir, f"{prefix}_Report")
#     os.makedirs(report_dir, exist_ok=True)

#     bowtie2_dir = os.path.join(outdir, "resfinder", "bowtie2")
#     bwa_dir = os.path.join(outdir, "resfinder", "bwa")

#     all_samples = sorted(os.listdir(bowtie2_dir))
#     merged_data = []

#     for sample in all_samples:
#         row = {"sample": sample}

#         # --- ResFinder (tabular) ---
#         for aligner, d in [("bowtie2", bowtie2_dir), ("bwa", bwa_dir)]:
#             resfinder_file = os.path.join(d, sample, "ResFinder_results_tab.txt")
#             if os.path.exists(resfinder_file) and os.path.getsize(resfinder_file) > 0:
#                 df = pd.read_csv(resfinder_file, sep="\t")
#                 # Flatten into strings (gene|identity|coverage…)
#                 resfinder_str = ";".join(df.astype(str).apply(lambda x: "|".join(x), axis=1))
#             else:
#                 resfinder_str = "NA"
#             row[f"ResFinder_{aligner}"] = resfinder_str

#         # --- PointFinder ---
#         for aligner, d in [("bowtie2", bowtie2_dir), ("bwa", bwa_dir)]:
#             pointfinder_file = os.path.join(d, sample, "PointFinder_results.txt")
#             if os.path.exists(pointfinder_file) and os.path.getsize(pointfinder_file) > 0:
#                 df = pd.read_csv(pointfinder_file, sep="\t", header=None)
#                 # Flatten into strings (mutation|aa-change|drug…)
#                 pointfinder_str = ";".join(df.astype(str).apply(lambda x: "|".join(x), axis=1))
#             else:
#                 pointfinder_str = "NA"
#             row[f"PointFinder_{aligner}"] = pointfinder_str

#         merged_data.append(row)

#     # Save merged summary
#     merged_df = pd.DataFrame(merged_data)
#     outfile = os.path.join(report_dir, f"{prefix}_combined_ResFinder_results.csv")
#     merged_df.to_csv(outfile, index=False)


def aggregate_resfinder_results(prefix):
    """
    Aggregate ResFinder + PointFinder results across bowtie2 and bwa into one long-format CSV.

    Output columns:
      sample, aligner, tool, determinant_info (expanded from each file)
    """
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "resfinder", "bowtie2")
    bwa_dir = os.path.join(outdir, "resfinder", "bwa")

    all_samples = sorted(os.listdir(bowtie2_dir))
    rows = []

    for sample in all_samples:
        for aligner, d in [("bowtie2", bowtie2_dir), ("bwa", bwa_dir)]:
            # --- ResFinder ---
            resfinder_file = os.path.join(d, sample, "ResFinder_results_tab.txt")
            if os.path.exists(resfinder_file) and os.path.getsize(resfinder_file) > 0:
                df = pd.read_csv(resfinder_file, sep="\t")
                for _, row in df.iterrows():
                    rows.append({
                        "sample": sample,
                        "aligner": aligner,
                        "tool": "ResFinder",
                        "Resistance_gene": row.get("Resistance gene", "NA"),
                        "Identity": row.get("Identity", "NA"),
                        "Alignment_Length/Gene_Length": row.get("Alignment Length/Gene Length", "NA"),
                        "Coverage": row.get("Coverage", "NA"),
                        "Position in reference": row.get("Position in reference", "NA"),
                        "Contig": row.get("Contig", "NA"),
                        "Position in contig": row.get("Position in contig", "NA"),
                        "Phenotype": row.get("Phenotype", "NA"),
                        "Accession_no": row.get("Accession no.", "NA")
                    })
            else:
                rows.append({
                    "sample": sample,
                    "aligner": aligner,
                    "tool": "ResFinder",
                    "Resistance_gene": "NA",
                    "Identity": "NA",
                    "Alignment_Length/Gene_Length": "NA",
                    "Coverage": "NA",
                    "Position in reference": "NA",
                    "Contig": "NA",
                    "Position in contig": "NA",
                    "Phenotype": "NA",
                    "Accession_no": "NA"
                })

            # --- PointFinder ---
            pointfinder_file = os.path.join(d, sample, "PointFinder_results.txt")
            if os.path.exists(pointfinder_file) and os.path.getsize(pointfinder_file) > 0:
                df = pd.read_csv(pointfinder_file, sep="\t")
                # Example format: mutation | nt-change | aa-change | phenotype | PMID
                for _, row in df.iterrows():
                    rows.append({
                        "sample": sample,
                        "aligner": aligner,
                        "tool": "PointFinder",
                        "Mutation": row[0] if len(row) > 0 else "NA",
                        "Nucleotide_change": row[1] if len(row) > 1 else "NA",
                        "Amino_acid_change": row[2] if len(row) > 2 else "NA",
                        "Resistance": row[3] if len(row) > 3 else "NA",
                        "PMID": row[4] if len(row) > 4 else "NA"
                    })
            else:
                rows.append({
                    "sample": sample,
                    "aligner": aligner,
                    "tool": "PointFinder",
                    "Mutation": "NA",
                    "Nucleotide_change": "NA",
                    "Amino_acid_change": "NA",
                    "Resistance": "NA",
                    "PMID": "NA"
                })

    # Convert to dataframe and write CSV
    merged_df = pd.DataFrame(rows)
    outfile = os.path.join(report_dir, f"{prefix}_ResFinder_results_long.csv")
    merged_df.to_csv(outfile, index=False)

import os
import glob
import pandas as pd

def aggregate_amrfinder_results(prefix):
    
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "amrfinder", "bowtie2")
    bwa_dir = os.path.join(outdir, "amrfinder", "bwa")

    all_results = []

    for aligner, indir in [("bowtie2", bowtie2_dir), ("bwa", bwa_dir)]:
        for tsv in glob.glob(os.path.join(indir, "*", "*_amrfinder.tsv")):
            sample = os.path.basename(tsv).replace("_amrfinder.tsv", "")
            try:
                df = pd.read_csv(tsv, sep="\t")
                df["sample"] = sample
                df["aligner"] = aligner
                all_results.append(df)
            except Exception as e:
                print(f"Skipping {tsv}: {e}")

    if not all_results:
        print("No AMRFinder results found.")
        return

    merged = pd.concat(all_results, ignore_index=True)

    # Write out combined file
    out_tsv = os.path.join(report_dir, f"{prefix}_amrfinder_aggregated.tsv")
    merged.to_csv(out_tsv, sep="\t", index=False)

    print(f"Aggregated AMRFinder results written to: {out_tsv}")

def parse_busco_output(prefix):
    # Handle prefix input
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "busco", "bowtie2")
    bwa_dir = os.path.join(outdir, "busco", "bwa")

    all_results = []

    for aligner, indir in [("bowtie2", bowtie2_dir), ("bwa", bwa_dir)]:
        for summary_file in glob.glob(os.path.join(indir, "*", "busco.txt")):
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
                    "Sample": sample_name,
                    "Aligner": aligner,
                    **stats
                })

    if not all_results:
        print("No BUSCO reports found.")
        return

    # Merge into DataFrame
    df = pd.DataFrame(all_results)

    # Write to CSV
    out_csv = os.path.join(report_dir, f"{prefix}_busco_stats.csv")
    df.to_csv(out_csv, index=False)

    print(f"Aggregated BUSCO stats written to: {out_csv}")

def aggregate_quast_results(prefix):
    prefix = prefix.pop() 

    outdir = f"/scratch/esnitkin_root/esnitkin1/dhatrib/NG-probe/results/{prefix}"
    report_dir = os.path.join(outdir, f"{prefix}_Report")
    os.makedirs(report_dir, exist_ok=True)

    bowtie2_dir = os.path.join(outdir, "quast", "bowtie2")
    bwa_dir = os.path.join(outdir, "quast", "bwa")

    # columns we care about
    keep_cols = [
        "# contigs (>= 1000 bp)",
        "# contigs",
        "Largest contig",
        "Total length",
        "GC (%)",
        "N50"
    ]

    all_results = []

    for aligner, indir in [("bowtie2", bowtie2_dir), ("bwa", bwa_dir)]:
        for tr_file in glob.glob(os.path.join(indir, "*", "transposed_report.tsv")):
            try:
                df = pd.read_csv(tr_file, sep="\t", engine="python")
                df.columns = df.columns.str.strip()  # Remove whitespace from column names

                # Check if all required columns exist
                if not all(col in df.columns for col in keep_cols):
                    raise ValueError(f"Columns {keep_cols} not all present in {tr_file}")

                # Add sample and aligner info
                df_subset = df[keep_cols].copy()
                df_subset["Sample"] = os.path.basename(os.path.dirname(tr_file))
                df_subset["Aligner"] = aligner
                
                # Reorder columns so Sample and Aligner come first
                df_subset = df_subset[["Sample", "Aligner"] + keep_cols]

                all_results.append(df_subset)

            except Exception as e:
                print(f"Skipping {tr_file} due to error: {e}")

    # print(all_results)
    if all_results:
        merged_df = pd.concat(all_results, ignore_index=True)
        outfile = os.path.join(report_dir, f"{prefix}_quast_results.csv")
        merged_df.to_csv(outfile, index=False)
        print(f"Saved QUAST summary to {outfile}")
    else:
        print("No transposed_report.txt files found.")

def Summary(outdir, prefix):
    prefix = prefix.pop() 

    outdir = f"results/{prefix}/{prefix}_Report"

    # File paths
    post_map_bwa_kraken_file = os.path.join(outdir, f"{prefix}_bwa_post_map_Final_Kraken_Report.csv")
    post_map_bowtie2_kraken_file = os.path.join(outdir, f"{prefix}_bowtie2_post_map_Final_Kraken_Report.csv")
    pre_map_kraken_file = os.path.join(outdir, f"{prefix}_pre_map_Final_Kraken_Report.csv")
    bwa_mapped_file = os.path.join(outdir, f"{prefix}_Final_bwa_mapped_reads.csv")
    bowtie2_mapped_file = os.path.join(outdir, f"{prefix}_Final_bowtie2_mapped_reads.csv")
    kneaddata_results = os.path.join(outdir, f"{prefix}_kneaddata_reads.csv")
    # total_reads_file = os.path.join(outdir, f"{prefix}_Total_reads.csv")

    # Load CSVs
    df_bwa_kraken = pd.read_csv(post_map_bwa_kraken_file)
    df_bowtie2_kraken = pd.read_csv(post_map_bowtie2_kraken_file)
    df_bwa_mapped = pd.read_csv(bwa_mapped_file)
    df_bowtie2_mapped = pd.read_csv(bowtie2_mapped_file)
    df_pre_map_kraken = pd.read_csv(pre_map_kraken_file)
    df_kneaddata_results = pd.read_csv(kneaddata_results)

    # Merge on sample
    # df = df_total_reads.rename(columns={'total_reads': 'bwa_total_reads'})
    df = df_bwa_mapped
    df = df.merge(df_bowtie2_mapped, on='sample', how='left')
    df = df.merge(df_bwa_kraken, on='sample', how='left')
    df = df.merge(df_bowtie2_kraken, on='sample', how='left')
    df = df.merge(df_pre_map_kraken, on='sample', how='left')
    df = df.merge(df_kneaddata_results, on='sample', how='left')
    # Calculate mapped percent columns
    # df['bwa_mapped_percent'] = (df['bwa_mapped_reads'] / df['bwa_total_reads'] * 100).round(2)
    # df['bowtie2_mapped_percent'] = (df['bowtie2_mapped_reads'] / df['bwa_total_reads'] * 100).round(2)

    # Reorder columns
    final_columns = [
        'sample',
        'initial_read_pair', 'after_trimming_read_pair', 'after_human_removal_read_pair',
        'bwa_mapped_reads', 'bowtie2_mapped_reads',
        'pre_map_top_species_1','pre_map_relative_abundance_1',
        'pre_map_top_species_2','pre_map_relative_abundance_2',
        'pre_map_top_species_3','pre_map_relative_abundance_3',
        'post_map_bwa_top_species_1','post_map_bwa_relative_abundance_1',
        'post_map_bwa_top_species_2','post_map_bwa_relative_abundance_2',
        'post_map_bwa_top_species_3','post_map_bwa_relative_abundance_3',
        'post_map_bowtie2_top_species_1','post_map_bowtie2_relative_abundance_1',
        'post_map_bowtie2_top_species_2','post_map_bowtie2_relative_abundance_2',
        'post_map_bowtie2_top_species_3','post_map_bowtie2_relative_abundance_3'
    ]

    # Add bowtie2 total reads column as same as bwa total (or compute separately if available)
    # df['bowtie2_total_reads'] = df['bwa_total_reads']

    # Save final summary
    qc_summary_file = os.path.join(outdir, f"{prefix}_QC_summary.csv")
    df = df[final_columns]
    df.to_csv(qc_summary_file, index=False)
    print(f"Saved QC summary to {qc_summary_file}")

rule generate_reports:
    input:
        report_bwa = expand("results/{prefix}/kraken2/post-binning/{sample}/{sample}_bwa.tsv",prefix=PREFIX, sample=SAMPLES),
        report_bowtie2 = expand("results/{prefix}/kraken2/post-binning/{sample}/{sample}_bowtie2.tsv",prefix=PREFIX, sample=SAMPLES),
        bowtie2 = expand("results/{prefix}/bowtie2/{sample}/{sample}_bowtie2_full_alignment_info.tsv",prefix=PREFIX, sample=SAMPLES),
        bwa = expand("results/{prefix}/bwa/{sample}/{sample}_bwa_full_mapping_stats.tsv",prefix=PREFIX, sample=SAMPLES),
    output:
        post_map_bwa_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_bwa_post_map_Final_Kraken_Report.csv",
        post_map_bowtie2_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_bowtie2_post_map_Final_Kraken_Report.csv",
        pre_map_bwa_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_pre_map_Final_Kraken_Report.csv",
        final_bwa_mapped_reads = "results/{prefix}/{prefix}_Report/{prefix}_Final_bwa_mapped_reads.csv",
        final_bowtie2_mapped_reads = "results/{prefix}/{prefix}_Report/{prefix}_Final_bowtie2_mapped_reads.csv",  
        final_kneaddata_results = "results/{prefix}/{prefix}_Report/{prefix}_kneaddata_reads.csv",  
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
        aggregate_resfinder_results({params.prefix})
        aggregate_amrfinder_results({params.prefix})
        parse_busco_output({params.prefix})
        aggregate_quast_results({params.prefix})

rule summary_report:
    input:
        post_map_bwa_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_bwa_post_map_Final_Kraken_Report.csv",
        post_map_bowtie2_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_bowtie2_post_map_Final_Kraken_Report.csv",
        pre_map_bwa_kraken_final_out = "results/{prefix}/{prefix}_Report/{prefix}_pre_map_Final_Kraken_Report.csv",
        final_bwa_mapped_reads = "results/{prefix}/{prefix}_Report/{prefix}_Final_bwa_mapped_reads.csv",
        final_bowtie2_mapped_reads = "results/{prefix}/{prefix}_Report/{prefix}_Final_bowtie2_mapped_reads.csv",   
        final_kneaddata_results = "results/{prefix}/{prefix}_Report/{prefix}_kneaddata_reads.csv",   
    output:
        "results/{prefix}/{prefix}_Report/{prefix}_QC_summary.csv"
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    benchmark:
        "benchmarks/{prefix}/summary_report/summary_report.benchmark.txt"
    run:
        Summary({params.outdir}, {params.prefix})
    