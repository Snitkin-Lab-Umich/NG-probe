# NG-probe
A Snakemake pipeline for culture-free Neisseria gonorrhoeae genome reconstruction and detection of putative resistance determinants from probe-based genomic sequencing.

### Summary

In short, it performs the following steps:
1. Input is paired end raw fastq reads which will be FastQC-ed  and you wil get fastqc reports.

2. The raw fastq reads will be trimmed and human read will be removed by kneaddata.

3. The resulting files will then be rid of PCR duplicates by clumpify (bbtools)

4. The trimmed and removed duplciates paired end reads will then go through kraken to see what are the 3 most abundant species hits and will go through FastQC as well to check there are no adapter contamination,e tc. Resulting output will be kraken reports and fastqc reports respectively. 

5.  The final output files from clumpify will go through BWA/Bowtie2 for read mappign to reference genome NG

6. The resulting paired end fastqs that should ideally map to NG only will go through kraken and resfinder to find out top 3 species hits post read mapping and find any amr genes respectively. 

7. The post read mapping paired end reads will be downsampled by seqtk and assembled using spades to get a fasta file. The fasta file will then go through bioawk to only retain scaffolds >1kb. The fasta file containing contigs > than 1kb will go through skani, BUSCO, QUAST, mlst and amrfinder to get species detection a sanity check, to find out genome completedness, assembly stats,  sequence typing and find more amr genes respectively. 

8. The final QC report will be split in two files:

    a. First file will contain: starting number of reads, post human depletion of reads, post NG mapping reads, kraken starting and kraken post mapping hits

    b. Second file will contain: NG reference mapping coverage and different depths, assembly metrics from quast, busco results and ST




The workflow generates all the output in the output prefix folder set in the config file (instructions on setup found [below](#config)). Each workflow steps gets its own individual folder as shown. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and some of their contents.**

## Installation 


> If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname. 

```

cd /scratch/esnitkin_root/esnitkin1/your_uniqname/

```

> Clone the github directory onto your system. 

```

git clone https://github.com/Snitkin-Lab-Umich/NG-probe.git

```

> Ensure you have successfully cloned NG-probe. Type `ls` and you should see the newly created directory **_NG-probe_**. Move to the newly created directory.

```

cd NG-probe

```

> Load Bioinformatics, snakemake and singularity modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity mamba

```

This workflow makes use of singularity containers available through [State Public Health Bioinformatics group](https://github.com/StaPH-B/docker-builds). If you are working on Great Lakes (umich cluster)—you can load snakemake and singularity modules as shown above. However, if you are running it on your local or other computing platform, ensure you have snakemake and singularity installed.


## Setup config, samples and cluster files

**_If you are just testing this pipeline, the config and sample files are already loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix (name of your output folder) in the config file to give you an idea of what variables need to be modified when running your own samples on NG-probe._**

### Config
As an input, the snakemake file takes a config file where you can set the path to `samples.csv`, path to your raw sequencing reads, path to adapter fasta file etc. Instructions on how to modify `config/config.yaml` is found in `config.yaml`. 

### Samples
Add samples to `config/samples.csv` following the explanation provided below. `samples.csv` should be a comma seperated file consisting of two columns—`sample_id` and `illumina_r1`.

* `sample_id` is the prefix that should be extracted from your FASTQ reads. For example, in  your raw FASTQ files directory, if you have a file called `Rush_KPC_110_R1.fastq.gz`, your sample_id would be `Rush_KPC_110`.

* `illumina_r1` is the name of the entire raw FASTQ file. In the same directory,  if your file is called `Rush_KPC_110_R1.fastq.gz`, your sample_id would be `Rush_KPC_110_R1.fastq.gz`. **_Only include forward reads._**

You can create samples.csv file using the following for loop. Replace *path_to_your_raw_reads* below with the actual path to your raw sequencing reads.

```

echo "sample_id,illumina_r1" > config/samples.csv

for read1 in path_to_your_raw_reads/*_R1.fastq.gz; do
    sample_id=$(basename $read1 | sed 's/_R1.fastq.gz//g')
    read1_basename=$(basename $read1)
    echo $sample_id,$read1_basename
done >> config/samples.csv

```

### Cluster file

Reduce the walltime in `config/cluster.json` to ensure the jobs are being submitted in a timely manner. 

## Quick start


### Run NG-probe on a set of samples.

> Preview the steps in NG-probe by performing a dryrun of the pipeline. 

```

snakemake -s workflow/Snakefile --dryrun 

```

>Run NG-probe directly on terminal (**_note: if you close your computer/shut down terminal, the pipeline will stop running. Terminal window has to be open until pipeline runs to completion._**)

```

snakemake -s workflow/Snakefile -p --use-conda --use-singularity --use-envmodules -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --conda-frontend mamba --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock

```
> Submit NG-probe as a batch job (**reccommended**)

Change these `SBATCH` commands: `--job-name` to a more descriptive name like run_NG-probe, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. Once you have made the necessary changes, save the below script as `bash_script_to_run_NG-probe.sbat`. Don't forget to submit the job to Slurm! `sbatch bash_script_to_run_NG-probe.sbat`.

```
#!/bin/bash

#SBATCH --job-name=run_NG-probe
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=6-00:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard

# Load necessary modules
module load Bioinformatics snakemake singularity mamba

# Run Snakemake
snakemake -s workflow/Snakefile -p --use-conda --use-singularity --conda-frontend mamba -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 100

```

<!--![Alt text](images/QCD_dag.svg) -->

### Tool stack used in workflow

* [kneaddata](https://github.com/biobakery/kneaddata)
* [clumpify](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/)
* [Bowtie2](https://github.com/BenLangmead/bowtie2)
* [kraken](https://ccb.jhu.edu/software/kraken/)
* [SPades](https://github.com/ablab/spades)
* [skani](https://github.com/bluenote-1577/skani)
* [bioawk](https://github.com/lh3/bioawk)
* [Prokka](https://github.com/tseemann/prokka)
* [mlst](https://github.com/tseemann/mlst)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Pandas](https://pandas.pydata.org/)
* [BUSCO](https://busco.ezlab.org/)
* [AMRFinder](https://github.com/ncbi/amr)
* [QUAST](https://quast.sourceforge.net/)
* [Fastqc](https://github.com/s-andrews/FastQC)