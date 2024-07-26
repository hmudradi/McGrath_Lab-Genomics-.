# Genome Assembly Pipeline with Snakemake

This repository contains a genome assembly pipeline using Snakemake. The pipeline includes quality control of reads, assembly with multiple tools (SKESA, SPAdes), quality assessment with REAPR, and reconciliation of assemblies with GAM-NGS.

## Prerequisites

Before running the pipeline, ensure the following tools are installed and available in your `PATH`:

- fastp
- SKESA
- SPAdes
- bwa
- samtools
- reapr
- GAM-NGS

You will also need to have Snakemake installed. You can install it using conda:

```
conda install -c bioconda snakemake
```

## Getting Started
### Clone the repository:

```
git clone https://github.com/yourusername/genome-assembly-pipeline.git
cd genome-assembly-pipeline
```

### Configure the pipeline:
Edit the config.yaml file to specify the paths to your input reads, output directory, and other parameters:

```
input_reads: "/path/to/input_reads"
output_folder: "/path/to/output"
tools_folder: "/path/to/tools"
threads: 6
memory: 10
cut_mean_quality: 28
average_qual: 28
cut_window_size: 20
kmer: 21
cores: 4
mem: 10
Bmin: 10
Tc: 0.75
isolates:
  - isolate1
  - isolate2
  - isolate3
```

### Run the pipeline:

```
snakemake --cores 4
```

Or, if you are using a SLURM cluster, use the following command:

```
snakemake --jobs 100 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -n {cluster.ntasks}" --cluster-config cluster.yaml
```

### Monitor the pipeline:

You can visualize the pipeline's progress using the following command:

```
snakemake --dag | dot -Tsvg > dag.svg
```

### Pipeline Overview
>Quality Control (fastp): This step involves filtering and trimming the reads based on quality scores.

>Assembly (SKESA, SPAdes): These steps perform genome assembly using three different tools.

>BAM Generation: Aligns reads back to assemblies and generates sorted BAM files.

>Quality Assessment (REAPR): This step assesses the quality of the assemblies.

>Assembly Reconciliation (GAM-NGS): Reconciles multiple assemblies to create a consensus meta-assembly.
