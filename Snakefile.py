import os

configfile: "config.yaml"

# Get input reads
READS = config["input_reads"]
OUTPUT_DIR = config["output_folder"]
TOOLS_DIR = config["tools_folder"]

# Rule to initialize directories
rule init_dirs:
    output:
        OUTPUT_DIR,
        TOOLS_DIR,
        os.path.join(OUTPUT_DIR, "read_QC", "fastp"),
        os.path.join(OUTPUT_DIR, "assemblies", "SKESA", "contigs"),
        os.path.join(OUTPUT_DIR, "assemblies", "SPAdes", "contigs"),
        os.path.join(OUTPUT_DIR, "assemblies", "SPAdes", "extra"),
        os.path.join(OUTPUT_DIR, "bam_files"),
        os.path.join(OUTPUT_DIR, "REAPR", "unreconciled"),
        os.path.join(OUTPUT_DIR, "assembly-reconciliation", "tmp")
    shell:
        """
        mkdir -p {output}
        """

# Rule to preprocess reads
rule preprocess_reads:
    input:
        r1=lambda wildcards: os.path.join(READS, wildcards.isolate, f"{wildcards.isolate}_1.fq.gz"),
        r2=lambda wildcards: os.path.join(READS, wildcards.isolate, f"{wildcards.isolate}_2.fq.gz")
    output:
        r1_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_1_fp.fq.gz"),
        r2_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_2_fp.fq.gz")
    params:
        cut_mean_quality=config["cut_mean_quality"],
        average_qual=config["average_qual"],
        cut_window_size=config["cut_window_size"]
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1_fp} -O {output.r2_fp} -M {params.cut_mean_quality} -W {params.cut_window_size} -e {params.average_qual} -c
        """



# Rule to assemble with SKESA
rule assemble_SKESA:
    input:
        r1_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_1_fp.fq.gz"),
        r2_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_2_fp.fq.gz")
    output:
        contig=lambda wildcards: os.path.join(OUTPUT_DIR, "assemblies", "SKESA", "contigs", f"{wildcards.isolate}_SKESA.fasta")
    params:
        cores=config["cores"],
        mem=config["mem"]
    shell:
        """
        skesa --reads {input.r1_fp},{input.r2_fp} --cores {params.cores} --memory {params.mem} > {output}
        """

# Rule to assemble with SPAdes
rule assemble_SPAdes:
    input:
        r1_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_1_fp.fq.gz"),
        r2_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_2_fp.fq.gz")
    output:
        contig=lambda wildcards: os.path.join(OUTPUT_DIR, "assemblies", "SPAdes", "contigs", f"{wildcards.isolate}_SPAdes.fasta")
    params:
        cores=config["cores"]
    shell:
        """
        spades.py -1 {input.r1_fp} -2 {input.r2_fp} -o {output} -t {params.cores}
        """

# Rule to index BAM files
rule index_bam:
    input:
        bam=lambda wildcards: os.path.join(OUTPUT_DIR, "bam_files", f"{wildcards.isolate}_sorted.{wildcards.assembly}.bam")
    output:
        bam_index=lambda wildcards: f"{input.bam}.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Rule to generate BAM files
rule generate_bam:
    input:
        contig=lambda wildcards: os.path.join(OUTPUT_DIR, "assemblies", wildcards.assembly, "contigs", f"{wildcards.isolate}_{wildcards.assembly}.fasta"),
        r1_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_1_fp.fq.gz"),
        r2_fp=lambda wildcards: os.path.join(OUTPUT_DIR, "read_QC", "fastp", wildcards.isolate, f"{wildcards.isolate}_2_fp.fq.gz")
    output:
        bam=lambda wildcards: os.path.join(OUTPUT_DIR, "bam_files", f"{wildcards.isolate}_sorted.{wildcards.assembly}.bam")
    shell:
        """
        bwa index {input.contig}
        bwa mem {input.contig} {input.r1_fp} {input.r2_fp} | samtools fixmate -O bam - - | samtools sort -O bam -o {output.bam}
        """

# Rule to run REAPR
rule run_reapr:
    input:
        bam=lambda wildcards: os.path.join(OUTPUT_DIR, "bam_files", f"{wildcards.isolate}_sorted.{wildcards.assembly}.bam"),
        contig=lambda wildcards: os.path.join(OUTPUT_DIR, "assemblies", wildcards.assembly, "contigs", f"{wildcards.isolate}_{wildcards.assembly}.fasta")
    output:
        qa_report=lambda wildcards: os.path.join(OUTPUT_DIR, "REAPR", "unreconciled", f"QA_{wildcards.assembly}", f"{wildcards.isolate}_QA", "05.summary.report.txt")
    shell:
        """
        reapr pipeline {input.contig} {input.bam} {output.qa_report}
        """

# Rule to reconcile assemblies using GAM-NGS
rule reconcile_assemblies:
    input:
        contigs_SKESA=lambda wildcards: os.path.join(OUTPUT_DIR, "assemblies", "SKESA", "contigs", f"{wildcards.isolate}_SKESA.fasta"),
        contigs_SPAdes=lambda wildcards: os.path.join(OUTPUT_DIR, "assemblies", "SPAdes", "contigs", f"{wildcards.isolate}_SPAdes.fasta"),
        bamlist=lambda wildcards: os.path.join(OUTPUT_DIR, "assembly-reconciliation", f"{wildcards.isolate}.bamlist")
    output:
        meta_assembly=lambda wildcards: os.path.join(OUTPUT_DIR, "assembly-reconciliation", f"{wildcards.isolate}_meta.gam.fasta")
    params:
        Bmin=config["Bmin"],
        Tc=config["Tc"]
    shell:
        """
        gam-create --master-bam {input.bamlist} --slave-bam {input.bamlist} --min-block-size {params.Bmin} --output {output}
        gam-merge --master-bam {input.bamlist} --slave-bam {input.bamlist} --blocks-file {output}.blocks --slave-fasta {input.contigs_SKESA} --output {output} --coverage-filter {params.Tc}
        """

# Define input and output for the workflow
rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "read_QC", "fastp", "{isolate}", "{isolate}_1_fp.fq.gz"), isolate=config["isolates"]),
        expand(os.path.join(OUTPUT_DIR, "read_QC", "fastp", "{isolate}", "{isolate}_2_fp.fq.gz"), isolate=config["isolates"]),
        expand(os.path.join(OUTPUT_DIR, "assemblies", "SKESA", "contigs", "{isolate}_SKESA.fasta"), isolate=config["isolates"]),
        expand(os.path.join(OUTPUT_DIR, "assemblies", "SPAdes", "contigs", "{isolate}_SPAdes.fasta"), isolate=config["isolates"]),
        expand(os.path.join(OUTPUT_DIR, "bam_files", "{isolate}_sorted.{assembly}.bam"), isolate=config["isolates"], assembly=["SKESA", "SPAdes"]),
        expand(os.path.join(OUTPUT_DIR, "REAPR", "unreconciled", "QA_{assembly}", "{isolate}_QA", "05.summary.report.txt"), isolate=config["isolates"], assembly=["SKESA", "SPAdes"]),
        expand(os.path.join(OUTPUT_DIR, "assembly-reconciliation", "{isolate}_meta.gam.fasta"), isolate=config["isolates"])
