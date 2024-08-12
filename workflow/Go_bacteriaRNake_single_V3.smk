# 20230922 (Heekuk Park)
# delete annotation.txt. only using gff

import pandas as pd
import glob

# Global variables
master_dir = config["project"] + "_RNAseq_output"
READ_DIR = config["read_dir"] # Get the directory of fastq.gz files from the command line
conda_shotgun = "/media/uhlemann/core4/DB/MAGs/shotgun.yaml"
conda_MAGs = "/media/uhlemann/core4/DB/MAGs/MAGs.yaml"


READS = glob.glob(READ_DIR + '/*R1*.fastq.gz')

SAMPLES = [read.split('/')[-1].split('_R1_001')[0] for read in READS]

shell("mkdir -p {master_dir}")

# Make required directories
directories = [master_dir, master_dir + "/1_trim",  master_dir + "/2_bowtie2_index", master_dir + "/3_bowtie2_files", master_dir + "/4_htseq-count"]

for dir in directories:
    shell("mkdir -p {dir}")

rule all:
    input:
        expand(f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt", sample=SAMPLES),
        f"{master_dir}/2_bowtie2_index/index_build.done",
        f"{master_dir}/merged_counts.csv",
        f"{master_dir}/merged_counts_with_gene_names.csv"  # Add this line



rule trim_reads:
    input:
        f"{READ_DIR}/{{sample}}_R1_001.fastq.gz"
    output:
        f"{master_dir}/1_trim/{{sample}}.R1.output.fastq.gz"
    conda:
        conda_shotgun
    shell:
        """
        echo "Trimming {wildcards.sample} ..."
        trimmomatic SE -threads 10 -summary {master_dir}/1_trim/{wildcards.sample}.statssummary -phred33 {input} {output} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
        """

rule index_genome:
    output:
        touch(f"{master_dir}/2_bowtie2_index/index_build.done")
    params:
        genome=config["genome"]
    conda:
        conda_MAGs
    shell:
        """
        bowtie2-build --threads 6 -f {params.genome} {master_dir}/2_bowtie2_index/index
        touch {output}  # Create a dummy file after bowtie2-build is done.
        """



rule map_reads:
    input:
        reads=f"{master_dir}/1_trim/{{sample}}.R1.output.fastq.gz",
        index=f"{master_dir}/2_bowtie2_index/index_build.done"  # Use the dummy file as an input
    output:
        f"{master_dir}/3_bowtie2_files/{{sample}}.aligned.sam"
    conda:
        conda_MAGs
    shell:
        """
        echo "Mapping {wildcards.sample} against indexed NR5452 WT SPADES assembly..."
        # Provide the prefix to the bowtie2 command
        bowtie2 --no-unal -p 12 -x {master_dir}/2_bowtie2_index/index -U {input.reads} -S {output} 2>{master_dir}/3_bowtie2_files/{wildcards.sample}.log
        """


rule sort_sam:
    input:
        f"{master_dir}/3_bowtie2_files/{{sample}}.aligned.sam"
    output:
        f"{master_dir}/3_bowtie2_files/{{sample}}.sorted.sam"
    conda:
        conda_MAGs
    shell:
        """
        echo "Processing {wildcards.sample} ..."
        samtools sort -n {input} -o {output}
        """

rule count_htseq:
    input:
        sam=f"{master_dir}/3_bowtie2_files/{{sample}}.sorted.sam"
    output:
        f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt"
    params:
        reference=config["gff"] # Get the reference filename from command line argument
    conda:
        conda_shotgun
    shell:
        """
        echo "Processing {wildcards.sample} ..."
        htseq-count --order=name --stranded=no --type=gene --idattr=ID -a 8 -o {master_dir}/4_htseq-count/{wildcards.sample}.htseq.sam {input.sam} {params.reference} > {output}
        """


rule merge_counts:
    input:
        files=expand(f"{master_dir}/4_htseq-count/{{sample}}.gene_id.minqual8.txt", sample=SAMPLES)
    output:
        f"{master_dir}/merged_counts.csv"
    run:
        print(f"Merging counts. Input files: {input.files}, Output file: {output}")

        # Create a dictionary where each key is a filename (without extension) and each value is a Series of counts
        dfs = {}
        for file in input.files:
            sample_name = file.split('/')[-1].split('.')[0] # Change this line to correctly extract sample name
            df = pd.read_csv(file, sep='\t', index_col=0, header=None, names=[sample_name])
            dfs[sample_name] = df[sample_name]

        # Concatenate all Series along the column axis into a DataFrame
        merged_df = pd.concat(dfs, axis=1)

        # Fill any missing values with 0
        merged_df.fillna(0, inplace=True)

        # Write the DataFrame to a new CSV file
        merged_df.to_csv(output[0])

rule map_gene_ids:
    input:
        count_file=master_dir + "/merged_counts.csv",
        annotation_file=config["gff"]  # Using GFF for gene annotation
    output:
        mapped_file=master_dir + "/merged_counts_with_gene_names.csv"
    shell:
        """
        python -c "
import pandas as pd
import numpy as np

# Read the count file and the GFF file
count_df = pd.read_csv('{input.count_file}', index_col=0)
annotation_df = pd.read_csv('{input.annotation_file}', sep='\\t', comment='#', header=None)

# Extract the gene names and locus tags
annotation_df = annotation_df[annotation_df[2] == 'gene']
annotation_df['locus_tag'] = annotation_df[8].str.extract('ID=([^;]+)')
annotation_df['symbol'] = annotation_df[8].str.extract('Name=([^;]+)')


name_to_locus = annotation_df[['symbol', 'locus_tag']].drop_duplicates()

# Remove 'gene-' prefix from index
# count_df.index = count_df.index.str.replace('gene-', '')

# Reset index to move gene ID to a column
count_df.reset_index(inplace=True)
count_df.columns.values[0] = 'locus_tag'

# Map locus tags in count file to gene names
count_df['symbol'] = count_df['locus_tag'].map(name_to_locus.set_index('locus_tag')['symbol'])

# If the gene name is missing, fill it with the locus tag
count_df['symbol'] = count_df['symbol'].where(pd.notnull(count_df['symbol']), count_df['locus_tag'])

# Rearrange the columns to place 'symbol' next to 'locus_tag'
cols = count_df.columns.tolist()
rearranged_cols = [cols[0]] + [cols[-1]] + cols[1:-1]
count_df = count_df[rearranged_cols]

# Save the result to a new CSV file
count_df.to_csv('{output.mapped_file}', index=False)"
        """
