#!/usr/bin/env python

"""
align_genome.snakefile

Snakemake workflow for downloading genome from NCBI and
RNA-Seq files from SRA
"""

__author__ = "David Levy-Booth"
__copyright__ = "Copyright 2019, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
from pathlib import Path
THREADS = config["threads"]
TRIMMER = config["TRIMMER"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]
if "salmon" in METHOD and len(METHOD) == 1:
    count_out = "RREP42/results/tables/salmon.{{trimmer}}.counts.tsv"
if "salmon" in METHOD and len(METHOD) > 1:
    METHOD.remove("salmon")
    count_out = ["RREP42/results/tables/salmon.{trimmer}.counts.tsv",
    "RREP42/results/tables/{method}.{aligner}.{trimmer}.counts.tsv"]


rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        # expand("genome/{genome_id}/{genome_id}_genomic_prokka.gbff",
        # genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1]),
        expand("RREP4/genome/{genome_id}/{genome_id}_genomic.fna",
        genome_id=config["genomes"][config["genome_id"]]["url"].split("/")[-1]),
        expand(count_out,
        method=METHOD, aligner=ALIGNER, trimmer=TRIMMER)


rule download_genome:
    params:
        genome_url = config["genomes"][config["genome_id"]]["url"]
    output:
        gff = "RREP4/genome/{genome_id}/{genome_id}_genomic.gff",
        gbff = "RREP4/genome/{genome_id}/{genome_id}_genomic.gbff",
        cd_fna = "RREP4/genome/{genome_id}/{genome_id}_cds_from_genomic.fna",
        prep_gbff = "RREP4/genome/{genome_id}/{genome_id}_genomic_prokka.gbff",
        prep_gff = "RREP4/genome/{genome_id}/{genome_id}_genomic_salmon.gff",
        prep_cd_fna = "RREP4/genome/{genome_id}/{genome_id}_cds_from_genomic_salmon.fna",
        fna = "RREP4/genome/{genome_id}/{genome_id}_genomic.fna",
    log:
        "RREP4/log/genome/{genome_id}.log"
    benchmark:
        "RREP4/benchmarks/{genome_id}.download.benchmark.txt"
    run:
        shell("rsync --copy-links --recursive --times --verbose "
        "rsync://{params.genome_url} genome --log-file={log}")
        shell("gunzip -r genome")
        shell("sed 's/product/protein/g' {output.gbff} > {output.prep_gbff}")
        shell("sed -i 's/locus_tag/product/g' {output.prep_gbff}")
        shell("sed 's/locus_tag/gene_id/g' {output.gff} > {output.prep_gff}")
        shell("sed 's/.*\[locus_tag=\([^]]*\)\].*/>\1/g' {output.cd_fna} > {output.prep_cd_fna}")
        touch("{output.fna}")

#sed 's/.*\[locus_tag=\([^]]*\)\].*/>\1/g' /home/david/RREP4/genome/GCF_003004765.2_ASM300476v2/GCF_003004765.2_ASM300476v2_cds_from_genomic.fna > /home/david/RREP4/genome/GCF_003004765.2_ASM300476v2/GCF_003004765.2_ASM300476v2_cds_from_genomic_salmon.fna


rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    max_reads: Maximal number of reads to download for each sample.
    """
    output:
        "RREP4/transcriptome/reads/{sra_id}_1.fastq",
        "RREP4/transcriptome/reads/{sra_id}_2.fastq"
    params:
        max_reads = config["max_reads"]
    threads: THREADS
    benchmark:
        "RREP4/benchmarks/{sra_id}.download.benchmark.txt"
    shell:
        """
        fastq-dump {wildcards.sra_id} -O transcriptome/reads --split-files
        # This clears a cache where SRA Tools reserves a lot of space
        cache-mgr --clear >/dev/null 2>&1
        """

rule rename_reads:
    input:
        "RREP4/transcriptome/reads/{sra_id}_1.fastq",
        "RREP4/transcriptome/reads/{sra_id}_2.fastq"
    output:
        "RREP4/transcriptome/reads/{sra_id}_1.untrimmed.fastq",
        "RREP4/transcriptome/reads/{sra_id}_2.untrimmed.fastq"
    run:
        shell("rename 's/.fastq/.untrimmed.fastq/' {input}")

rule trimmomatic:
    input:
        read1 = "RREP4/transcriptome/reads/{sra_id}_1.untrimmed.fastq",
        read2 = "RREP4/transcriptome/reads/{sra_id}_2.untrimmed.fastq"
    output:
       paired1 = "RREP4/transcriptome/reads/{sra_id}_1.trimmomatic.fastq",
       paired2 = "RREP4/transcriptome/reads/{sra_id}_2.trimmomatic.fastq",
       unpaired1 = "RREP4/transcriptome/reads/{sra_id}_1.unpaired.fastq",
       unpaired2 = "RREP4/transcriptome/reads/{sra_id}_2.unpaired.fastq"
    threads: THREADS
    log:
        "RREP4/log/{sra_id}.trimmomatic.log"
    benchmark:
        "RREP4/benchmarks/{sra_id}.trimmomatic.benchmark.txt"
    shell:
      "trimmomatic PE -threads {threads} -phred33 "
      "{input.read1} {input.read2} {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} "
      "ILLUMINACLIP:RREP4/utils/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> {log}"

rule bbduk:
    input:
        read1 = "RREP4/transcriptome/reads/{sra_id}_1.untrimmed.fastq",
        read2 = "RREP4/transcriptome/reads/{sra_id}_2.untrimmed.fastq"
    output:
       paired1 = "RREP4/transcriptome/reads/{sra_id}_1.bbduk.fastq",
       paired2 = "RREP4/transcriptome/reads/{sra_id}_2.bbduk.fastq"
    threads: THREADS
    log:
        "RREP4/log/{sra_id}.bbduk.log"
    benchmark:
        "RREP4/benchmarks/{sra_id}.bbduk.benchmark.txt"
    shell:
        """
        bbduk.sh -Xmx20g t={threads} in1={input.read1} in2={input.read2} out1={output.paired1} out2={output.paired2} \
        ref=RREP4/utils/adapters/TruSeq3-PE-2.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 2> {log}
        """



rule bt2_index_genome:
    """
    Index a genome using Bowtie 2.
    """
    params:
        genome_base = config["genomes"][config["genome_id"]]["url"].split("/")[-1]
    input:
        "RREP4/genome/{params.genome_base}/{params.genome_base}_genomic.fna"
    output:
        expand("RREP4/intermediate/genome.{n}.bt2", n = ["1","2","3","4"]),
        expand("RREP4/intermediate/genome.rev.{n}.bt2", n = ["1","2"])
    benchmark:
        "RREP4/benchmarks/{params.genome_base}.bt2.index.benchmark.txt"
    shell:
        """
        bowtie2-build {input} RREP4/intermediate/genome 2> {log}
        """

rule bt2_align:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    input:
        fastq_1 = "RREP4/transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        fastq_2 = "RREP4/transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
        index = expand("RREP4/intermediate/genome.{n}.bt2", n = ["1","2","3","4"]),
        index_rev = expand("RREP4/intermediate/genome.rev.{n}.bt2", n = ["1","2"])
    output:
        sam = temp("RREP4/intermediate/{sra_id,\w+}.{trimmer,\w+}.bt2.sam"),
        bam = temp("RREP4/intermediate/{sra_id,\w+}.{trimmer,\w+}.bt2.bam")
    threads: THREADS
    log:
        "RREP4/log/{sra_id}.{trimmer}.bt2.align.log"
    benchmark:
        "RREP4/benchmarks/{sra_id}.{trimmer}.bt2.align.benchmark.txt"
    run:
        # This gives the base name for the genome index, i.e. "intermediate/some_id"
        # rather than "intermediate/some_id.*.bt2"
        indexBase = "RREP4/intermediate/genome"
        shell("bowtie2 -x " + indexBase + " --threads {threads} -1 {input.fastq_1} -2 {input.fastq_2} -S {output.sam} 2> {log}")
        shell("samtools view -bS {output.sam} > {output.bam}")

rule bbmap_align:
    params:
        genome_base = config["genomes"][config["genome_id"]]["url"].split("/")[-1]
    input:
        fastq_1 = "RREP4/transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        fastq_2 = "RREP4/transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
        genome = expand("RREP4/genome/{genome_base}/{genome_base}_genomic.fna", genome_base=config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    output:
        sam = temp("RREP4/intermediate/{sra_id,\w+}.{trimmer,\w+}.bbmap.sam"),
        bam = temp("RREP4/intermediate/{sra_id,\w+}.{trimmer,\w+}.bbmap.bam")
    threads: THREADS
    log:
        "RREP4/log/{sra_id}.{trimmer}.bbmap.align.log"
    benchmark:
        "RREP4/benchmarks/{sra_id}.{trimmer}.bbmap.align.benchmark.txt"
    run:
        shell("bbmap.sh -Xmx20g threads={threads} in1={input.fastq_1} in2={input.fastq_2} \
        out={output.sam} ref={input.genome}  path=genome/bbmap_index/ 2> {log}")
        shell("samtools view -bS {output.sam} > {output.bam}")


rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        "{prefix}.bam"
    output:
        "{prefix}.sorted.bam"
    threads: THREADS
    shell:
        """
        samtools sort {input} > {output}
        """

rule htseq_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams=expand("RREP4/intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam",
        sra_id = config["sample_ids"], trimmer = config["TRIMMER"], aligner = config["ALIGNER"]),
        gff=expand("RREP4/genome/{base}/{base}_genomic.gff", base = config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    output:
        expand("RREP42/results/tables/htseq.{aligner}.{trimmer}.counts.tsv",
        trimmer = config["TRIMMER"], aligner = config["ALIGNER"])
    shadow: "minimal"
    threads: THREADS
    benchmark:
        "RREP42/benchmarks/htseq.{aligner}.{trimmer}.benchmark.txt".format(trimmer=TRIMMER, aligner = ALIGNER)
    shell:
        """
        # Save the count table as a temporary file and then prepend a header line
        # with the sample names
        htseq-count --format bam --type gene --idattr locus_tag -r pos {input.bams} {input.gff} > {output}
        #echo '{input.bams}' | tr ' ' '\t' | cat - tempfile2 > {output}
        """

rule feature_counts_table:
    """
    Generate a count table with featureCounts
    """
    input:
        bams=expand("RREP4/intermediate/{sra_id}.{trimmer}.{aligner}.sorted.bam",
        sra_id = config["sample_ids"], trimmer = config["TRIMMER"], aligner = config["ALIGNER"]),
        gff=expand("RREP4/genome/{base}/{base}_genomic_salmon.gff", base = config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    output:
        expand("RREP42/results/tables/featureCounts.{aligner}.{trimmer}.counts.tsv",
        trimmer = config["TRIMMER"], aligner = config["ALIGNER"])
    threads: THREADS
    log:
        "RREP42/log/featureCounts.{aligner}.{trimmer}.log".format(trimmer=TRIMMER, aligner = ALIGNER)
    benchmark:
        "RREP42/benchmarks/featureCounts.{aligner}.{trimmer}.benchmark.txt".format(trimmer=TRIMMER, aligner = ALIGNER)
    shell:
        """
        featureCounts -O -p -T {threads} -t gene -F GTF -a {input.gff} -o {output} {input.bams}
        """

rule salmon_index:
    input:
        fasta=expand("RREP4/genome/{base}/{base}_cds_from_genomic_salmon.fna", base = config["genomes"][config["genome_id"]]["url"].split("/")[-1])
    output:
        directory("RREP42/genome/salmon_quasi")
    threads: 1
    run:
        shell("salmon index -t {input.fasta} -i {output} --type quasi -k 31")


rule salmon_quant:
    """
    Generate directories containing count files with salmon (quasi mode)
    """
    input:
        fastq_1="RREP4/transcriptome/reads/{sra_id}_1.{trimmer}.fastq",
        fastq_2="RREP4/transcriptome/reads/{sra_id}_2.{trimmer}.fastq",
        salmon_dir = directory("genome/salmon_quasi")
    output:
        directory("RREP42/transcriptome/salmon/{sra_id}_{trimmer}")
    threads: THREADS
    log:
        "RREP42/log/salmon.{sra_id}.{trimmer}.log"
    benchmark:
        "RREP42/benchmarks/salmon.quant.{sra_id}.{trimmer}.benchmark.txt"
    run:
        shell("salmon quant -i {input.salmon_dir} -l A -p {threads} --validateMappings \
        -1 {input.fastq_1} -2 {input.fastq_2} -o {output} 2> {log}")


rule salmon_quant_table:
    input:
        expand(directory("RREP42/transcriptome/salmon/{sra_id}_{trimmer}"), sra_id = config["sample_ids"], trimmer = TRIMMER)
    output:
        "RREP42/results/tables/salmon.{trimmer}.counts.tsv"
    run:
        shell("salmon quantmerge --quants {input} --column numreads  -o {output}")
