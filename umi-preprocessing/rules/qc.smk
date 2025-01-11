#creates fastqc files to get fastq statistics (fastqc)
rule get_fastqc_input:
    input:
        "unmapped/{sample}.unmapped.bam"
    output:
        read1=temp("reads/{sample}.R1.fastq"),
        read2=temp("reads/{sample}.R2.fastq")
    conda:
            "../envs/mapping.yaml"
    params:
        musage=config["picard"]["memoryusage"]
    benchmark:
        "benchmarks/get_fastqc_input/{sample}.tsv"
    log:
        "logs/reads/{sample}.samtofastq.log"
    threads:
        8
    resources:
        mem_mb=get_mem_40_10,
        runtime="1h"
    shell:
        r"""
        picard {params.musage} SamToFastq I={input} F={output.read1} SECOND_END_FASTQ={output.read2} &> {log}
        """

#produces fastq statistics/quality metrics
rule fastqc:
    input:
        "reads/{sample}.{read}.fastq"
    output:
        html=temp("qc/fastqc/{sample}.{read}.html"),
        zip=temp("qc/fastqc/{sample}.{read}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    benchmark:
        "benchmarks/fastqc/{sample}.{read}.tsv"
    log:
        "logs/fastqc/{sample}.{read}.log"
    threads: 1
    resources:
        mem_mb=1024,
        runtime="3h"
    wrapper:
        "v4.7.2/bio/fastqc"

#produces comprehensive statistics from alignment file
rule samtools_stats:
    input:
        "mapped/{sample}.{type}.bam"
    params:
        extra=" ".join(["-t",config["reference"]["region_file"]]),
        region=""
    output:
        temp("qc/samtools-stats/{sample}.{type}.txt")
    resources:
        mem="10G",
        runtime="1h"
    benchmark:
        "benchmarks/samtools_stats/{sample}.{type}.tsv"
    log:
        "logs/samtools-stats/{sample}.{type}.log"
    wrapper:
        "v3.3.3/bio/samtools/stats"

#Collects hybrid-selection (HS) metrics for a SAM or BAM file.
#This tool takes a SAM/BAM file input and collects metrics that are specific for sequence datasets generated through hybrid-selection
#Hybrid-selection (HS) is the most commonly used technique to capture exon-specific sequences for targeted sequencing experiments such as exome sequencing
rule picard_collect_hs_metrics:
    input:
        bam="mapped/{sample}.{type}.bam",
        reference="refs/genome.fasta",
        # Baits and targets should be given as interval lists. These can
        # be generated from bed files using picard BedToIntervalList.
        bait_intervals="refs/region.intervals",
        target_intervals="refs/region.intervals"
    output:
        temp("qc/hs_metrics/{sample}.{type}.txt")
    params:
        # Optional extra arguments. Here we reduce sample size
        # to reduce the runtime in our unit test.
        extra="--SAMPLE_SIZE 1000"
    benchmark:
        "benchmarks/picard_collect_hs_metrics/{sample}.{type}.tsv"
    log:
        "logs/picard/collect_hs_metrics/{sample}.{type}.log"
    resources:
        mem_mb=1024,
        runtime="2h"
    wrapper:
        "v5.5.0/bio/picard/collecthsmetrics"

#combines all statistics in a single report
rule multiqc:
    input:
        expand("qc/fastqc/{sample}.{reads}_fastqc.zip", sample=SAMPLES, reads=["R1","R2"]),
        expand("qc/{ctype}/{sample}.{ftype}.txt", sample=SAMPLES, ctype=["samtools-stats","hs_metrics"], ftype=["woconsensus", "realigned"])
    output:
        temp("qc/multiqc.html")
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    log:
        "logs/multiqc/multiqc.log"
    params:
        "--interactive --force --cl_config 'max_table_rows: 10000'"
    resources:
        mem="40G",
        runtime="2h"
    wrapper:
        "v3.7.0-10-g491d5b6/bio/multiqc"
