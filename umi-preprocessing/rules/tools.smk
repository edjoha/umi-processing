#produces reference files
rule get_refs:
    input:
        genome=config["reference"]["genome"],
        known_indels=config["reference"]["known_indels"]
    output:
        genome="refs/genome.fasta",
        known_indels="refs/known_indels.vcf.gz",
        index="refs/known_indels.vcf.gz.tbi"
    conda: "../envs/refs.yaml"
    benchmark:
        "benchmarks/tools/get_refs.tsv"
    resources:
        mem="2G",
        runtime="1h"
    shell:
        r"""
        mkdir -p refs
        cp {input.genome} {output.genome}
        cp {input.known_indels} {output.known_indels}
        tabix -p vcf {output.known_indels}
        """

#produces index file
rule samtools_faidx:
    input:
        "refs/genome.fasta"
    output:
        "refs/genome.fasta.fai"
    params:
        "" # optional params string
    benchmark:
        "benchmarks/tools/samtools_faidx.tsv"
    resources:
        mem="2G",
        runtime="1h"
    wrapper:
        "v5.0.1/bio/samtools/faidx"

#sorts bam by queryname
rule query_bam_sort:
    input:
        "{file}.bam"
    output:
        "{file}_qsorted.bam"
    conda: "../envs/mapping.yaml"
    resources:
        mem="20G",
        runtime="30m"
    benchmark:
        "benchmarks/tools/query_bam_sort/{file}.tsv"
    log:
        "logs/picard/query_bam_sort/{file}.log"
    shell:
        r"""
        picard SortSam I={input} \
        SORT_ORDER=queryname \
        o={output} &> {log}
        """

#sorts bam by coordinate
rule coordinate_bam_sort:
    input:
        "fixed-rg/{sample}.filtered.bam"
    output:
        "mapped/{sample}.filtered_csorted.bam"
    conda: "../envs/mapping.yaml"
    resources:
        mem="20G",
        runtime="30m"
    benchmark:
        "benchmarks/tools/coordinate_bam_sort/{sample}.filtered.tsv"
    log:
        "logs/picard/coordinate_bam_sort/{sample}.filtered.log"
    shell:
        r"""
        picard SortSam I={input} \
        SORT_ORDER=coordinate \
        o={output} &> {log}
        """

#produces index/helper files for bwa_mem2 (map_reads1b and 2b)
rule bwa_mem2_index:
    input:
        "refs/genome.fasta",
    output:
        multiext("refs/genome.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    benchmark:
        "benchmarks/tools/bwa_mem2_index.tsv"
    log:
        "logs/bwa-mem2_index/genome.log",
    resources:
        mem="100G",
        runtime="2h"
    wrapper:
        "v5.5.0/bio/bwa-mem2/index"

#creates dictionary
rule create_dict:
    input:
        "refs/genome.fasta"
    output:
        "refs/genome.dict"
    benchmark:
        "benchmarks/tools/create_dict.tsv"
    log:
        "logs/picard/create_dict.log"
    params:
        extra=""  # optional: extra arguments for picard.
    resources:
        mem_mb=1024,
        runtime="1h"
    wrapper:
        "v5.0.1/bio/picard/createsequencedictionary"

#creates index files
rule samtools_index:
    input:
        "mapped/{file}.bam"
    output:
        "mapped/{file}.bam.bai"
    benchmark:
        "benchmarks/tools/samtools_index/{file}.tsv"
    resources:
        mem="2G",
        runtime="1h"
    wrapper:
        "v3.3.3/bio/samtools/index"

#creates interval list
rule bed_to_interval_list:
    input:
        bed=config["reference"]["region_file"],
        dict="refs/genome.dict"
    output:
        "refs/region.intervals"
    benchmark:
        "benchmarks/tools/bed_to_interval_list.tsv"
    log:
        "logs/picard/bedtointervallist.log"
    params:
        # optional parameters
        extra="--SORT true", # sort output interval list before writing
    resources:
        mem_mb=1024,
        runtime="1h"
    wrapper:
        "v1.3.2/bio/picard/bedtointervallist"

#saves all files to keep (.realigned.bam, .realigned.bai, multiqc.html) in a longterm storage location to avoid being deleted on scratch
rule save:
    input:
        bam=expand("mapped/{sample}.realigned.bam", sample=SAMPLES),
        index=expand("mapped/{sample}.realigned.bam.bai", sample=SAMPLES),
        multiqc="qc/multiqc.html"
    output:
        multiqc=config["longterm_storage"]["multiqc"],
        bam=directory(config["longterm_storage"]["realignedbamdir"]),
        flag=touch("logs/save/preprocess/done.txt")
    benchmark:
        "benchmarks/save/demux.tsv"
    resources:
        runtime="12h"
    log:
        "logs/save/save_realigned.log"
    shell:
        """
        rsync -a {input.bam} {input.index} {output.bam}
        rsync -a {input.multiqc} {output.multiqc} &> {log}
        """
