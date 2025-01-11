rule basecalls_to_fastq:
    input:
        basecalls=config["illumina"]["basecall_dir"],
        SampleSheet=config["general"]["SampleSheet"]
    output:
        flag = touch("logs/demux/done.txt")
    conda: "../envs/bcl2fastq.yaml"
    params:
        rstructure=config["illumina"]["readstructure"],
    threads: workflow.cores * 0.75 #was 30 before ... maybe change back if problems occur
    benchmark:
        "benchmarks/basecalls_to_fastq/basecalls.tsv"
    resources:
        mem=get_mem_70_10,
        runtime="24h"
    shell:
        """
        bcl2fastq \
                    --input-dir {input.basecalls} \
                    --runfolder-dir {input.basecalls}/../../../ \
                    --output-dir reads \
                    --sample-sheet {input.SampleSheet} \
                    --barcode-mismatches 1 \
                    --loading-threads 6 \
                    --processing-threads 30 \
                    --writing-threads 6 \
                    --ignore-missing-bcl \
                    --ignore-missing-filter \
                    --mask-short-adapter-reads 0 \
                    --create-fastq-for-index-reads \
                    --use-bases-mask {params.rstructure}
        """

def get_fq(wildcards):
    # code that returns a list of fastq files  based on *wildcards.sample* e.g.
    return sorted(glob.glob("reads/" + wildcards.sample + '_S*.fastq.gz'))

rule fastq_to_bam:
    input:
        flag = "logs/demux/done.txt",
        fastq=get_fq
    output:
        "unmapped/{sample}.unmapped.bam"
    conda: "../envs/fgbio.yaml"
    params:
        sampleName = "{sample}",
        readstructure = config["fgbio"]["readstructure"]
    log:
        "logs/fgbio/{sample}.log"
    benchmark:
        "benchmarks/fastq_to_bam/{sample}.tsv"
    resources:
        mem=get_mem_70_10
        runtime="5h" #ist das nötig? auch dynamic? runtime=lambda wc, input: max(0.3 * input.size_mb) --> wie komme ich von size zu time (Einheit)?
    shell:
        """
        fgbio FastqToBam \
        --input {input.fastq} \
        --read-structures {params.readstructure} \
        --output {output} \
        --sort true \
        --sample {params.sampleName} \
        --library=illumina &> {log}
        """

rule save:
    input:
        bam=expand("unmapped/{sample}.unmapped.bam", sample=SAMPLES)
    output:
        dir=directory(config["longterm_storage"]["unmappedbamdir"]),
        flag=touch("logs/save/demux/done.txt")
    benchmark:
        "benchmarks/save/demux.tsv"
    resources:
        runtime="4h" #partition short
    log:
        "logs/save/save_unmapped.log"
    shell:
        """
        rsync -a {input.bam} {output.dir} &> {log}
        """
