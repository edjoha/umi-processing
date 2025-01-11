#annotates the variants
rule annovar:
    input:
        "table/{sample}.csv"
    output:
        anno = temp("table/{sample}.anno.csv"),
        tmp = temp("table/{sample}.annotmp.csv")
    benchmark:
        "benchmarks/annovar/{sample}.tsv"
    log:
        "logs/annovar/{sample}.log"
    threads:
        8
    resources:
        runtime="30m",
        mem="20G"
    params:
        annovar = config['tools']['annovar'],
        anno = get_anno_params,
        sd = config["general"]['snakedir'],
        anno_info = f"scripts/shell/sed.cmd"
    shell:
        r"""
        tail -n +2 {input} > {output.tmp} && \
        {params.annovar}/table_annovar.pl {output.tmp} {params.anno} -thread {threads} && \
        sed --file={params.sd}/{params.anno_info} {output.tmp}.hg38_multianno.txt > {output.anno} && \
        rm {output.tmp}.hg38_multianno.txt
        """

#ebfilter: is useful when a "normal reference" is available (e.g. germline reference vs. cancer biopsy)
if config["general"]["control"]:
    rule ebfilter:
        input:
            sample = lambda wildcards: input_bam[wildcards.sample],
            vcf = "vardict/{sample}.vcf"
        output:
            vcf = temp("vardict/{sample}_EB.vcf"),
            txt = temp("vardict/{sample}_EB.txt")
        log:
            "logs/EBFilter/{sample}.log"
        threads:
            4
        resources:
            time=get_time_3_1,
            mem=get_mem_30_10
        benchmark:
            "benchmarks/ebfilter/{sample}.tsv"
        conda:
            "../envs/EBFilter-env.yaml"
        params:
            normals = config['edit']['normals']
        shell:
            r"""
            EBFilter -f vcf -t {threads} {input.vcf} {input.sample} {params.normals} {output.vcf}
            bcftools query -f '[%EB]\n' {output.vcf} > {output.txt} 2>/dev/null
            """
else:
    rule fake_ebfilter:
        input:
            vcf = "vardict/{sample}.vcf"
        output:
            txt = temp("vardict/{sample}_EB.txt")
        log:
            "logs/EBFilter/{sample}.log"
        threads:
            1
        resources:
            time=get_time_1_1
        benchmark:
            "benchmarks/fake_ebfilter/{sample}.tsv"
        run:
           with open(input.vcf, "r") as input_file:
              vcf_lines = sum(1 for line in input_file if not line.startswith("#"))
           
           with open(output.txt, "w") as output_file:
              output_file.write("\n".join(["NaN"] * vcf_lines))

rule add_ebfilter:
    input:
        anno = "table/{sample}.anno.csv",
        ebfilter = "vardict/{sample}_EB.txt"
    output:
        temp("table/{sample}.edit.csv")
    conda:
        "../envs/Renv.yaml"
    benchmark:
        "benchmarks/add_ebfilter/{sample}.tsv"
    threads:
        1
    resources:
        runtime="30m"
    params:
        # wd = config["general"]["work_dir"] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/AddParameters.R",
        candidate = config['edit']['candidate_list'],
        driver = config['edit']['driver_list'],
        CHIP = config['edit']['CHIP_list']
    shell:
        r"""
        Rscript {params.sd}/{params.rs} {input.anno} {output} {params.candidate} {params.driver} {params.CHIP} {input.ebfilter}
        """
        #Rscript {params.sd}/{params.rs} {params.wd}/{input.anno} {params.wd}/{output} {params.candidate} {params.driver} {params.CHIP} {params.wd}/{input.ebfilter}

rule primer3:
    input: "table/{sample}.edit.csv"
    output: temp("table/{sample}.edit.primer.csv")
    conda:
        "../envs/primer3-env.yaml"
    threads: 1
    resources:
        runtime="30m",
        mem_mb=get_mem_3_1
    benchmark:
        "benchmarks/primer3/{sample}.tsv"
    params:
        genome_split = config["primer3"]["split"]
    script:
        "../scripts/primer3.py"


rule detect_HDR:
    input:
        filter_file = "table/{sample}.edit.csv",
        bam = "filterbam/{sample}.bam",
        index = "filterbam/{sample}.bai",
        pileup = "pileup/{sample}.pileup"
    output:
        HDR = temp("table/{sample}.edit.HDR.csv")
    benchmark:
        "benchmarks/detect_HDR/{sample}.tsv"
    conda:
        f"../envs/HDR-env.yaml"
    threads: 1
    resources:
        time=get_time_3_1,
        mem=get_mem_40_10,
        mem_mb=get_mem_40_10
    params:
        min_sim = config['HDR']['min_similarity'],
        min_q = config['HDR']['min_q'],
        min_HDR_count = config['HDR']['min_HDR_count']
    script:
        "../scripts/HDR.py"


rule combine_annotations:
    input:
        anno_edit = "table/{sample}.edit.csv",
        primer = "table/{sample}.edit.primer.csv",
        hdr = "table/{sample}.edit.HDR.csv"
    output: 
        temp("filter/{sample}.edit.csv")
    benchmark:
        "benchmarks/combine_annotations/{sample}.tsv"
    conda:
        "../envs/Renv.yaml"
    params: 
        wd = config["general"]["work_dir"] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/CombineAnno.R"
    threads: 1
    resources:
        runtime="1h"
    shell:
        r"""
        Rscript {params.sd}/{params.rs} {input.anno_edit} {input.primer} {input.hdr} {output}
        """

#combines all annotates variants of all samples within one csv file
rule combine_samples:
    input:
        expand("filter/{sample}.edit.csv", sample = samples)
    output: 
        update("filter/variantcalls.csv")
    conda:
        "../envs/Renv.yaml"
    benchmark:
        "benchmarks/combine_samples.tsv"
    params: 
        wd = config["general"]["work_dir"] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/MergeCalls.R"
    threads: 1
    resources:
        runtime="1h",
        mem_mb=get_mem_30_10
    shell:
        r"""
        Rscript {params.sd}/{params.rs} filter {output}
        """
        #Rscript {params.sd}/{params.rs} {params.wd}/filter {params.wd}/{output}

#saves all files to keep (variantcalls.csv, .vcf.gz) in a longterm storage location to avoid being deleted on scratch
rule save:
    input:
        csv="filter/variantcalls.csv",
        vcf=expand("vardict/{sample}.vcf.gz", sample=samples),
    output:
        csv=config["longterm_storage"]["csv"],
        vcf=directory(config["longterm_storage"]["vcf"]),
        flag=touch("logs/save/varcall/done.txt")
    benchmark:
        "benchmarks/save/varcall.tsv"
    resources:
        runtime="12h"
    log:
        "logs/save/save_variantcalls.log"
    shell:
        """
        rsync -a {input.csv} {output.csv}
        rsync -a {input.vcf} {output.vcf} &> {log}
        """
