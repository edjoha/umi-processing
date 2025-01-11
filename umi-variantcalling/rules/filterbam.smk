#creates a mutation bed file
rule mutation_bed:
    input:
        "table/{sample}.edit.csv"
    output:
        temp("filterbam/{sample}.bed")
    params:
        config['filter_bam']['padding']
    benchmark:
        "benchmarks/mutation_bed/{sample}.tsv"
    resources:
        runtime="30m"
    run:
        dir_path = config["general"]["work_dir"] + "variantcalling"
        conf = config['filter_bam']
        print(f"{input}")
        anno_df = pd.read_csv(f"{input}", sep='\t').sort_values(['Chr', 'Start']).iloc[:,:5]
        if not len(anno_df.index):
            anno_df.to_csv("filterbam/{sample}.bed", index=False, sep='\t', header=False)
        else:
            # get the bedfile with padded and collapsed regions
            bed_df = anno_df.groupby('Chr').apply(reduce_regions, conf['padding'])
            # remove Chr index
            bed_df = bed_df.reset_index().drop(columns='level_1')
            # write bed_df to file
            print('MUTFILE', anno_df, '\n', 'BED', bed_df)
            bed_df.to_csv(f"{output}", index=False, sep='\t', header=False)
            #print('MUTFILE', anno_df, '\n', 'BED', bed_df)

#filters bam file on bed file mutations
rule samtools_view:
    input:
        lambda wildcards: input_bam[wildcards.sample], "filterbam/{sample}.bed"
    output:
        temp("filterbam/{sample}.bam")
    params:
        extra="-bhL filterbam/{sample}.bed" # optional params string
    benchmark:
        "benchmarks/samtools_view/{sample}.tsv"
    resources:
        runtime="30m"
    wrapper:
        "v5.0.2/bio/samtools/view"

#produces index file
rule samtools_index_vc:
    input:
        "filterbam/{sample}.bam"
    output:
        temp("filterbam/{sample}.bai")
    benchmark:
        "benchmarks/samtools_index/{sample}.tsv"
    threads:
        4
    params:
        "" # optional params string
    resources:
        runtime="30m"
    wrapper:
        "v5.0.2/bio/samtools/index"

#creates "pieces" of bam files to facilitate high dynamic reagion (HDR) detection
rule mpilup:
    input:
        # single or list of bam files
        bam=lambda wildcards: input_bam[wildcards.sample],
        reference_genome=config["reference"]["genome"]
    output:
        temp("mpileup/{sample}.raw.pileup.gz")
    log:
        "logs/samtools/mpileup/{sample}.log"
    benchmark:
        "benchmarks/samtools_mpilup/{sample}.tsv"
    params:
        extra="".join([" -q ", str(config['mpileup']['MAPQ']), " -Q ", str(config['mpileup']['Q'])]),  # optional
    resources:
        runtime="30m"
    wrapper:
        "v5.0.2/bio/samtools/mpileup"

rule filter_pileup:
    input:
        "mpileup/{sample}.raw.pileup.gz"
    output:
        temp("pileup/{sample}.pileup")
    conda:
        "../envs/mawk.yaml"
    threads:
        2
    benchmark:
        "benchmarks/filter_pileup/{sample}.tsv"
    resources:
        runtime="30m"
    params:
        ''.join([config["general"]['snakedir'],'/scripts/shell/cleanpileup.mawk'])
    shell:
        r"""
        mkdir -p pileup
        zcat {input} | {params} > {output}
        """

#rule IGVnav:
#    input:
#        filter = "table/{sample}.filter1.csv",
#        bam = "filterbam/{sample}.bam"
#    output:
#        IGVNav = "filterbam/{sample}.filter1.txt"
#    threads:
#        1
#    benchmark:
#        "benchmarks/IGVnav/{sample}.tsv"
#    resources:
#        runtime="30m"
#    run:
#        # selectinng the right filter2 file from the configs
#        df = pd.read_csv(input.filter, sep='\t', index_col=False).iloc[:, :5]
#        print(f'Loaded {input.filter}')
#        for col in ['Call', 'Tags', 'Notes']:
#            df[col] = ''
#        df.loc[:, 'Chr'] = df['Chr'].str.replace('chr', '')
#        df.to_csv(str(output), sep='\t', index=False)
#        print(f"Written to {output.IGVNav}")
