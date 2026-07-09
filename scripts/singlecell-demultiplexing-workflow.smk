import os
import pandas as pd

configfile: "config.yaml"

inputpath = config["BASE_DIR"]
outputpath = config["OUT_DIR"]
metadata_file = config["METADATA_FILE"]
singlecellfile = config["SINGLECELL_PATH_PREFIX"]
actacseqfile = config["SINGLECELLATACSEQ_PATH_PREFIX"]
metadata = pd.read_csv(metadata_file, sep="\t")
#callrates = ["0.50", "0.30", "0.80"]
callrates = ["0.50"]


rule all:
    input:
        outputpath + "/data/merged_vcf/renamed.vcf.gz",
        subset_vcfs,
        filter_subset_vcfs,
        Maf_filter_subset_vcfs,
	cellsnp_outputs,
	viero_outputs,
        cellsnp_atacseq_outputs,
        viero_atacseq_outputs


rule Subset_Pool:
    input:
        vcf = outputpath + "/data/merged_vcf/renamed.vcf.gz",
        samples = outputpath + "/temp_samples/pool{pool}_samples.txt"
    output:
        vcf = protected(outputpath + "/data/pools/pool{pool}.vcf.gz"),
        modulindex = protected(outputpath + "/data/pools/pool{pool}.vcf.gz.tbi")
    params:
        module_load = config["modules"]["bcftools_module"]
    shell:
        """
        module load {params.module_load}
        echo "[INFO] Pool {wildcards.pool}: samples ="
        cat {input.samples}
        echo "[INFO] Running bcftools..."
        bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf} || echo "bcftools view failed"
        bcftools index --tbi -f {output.vcf}
        ls -lh {output.vcf} {output.modulindex}
        """


rule Filter_Subset_Pool:
    input:
        vcf = outputpath + "/data/pools/pool{pool}.vcf.gz"
    output:
        vcf = protected(outputpath + "/data/filter_pools/pool{pool}.vcf.gz"),
        index = protected(outputpath + "/data/filter_pools/pool{pool}.vcf.gz.tbi"),
        sortedvcf = protected(outputpath + "/data/filter_pools/pool{pool}.sorted.vcf.gz"),
        sortedindex = protected(outputpath + "/data/filter_pools/pool{pool}.sorted.vcf.gz.tbi")
    params:
        module_load = config["modules"]["bcftools_module"],
	module_load2 = config["modules"]["vcfsort_module"],
        script = "workflow/scripts/Filter_MergedVCF2.pl",
	chromorderfile=config["CHROMORDER"],
        tmp_out = temp(outputpath + "/data/filter_pools/pool{pool}.vcf"),
        tmp_out2 = temp(outputpath + "/data/filter_pools/pool{pool}.sorted.vcf")
    shell:
        """
        module load {params.module_load}

        # Run Perl script to produce uncompressed VCF
        perl {params.script} {input.vcf} {params.tmp_out}

        # Compress the result
        bgzip -c {params.tmp_out} > {output.vcf}

        # Index it
        bcftools index --tbi -f {output.vcf}
	
	module purge
	{params.module_load2}

	vcf-sort -p 22  {output.vcf}  >{params.tmp_out2}
	bgzip  -c {params.tmp_out2} > {output.sortedvcf}
	tabix -p vcf  {output.sortedvcf}


        #bcftools index --tbi -f {output.sortedvcf}

        # Cleanup if needed
        rm -f {params.tmp_out}
        rm -f {params.tmp_out2}

        ls -lh {output.vcf} {output.index} {output.sortedvcf} {output.sortedindex}
        """


rule Filter_MAF:
    input:
        vcf = outputpath + "/data/filter_pools/pool{pool}.sorted.vcf.gz"
    output:
        vcf = protected(outputpath + "/data/filter_pools_MAF/pool{pool}.sorted.vcf.gz"),
    params:
        module_load = config["modules"]["bcftools_module"]
    shell:
        """
        module load {params.module_load}
        bcftools view -m2 -M2 -v snps {input.vcf} -Ou | \
        bcftools +fill-tags -- -t MAF | \
        bcftools +fill-tags --threads 10 -Oz -o {output.vcf} -- -t AC,AN,AF 
        tabix -p vcf {output.vcf}
        """


rule cellsnp:
    input:
        vcf = outputpath + "/data/filter_pools_MAF/pool{pool}.sorted.vcf.gz",
        bam = lambda wildcards: f"{singlecellfile}/{pool_to_cellranger_dir[wildcards.pool]}/outs/possorted_genome_bam.bam",
        barcodes = lambda wildcards: f"{singlecellfile}/{pool_to_cellranger_dir[wildcards.pool]}/outs/filtered_feature_bc_matrix/barcodes.tsv.g
z"
    output:
        out = protected(f"{outputpath}/data/cellsnp/Pool{{pool}}/cellSNP.base.vcf.gz"),
        out2 = directory(f"{outputpath}/data/cellsnp/Pool{{pool}}")
    params:
        module_load = config['modules']['cellsnp_module'],
    threads: 22
    shell:
        """
        module purge

        {params.module_load}

        mkdir -p {outputpath}/data/cellsnp/


        cellsnp-lite -s {input.bam} \
            -R {input.vcf} \
            -b {input.barcodes}  \
            --minMAF 0.05  --minCOUNT 10 -p 22 --gzip \
            -O {output.out2} 
        """



rule viero:
    input:
        vcf = outputpath + "/data/filter_pools_MAF/pool{pool}.sorted.vcf.gz",
        cell = outputpath + "/data/cellsnp/Pool{pool}/cellSNP.base.vcf.gz",
        samples = outputpath + "/temp_samples/pool{pool}_samples.txt",
        inf = f"{outputpath}/data/cellsnp/Pool{{pool}}"
    output:
        summary = f"{outputpath}/data/viero/Pool{{pool}}/summary.tsv",
        out2 = directory(f"{outputpath}/data/viero/Pool{{pool}}")
    params:
        module_load = config['modules']['viero_module'],
    threads: 22
    shell:
        """
        module purge

        {params.module_load}


        NCOUNT=$(wc -l < {input.samples})
        echo "Running vireo on pool {wildcards.pool} with N=$NCOUNT"

        vireo \
          -c {input.inf} \
          -N $NCOUNT  \
          -p {threads} \
          -d {input.vcf} \
          -o {output.out2}

        """


rule atacseq_vcf:
    input:
        vcf    = outputpath + "/data/filter_pools_MAF/pool{pool}.sorted.vcf.gz",
        chrmap = outputpath + "/temp_samples/to_chr.map"
    output:
        vcf = protected(outputpath + "/data/filter_pools_MAF_atacseq/pool{pool}.sorted.vcf.gz"),
        tbi = protected(outputpath + "/data/filter_pools_MAF_atacseq/pool{pool}.sorted.vcf.gz.tbi"),
        vcf2 = protected(outputpath + "/data/filter_pools_MAF_atacseq/pool{pool}.sorted2.vcf.gz"),
        tbi2 = protected(outputpath + "/data/filter_pools_MAF_atacseq/pool{pool}.sorted2.vcf.gz.tbi")
    threads: 8
    params:
        module_load = config["modules"]["bcftools_module"]
    log:
        outputpath + "/logs/rename_chr_pool{pool}.log"
    shell:
        r"""
        module purge
        module load {params.module_load}
        mkdir -p $(dirname {output.vcf})
        bcftools annotate --rename-chrs {input.chrmap}  -Oz --threads {threads} -o {output.vcf2} {input.vcf} 2> {log}
        tabix -f -p vcf {output.vcf2}
        bcftools sort -Oz -o {output.vcf} {output.vcf2}
        tabix -f -p vcf {output.vcf}

        """


rule actaseq_cellsnp:
    input:
        vcf = outputpath + "/data/filter_pools_MAF_atacseq/pool{pool}.sorted.vcf.gz",
        bam = lambda wildcards: f"{actacseqfile}/{atacseq_cellranger_dir[wildcards.pool]}/outs/possorted_bam.bam",
        barcodes = lambda wildcards: f"{actacseqfile}/{atacseq_cellranger_dir[wildcards.pool]}/outs/filtered_peak_bc_matrix/barcodes.tsv"
    output:
        out = protected(f"{outputpath}/data/cellsnp_atacseq/Pool{{pool}}/cellSNP.base.vcf.gz"),
        out2 = directory(f"{outputpath}/data/cellsnp_atacseq/Pool{{pool}}")
    params:
        module_load = config['modules']['cellsnp_module'],
    threads: 22
    shell:
        """
        module purge

        {params.module_load}

        mkdir -p {outputpath}/data/cellsnp/


        cellsnp-lite -s {input.bam} \
            -R {input.vcf} \
            -b {input.barcodes} --cellTAG CB --UMItag None \
            --minMAF 0.05  --minCOUNT 10 -p 22 --gzip \
            -O {output.out2} 
        """



rule actaseq_viero:
    input:
        vcf = outputpath + "/data/filter_pools_MAF_atacseq/pool{pool}.sorted.vcf.gz",
        cell = outputpath + "/data/cellsnp_atacseq/Pool{pool}/cellSNP.base.vcf.gz",
        samples = outputpath + "/temp_samples/pool{pool}_samples.txt",
        inf = f"{outputpath}/data/cellsnp_atacseq/Pool{{pool}}"
    output:
        summary = f"{outputpath}/data/viero_atacseq/Pool{{pool}}/summary.tsv",
        out2 = directory(f"{outputpath}/data/viero_atacseq/Pool{{pool}}")
    params:
        module_load = config['modules']['viero_module'],
    threads: 22
    shell:
        """
        module purge

        {params.module_load}


        NCOUNT=$(wc -l < {input.samples})
        echo "Running vireo on pool {wildcards.pool} with N=$NCOUNT"

        vireo \
          -c {input.inf} \
          -N $NCOUNT  \
          -p {threads} \
          -d {input.vcf} \
          -o {output.out2}

        """

       

