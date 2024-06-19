import os

# Load configuration
configfile: "config.yaml"

# Define the sample list from the config file
BASE_SAMPLES = config["base"]["samples"]
_keys=list(config["mapping"]["AN"].keys())
_values=list(config["mapping"]["AN"].values())

# Define the final target rule
rule all:
    input:
        # expand("{output}", output=[config["base"]["output"][sample] for sample in BASE_SAMPLES]),
        # expand("results/{sample}/target/QC/target.QC.fam", sample=BASE_SAMPLES),
        # expand("results/{sample}/target/QC/target.QC.snplist", sample=BASE_SAMPLES),
        # expand("results/{sample}/target/QC/target.QC.nodup.snplist", sample=BASE_SAMPLES),
        # expand("results/{sample}/target/QC/target.valid.sample", sample=BASE_SAMPLES),
        # expand("results/{sample}/target/QC/target.a1", sample=BASE_SAMPLES),
        # expand("results/{sample}/target/QC/target.mismatch", sample=BASE_SAMPLES),
        # expand("results/{sample}/target/QC/target.QC.valid", sample=BASE_SAMPLES), # used in sexcheck
        config["results"]["AN"]["final_target_prefix"] + "/target.QC.bed",
        config["results"]["AN"]["final_target_prefix"] + "/target.QC.bim",
        config["results"]["AN"]["final_target_prefix"] + "/target.QC.fam",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.5.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.4.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.3.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.2.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.1.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.05.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.001.profile"


# Rule: Quality Control for base data using R
rule qc_base:
    input:
        file=config["base"]["input"]["AN"]
    output:
        config["base"]["output"]["AN"]
    shell:
        """
        module load r/4.3.2 && Rscript scripts/qc_base.R --input {input.file} --output {output} --keys {_keys} --values {_values}
        """

# target data
# Step 2: Initial QC on target data
rule qc_target:
    input:
        bed=config["target"]["genotype"]["bed"],
        bim=config["target"]["genotype"]["bim"],
        fam=config["target"]["genotype"]["fam"]
    output:
        "results/AN/target/QC/target.QC.fam",
        "results/AN/target/QC/target.QC.snplist"
    params:
        input_prefix=config["target"]["genotype"]["prefix"],
        output_prefix=config["results"]["AN"]["qc_prefix_prefix"]+"/target.QC"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
          --bfile {params.input_prefix} \
          --maf 0.01 \
          --hwe 1e-6 \
          --geno 0.01 \
          --mind 0.01 \
          --write-snplist \
          --make-just-fam \
          --out {params.output_prefix}
        """

# Step 3: Pruning for LD
rule prune_ld:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.QC.fam",
        snplist="results/AN/target/QC/target.QC.snplist"
    output:
        "results/AN/target/QC/target.QC.prune.in",
        "results/AN/target/QC/target.QC.prune.out"
    params:
        input_prefix=config["target"]["genotype"]["prefix"],
        output_prefix=config["results"]["AN"]["qc_prefix_prefix"]+"/target.QC"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
          --bfile {params.input_prefix} \
          --keep {input.fam} \
          --extract {input.snplist} \
          --indep-pairwise 200 50 0.25 \
          --out {params.output_prefix}
        """

# Step 4: Calculate heterozygosity rate
rule calculate_het:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.QC.fam",
        prune_in="results/AN/target/QC/target.QC.prune.in"
    output:
        "results/AN/target/QC/target.QC.het"
    params:
        input_prefix=config["target"]["genotype"]["prefix"],
        output_prefix=config["results"]["AN"]["qc_prefix_prefix"]+"/target.QC"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
          --bfile {params.input_prefix} \
          --extract {input.prune_in} \
          --keep {input.fam} \
          --het \
          --out {params.output_prefix}
        """

rule filter_het:
    input:
        "results/AN/target/QC/target.QC.het"
    output:
        "results/AN/target/QC/target.valid.sample"
    shell:
        """
        Rscript scripts/filter_het.R {input} {output}
        """

# Step 5: Strand flipping and recoding
rule strand_flipping:
    input:
        bim=config["target"]["genotype"]["bim"],
        AN=config["base"]["output"]["AN"],
        snplist="results/AN/target/QC/target.QC.snplist"
    output:
        a1="results/AN/target/QC/target.a1",
        mismatch="results/AN/target/QC/target.mismatch"
    shell:
        """
        mkdir -p results/AN/target/QC && Rscript scripts/strand_flipping.R --bim {input.bim} --target {input.AN} --snp {input.snplist} --output_a1 {output.a1} --output_mismatch {output.mismatch}
        """

# Step 6: De-duplicate rs-ID
rule dedeuplicate:
    input:
        "results/AN/target/QC/target.QC.snplist"     
    output:
        "results/AN/target/QC/target.QC.nodup.snplist"
    shell:
        """
        Rscript scripts/deduplicate_snplist_rs_id.R {input} {output}
        """
# Step 7: Checking sex-related quality issues
# rule check_sex:
#    input:
#        bed=config["target"]["genotype"]["bed"],
#        fam="results/AN/target/QC/target.valid.sample",
#        prune_in="results/AN/target/QC/target.QC.prune.in"
#    output:
#        "results/AN/target/QC/target.QC.sexcheck"
#    params:
#        input_prefix=config["target"]["genotype"]["prefix"],
#        output_prefix=config["results"]["AN"]["qc_prefix_prefix"]+"/target.QC"
#    shell:
#        """
#        /proj/htzhu/UKB_GWAS/phase1and2/plink \
#          --bfile {params.input_prefix} \
#          --extract {input.prune_in} \
#          --keep {input.fam} \
#          --check-sex \
#          --out {params.output_prefix}
#        """

# rule filter_sex:
#    input:
#        sexcheck="results/AN/target/QC/target.QC.sexcheck",
#        valid="results/AN/target/QC/target.valid.sample"
#    output:
#        "results/AN/target/QC/target.QC.valid"
#    shell:
#        """
#        Rscript scripts/filter_sex.R {input.valid} {input.sexcheck} {output}
#        """

# Step 8: Checking relatedness between samples
rule check_relatedness:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.valid.sample", # or target.QC.valid produced by sexcheck
        prune_in="results/AN/target/QC/target.QC.prune.in"
    output:
        "results/AN/target/QC/target.QC.rel.id"
    params:
        input_prefix=config["target"]["genotype"]["prefix"],
        output_prefix=config["results"]["AN"]["qc_prefix_prefix"]+"/target.QC"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
          --bfile {params.input_prefix} \
          --extract {input.prune_in} \
          --keep {input.fam} \
          --rel-cutoff 0.125 \
          --out output_prefix
        """

# Step 8: Generate final QCed target data
rule generate_final_data:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.QC.rel.id",
        snplist="results/AN/target/QC/target.QC.nodup.snplist",
        mismatch="results/AN/target/QC/target.mismatch",
        a1="results/AN/target/QC/target.a1"
    output:
        bed=config["results"]["AN"]["final_target_prefix"] + "/target.QC.bed",
        bim=config["results"]["AN"]["final_target_prefix"] + "/target.QC.bim",
        fam=config["results"]["AN"]["final_target_prefix"]+ "/target.QC.fam"
    params:
        input_prefix=config["target"]["genotype"]["prefix"],
        output_prefix=config["results"]["AN"]["final_target_prefix"]+"/target.QC"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
          --bfile {params.input_prefix} \
          --make-bed \
          --keep {input.fam} \
          --out {params.output_prefix} \
          --extract {input.snplist} \
          --exclude {input.mismatch} \
          --a1-allele {input.a1}
        """

# step 9: update effect size

rule update_effect_size:
    input:
        summary_statistic=config["base"]["output"]["AN"]
    output:
        "results/AN/base/QC/target.QC.transformed"

    shell:
        # """
        # Rscript scripts/update.R --input {input.summary_statistic} --output {output}
        # """
        # already preprocessed
        """
        mkdir -p results/AN/base/QC/ &&
        cp {input.summary_statistic} {output}
        """

# step 10 clumping

rule clumping:
    input:
        base="results/AN/base/QC/target.QC.transformed",
        bed=config["results"]["AN"]["final_target_prefix"]+"/target.QC.bed",
        bim=config["results"]["AN"]["final_target_prefix"]+"/target.QC.bim",
        fam=config["results"]["AN"]["final_target_prefix"]+"/target.QC.fam"
    output:
        clumped_output=config["prs"]["AN"]["output_prefix"]+"/prs.clumped",
        valid_snp=config["prs"]["AN"]["output_prefix"]+"/prs.valid.snp"
    params:
        input_prefix=config["results"]["AN"]["final_target_prefix"]+"/target.QC",
        output_prefix=config["prs"]["AN"]["output_prefix"]+"/prs"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
          --bfile {params.input_prefix} \
          --clump-p1 1 \
          --clump-r2 0.1 \
          --clump-kb 250 \
          --clump {input.base} \
          --clump-snp-field SNP \
          --clump-field P \
          --out {params.output_prefix} && 
        awk 'NR!=1{{print $3}}' {output.clumped_output} > {output.valid_snp} 
        """
    
# step 11 generate valid 
rule generate_PRS_step1:
    input:
        base="results/AN/base/QC/target.QC.transformed"
    output:
        rangel=config["prs"]["AN"]["output_prefix"]+"/range_list",
        snp_value=config["prs"]["AN"]["output_prefix"]+"/SNP.pvalue"
    shell:
        """
        awk '{{print $3,$8}}' {input.base} > {output.snp_value} && # $3 and $8 correspond to SNPID and P
        echo "0.001 0 0.001" > {output.rangel} &&
        echo "0.05 0 0.05" >> {output.rangel} &&
        echo "0.1 0 0.1" >> {output.rangel} &&
        echo "0.2 0 0.2" >> {output.rangel} &&
        echo "0.3 0 0.3" >> {output.rangel} &&
        echo "0.4 0 0.4" >> {output.rangel} &&
        echo "0.5 0 0.5" >> {output.rangel} &&
        """

rule generate_PRS_step2:
    input:
        base="results/AN/base/QC/target.QC.transformed",
        bed=config["results"]["AN"]["final_target_prefix"]+"/target.QC.bed",
        bim=config["results"]["AN"]["final_target_prefix"]+"/target.QC.bim",
        fam=config["results"]["AN"]["final_target_prefix"]+"/target.QC.fam",
        rangel=config["prs"]["AN"]["output_prefix"]+"/range_list",
        snp_value=config["prs"]["AN"]["output_prefix"]+"/SNP.pvalue",
        valid_snp=config["prs"]["AN"]["output_prefix"]+"/prs.valid.snp"
    output:
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.5.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.4.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.3.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.2.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.1.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.05.profile",
        config["prs"]["AN"]["output_prefix"]+"/prs"+".0.001.profile"
    params:
        input_prefix=config["results"]["AN"]["final_target_prefix"]+"/target.QC",
        output_prefix=config["prs"]["AN"]["output_prefix"]+"/prs"
    shell:
        """
        /proj/htzhu/UKB_GWAS/phase1and2/plink \
        --bfile {params.input_prefix} \
        --score {input.base} 3 4 6 header \
        --q-score-range {input.rangel} {input.snp_value} \
        --extract {input.valid_snp} \
        --out {params.output_prefix}
        """

rule population_accounting:
    input:
        bed=config["results"]["AN"]["final_target_prefix"]+"/target.QC.bed",
        bim=config["results"]["AN"]["final_target_prefix"]+"/target.QC.bim",
        fam=config["results"]["AN"]["final_target_prefix"]+"/target.QC.fam"
    output:
        prune_in=config["prs"]["AN"]["output_prefix"]+"/prs.prune.in",
        eigenvec=config["prs"]["AN"]["output_prefix"]+"/prs.eigenvec"
    params:
        input_prefix=config["results"]["AN"]["final_target_prefix"]+"/target.QC",
        output_prefix=config["prs"]["AN"]["output_prefix"]+"/prs"
    shell:
        """
        # First, we need to perform prunning
        plink \
        --bfile {params.input_prefix} \
        --indep-pairwise 200 50 0.25 \
        --out {params.output_prefix}
        &&
        # Then we calculate the first 6 PCs
        plink \
        --bfile EUR.QC \
        --extract {output.prune_in} \
        --pca 6 \
        --out {params.output_prefix}
        """
