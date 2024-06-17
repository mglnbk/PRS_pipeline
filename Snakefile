import os

# Load configuration
configfile: "config.yaml"

# Define the sample list from the config file
BASE_SAMPLES = config["base"]["samples"]

# Define the final target rule
rule all:
    input:
        expand("{output}", output=[config["base"]["output"][sample] for sample in BASE_SAMPLES]),
        expand("results/{sample}/target/QC/target.QC.fam", sample=BASE_SAMPLES),
        expand("results/{sample}/target/QC/target.QC.snplist", sample=BASE_SAMPLES),
        expand("results/{sample}/target/QC/target.QC.nodup.snplist", sample=BASE_SAMPLES),
        expand("results/{sample}/target/QC/target.valid.sample", sample=BASE_SAMPLES),
        expand("results/{sample}/target/QC/target.a1", sample=BASE_SAMPLES),
        expand("results/{sample}/target/QC/target.mismatch", sample=BASE_SAMPLES),
        expand("results/{sample}/target/QC/target.QC.valid", sample=BASE_SAMPLES),
        expand(config["target"]["final_prefix_prefix"] + "/target.QC.bed", sample=BASE_SAMPLES),
        expand(config["target"]["final_prefix_prefix"] + "/target.QC.bim", sample=BASE_SAMPLES),
        expand(config["target"]["final_prefix_prefix"] + "/target.QC.fam", sample=BASE_SAMPLES)

# Rule: Quality Control for base data using R
rule qc_base:
    input:
        config["base"]["input"]["AN"]
    output:
        config["base"]["output"]["AN"]
    shell:
        """
        module load r/4.3.2 && Rscript scripts/qc_base.R --input {input} --output {output}
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
    shell:
        """
        plink \
          --bfile {input.bed[:-4]} \
          --maf 0.01 \
          --hwe 1e-6 \
          --geno 0.01 \
          --mind 0.01 \
          --write-snplist \
          --make-just-fam \
          --out results/AN/target/QC/target.QC
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
    shell:
        """
        plink \
          --bfile {input.bed} \
          --keep {input.fam} \
          --extract {input.snplist} \
          --indep-pairwise 200 50 0.25 \
          --out results/AN/target/QC/target.QC
        """


# Step 4: Calculate heterozygosity rate
rule calculate_het:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.QC.fam",
        prune_in="results/AN/target/QC/target.QC.prune.in"
    output:
        "results/AN/target/QC/target.QC.het"
    shell:
        """
        plink \
          --bfile {input.bed} \
          --extract {input.prune_in} \
          --keep {input.fam} \
          --het \
          --out results/AN/target/QC/target.QC
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
        "results/AN/target/QC/target.a1",
        "results/AN/target/QC/target.mismatch"
    shell:
        """
        mkdir -p results/AN/target/QC && \
        Rscript scripts/strand_flipping.R {input.bim} {input.height} {input.snplist}
        """

# Step 6: De-duplicate rs-ID
rule dedeuplicate:
    input:
        "results/AN/target/QC/target.QC.snplist"     
    output:
        "results/AN/target/QC/target.QC.nodup.snplist"
    shell:
        """
        Rscript scripts/dedeuplicate_snplist_rs_id.R {input} {output}
        """

# Step 7: Checking sex-related quality issues
rule check_sex:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.valid.sample",
        prune_in="results/AN/target/QC/target.QC.prune.in"
    output:
        "results/AN/target/QC/target.QC.sexcheck"
    shell:
        """
        plink \
          --bfile {input.bed} \
          --extract {input.prune_in} \
          --keep {input.fam} \
          --check-sex \
          --out results/AN/target/QC/target.QC
        """

rule filter_sex:
    input:
        sexcheck="results/AN/target/QC/target.QC.sexcheck",
        valid="results/AN/target/QC/target.valid.sample"
    output:
        "results/AN/target/QC/target.QC.valid"
    shell:
        """
        Rscript scripts/filter_sex.R {input.valid} {input.sexcheck} {output}
        """

# Step 8: Checking relatedness between samples
rule check_relatedness:
    input:
        bed=config["target"]["genotype"]["bed"],
        fam="results/AN/target/QC/target.QC.valid",
        prune_in="results/AN/target/QC/target.QC.prune.in"
    output:
        "results/AN/target/QC/target.QC.rel.id"
    shell:
        """
        plink \
          --bfile {input.bed} \
          --extract {input.prune_in} \
          --keep {input.fam} \
          --rel-cutoff 0.125 \
          --out results/AN/target/QC/target.QC
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
        bed=config["target"]["final_prefix_prefix"] + "/target.QC.bed",
        bim=config["target"]["final_prefix_prefix"] + "/target.QC.bim",
        fam=config["target"]["final_prefix_prefix"] + "/target.QC.fam"
    shell:
        """
        plink \
          --bfile {input.bed} \
          --make-bed \
          --keep {input.fam} \
          --out {output.bed[:-4]} \
          --extract {input.snplist} \
          --exclude {input.mismatch} \
          --a1-allele {input.a1}
        """



