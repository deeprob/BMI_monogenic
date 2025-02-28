# Discovery of novel obesity genes through cross-ancestry analysis
Codebase for manuscript titled "Discovery of novel obesity genes through cross-ancestry analysis" 

[DOI](https://doi.org/10.1101/2024.10.13.24315422) | [Supplementary Tables](https://www.medrxiv.org/content/medrxiv/early/2025/02/28/2024.10.13.24315422/DC2/embed/media-2.xlsx?download=true) | [Supplementary Note](https://www.medrxiv.org/content/medrxiv/early/2025/02/28/2024.10.13.24315422/DC1/embed/media-1.docx?download=true])

# Repository guidelines
This repo contains step-by-step guide to the analysis carried out for the project. A description of each sub dir is as follows:

## notebooks/
All plotting and supplementary material preparation scripts are present here.

- 0_all_gene_forestplot.ipynb: Forest plot with cross ancestry results for BMI-associated genes
- 1_per_ancestry_forestplot.ipynb: Forest plot with per ancestry results of BMI associated genes
- 2_study_specific_qq_plots.ipynb: QQ plots of independent gene burden tests
- 3_metal_statistics.ipynb: Cross-ancestry results using METAL with genomic inflation control
- 4_bmi_cat_enrichment.ipynb: Enrichment of obesity clinical categories in individuals with PTVs in discovered genes
- 5_comorbidity_enrichment.ipynb: Enrichment of obesity related comorbidity in individuals with PTVs in discovered genes
- 6_phewas_manhattan.ipynb: Manhattan plots of YLPM1 and GIGYF1
- 7_pgs_interaction.ipynb: PGS interaction forest plot
- 8_proteomics.ipynb: Proteomics forest plot
- 9_known_gene_forestplot.ipynb: Forest plot with cross-ancestry results of known genes
- 10_known_gene_sex_strat_forestplot.ipynb: Forest plot with cross-ancestry results of sex-specific genes
- 11_supplementary_tables.ipynb: Supplementary tables organized and generated

## src/
All biobank-specific preprocessing, data preparation, association testing and meta analysis scripts are present here.

- genetic_data_processing/
    - snp_qc/
        - 0_qc_snp_data.sh: Variant and sample qc for SNP array data
        - 1_qc_snp_data_ld_prune.sh: LD pruning of the qc'd array variants
        - 2_snp_anc_helper.sh: Helper script for filtering variants with at least 5 samples per ancestry
        - 2_snp_anc.sh: Filtering SNP array variants with at least 5 samples per ancestry
    - variant_qc:
        - run_vqc.sh: Calls the notebook implementing variant quality control
        - dnanexus_notebooks/0_initial_variant_qc_chr1.ipynb: Converting the pVCF files to Hail Matrix Table format and implementing genotype and variant quality control filters for chromosome 1
    - sample_qc/
        - high_quality_variants/
            - run_autosome_hqc.sh: Get a set of high quality variants in the autosome
            - run_autosome_hqc_pruned.sh: LD prune the set of high quality variants
            - merge_autosome_hqc.sh: Merge high quality variants per chromosome into a single file
            - merge_autosome_hqc_pruned.sh: Merge LD pruned high quality variants per chromosome into a single file
            - run_X_hqc.sh: Get a set of high quality variants for X chromosome
            - filter_X_hqc.sh: Filter pseudoautosomal regions in X chromosome
            - run_X_hqc_pruned.sh: LD prune high quality variants in X chromosome
        - relatedness/
            - create_autosome_subsets.sh: Break samples with high quality autosomal variants into multiple subsets
            - get_king_estimates_diff_subset.sh: Estimate kinship coefficients of samples in different subsets
            - get_king_estimates_same_subset.sh: Estimate kinship coefficients of samples in same subsets
            - get_king_estimates.sh: Slower script to estimate kinship coefficients of samples without subsetting
            - dnanexus_notebooks/kinship.ipynb: Organize kinship estimates to infer related individuals
        - impute_sex/
            - dnanexus_notebooks/impute_sex.ipynb: Infer sex of individuals from high quality X chromosome variants
        - flag_samples/
            - run_repartition.sh: Calls the notebook that implements autosomal qc'd variant repartitioning
            - run_sample_qc.sh: Calls the notebook implementing sample qc
            - dnanexus_notebooks/chr1.ipynb: Repartition the variant qc'd matrix tables into smaller number of partitions for easier computing
            - dnanexus_notebooks/sample_annot.ipynb: Calculates PCs of unrelated and relateed individuals using high quality autosomal variants
            - dnanexus_notebooks/sample_qc.ipynb: Calculates ancestry residualized genomic qc metrics 
            - dnanexus_notebooks/flag_samples.ipynb: Marks all samples that failed any sample qc metrics
    - variant_annot/
        - run_vannot.sh: Calls the variant annotation script 
        - dnanexus_notebooks/chr1.ipynb: Annotates rare, deleterious variants using vep and dbnsfp and collects qc passed samples that possess those variants
    - burden_preparation/
        - dnanexus_notebooks/0_prepare_burden_file.ipynb: Categorized the variants by their transcript consequence annotations
        - dnanexus_notebooks/1_add_gnomad_annot.ipynb: Adds gnomad observed population specific MAF to each variant
        - dnanexus_notebooks/2_create_assoc_files.ipynb: Creates association files with variant consequence, maximum minor allele frequency and gene annotation for gene burden tests

- ancestry_inference/
    - dnanexus_notebooks/sample_qc.ipynb: Removes previously flagged samples that failed quality control
    - dnanexus_notebooks/pca_computation.ipynb: Calculates PC of 1000 Genome samples and projects the UKB samples in the same PC space
    - dnanexus_notebooks/pca_computation_ni.ipynb: Copy of the above file running the final steps
    - dnanexus_notebooks/ancestry_inference.ipynb: Trains a random forest classifier on the 1000 Genome samples and predicts ancestry of UKB samples

- phenotypic_data_processing/
    - dnanexus_notebooks/05a_extract_bmi_covariate_data.ipynb: Acquiring BMI and other covariates from UKB
    - dnanexus_notebooks/05b_extract_icd_codes_data.ipynb: Acquiring ICD diagnosis codes for each individual from UKB
    - dnanexus_notebooks/05c_extract_lifestyle_data.ipynb: Acquiring lifestyle factors for each individual from UKB
    - dnanexus_notebooks/05d_prepare_pheno_files.ipynb: Prepare phenotype file with all relevant information
    - dnanexus_notebooks/05e_normalize_pheno.ipynb: Normalize BMI, add additional covariates and create ancestry specific phenotype files

- gene_burden_association_tests/
    - 1_run_regenie_step1_quant_helper.sh: Regenie step 1 helper script per ancestry
    - 1_run_regenie_step1.sh: Regenie step 1 run
    - 2_run_regenie_step2_quant_helper.sh: Regenie step 2 helper script per ancestry
    - 2_run_regenie_step2.sh: Regenie step 2 run

- downstream/
    - 0_process_meta_files.py: Process meta analysis output genes to hgnc updated symbols
    - 1_conduct_meta_analysis.py: Inverse variance weighted fixed effect meta analysis script
    - 2_create_meta_files.py: Creates supplementary tables for meta analysis results
    - 3_enrichment_cat.py: BMI category and comorbidity enrichment script
    - 4_sem.py: Mediation meta analysis script
    - 5_phewas.py: Phewas meta analysis script
    - download_hgnc.sh: Download the latest hgnc symbols for each gene
    - run_proteomics_protein_assoc.sh: Calls script to find differentially expressed proteins
    - dnanexus_notebooks/10a_shadow_effect_nearby_common_variants.ipynb: Finds common variants near BMI-associated genes and overlaps the list with GWAS hits for BMI
    - dnanexus_notebooks/11a_medication_extract_data.ipynb: Extracts medication data for UKB samples
    - dnanexus_notebooks/11b_medication_prepare_pheno.ipynb: Prepare phenotype file with medication as a covariate
    - dnanexus_notebooks/12a_enrichment_obesity_category.ipynb: Create contingency table of PTV carriers in BMI clinical category
    - dnanexus_notebooks/12b_enrichment_obesity_comorbidity: Creates contingency table fo PTV carriers carrying obesity related comorbidities
    - dnanexus_notebooks/12c_phewas_map_samples_to_phecodes.ipynb: Maps ICD codes of UKB samples to phecodes using R PheWAS package
    - dnanexus_notebooks/12d_phewas_run.ipynb: Runs phewas of PTV carriers per ancestry
    - dnanexus_notebooks/12e_sem_conduct.ipynb: Conducts mediation analysis of BMI-associated genes using lavaan
    - dnanexus_notebooks/13b_pgs_interaction.ipynb: Create PTV-PGS interaction models for BMI-associated genes
    - dnanexus_notebooks/14a_proteomics_protein_assoc.ipynb: Finds differentially expressed proteins in PTV carriers of BMI associated genes
    - dnanexus_notebooks/14b_proteomics_bmi_assoc.ipynb: Create association models of differentially expressed genes with BMI

**Note**: The AoU scripts are very similar to the UKB scripts except certain cohort specific difference in global variables. We intend to make our AoU workspace a community workspace soon to make the scripts accessible to all registered AoU researchers.


# Protocol description
A detailed description of the protocols can be obtained using in the [Supplementary Note](https://www.medrxiv.org/content/medrxiv/early/2025/02/28/2024.10.13.24315422/DC1/embed/media-1.docx?download=true]) of the manuscript: 

# Demo dataset
Input data for scripts in src/ folder can be accessed by registering to UK Biobank or All of Us portal. The scripts are parallelized by chromosomes and chrm 21 can be used as a demo dataset to test the scripts. The time estimates along with the cluster configuration of variant quality control, annotation and gene association test scripts are provided in [Supplementary Note](https://www.medrxiv.org/content/medrxiv/early/2025/02/28/2024.10.13.24315422/DC1/embed/media-1.docx?download=true]) of the manuscript.


# System requirements and installation guide
UK Biobank and All of Us both offer linux based cloud compute environment with most packages pre-installed. Python packages used in the source code which are not preinstalled, such as statsmodels, can be directly installed through pip. The first cell of each notebook contains a bash magic line to install packages that are not installed in the default jupyter environment provided by the respective cloud compute platform. It takes less than a minute to install these python packages. More details about the cloud compute platforms along with the version number of their installed packages is available here:

1. UK Biobank: https://dnanexus.gitbook.io/uk-biobank-rap
2. All of Us: https://www.researchallofus.org/about-the-research-hub/

UK Biobank exome annotation, association testing and downstream analysis were performed in Jupyter Lab with Spark Cluster v2.2.0, Swiss Army Knife v4.12.0 and Jupyter Lab v2.3.0 respectively. All of Us data analysis was performed with Controlled Tier Dataset v8.

# Acknowledgements
We thank the participants and investigators in the UK Biobank and All of US research studies who made this work possible. This research has been conducted using the UK Biobank Resource under Application Number 45023 and All of Us Controlled Tier Dataset v8. We also thank the National Institutes of Health’s All of Us Research Program for making available the participant data examined in this study.
