# Cross-ancestry analysis identifies genes associated with obesity risk and protection
Codebase for manuscript titled "Cross-ancestry analysis identifies genes associated with obesity risk and protection";

# Repository guidelines
This repo contains step-by-step guide to the analysis carried out for the project. A description of each sub dir is as follows:

1. src/: All biobank-specific preprocessing, variant annotation and association testing scripts are present here.
2. notebooks/: All downstream analysis scripts are present here.

## src/exome_annot/ukb
- 0_MT.ipynb: Converting the pVCF files to Hail Matrix Table format
- 1_variant_annot_vep109.ipynb: Annotating the variants using VEP


## src/phenotype_prep/ukb
- 0_bmi_covariates.ipynb: Acquiring BMI and other covariates from UKB
- 1_ICD10.ipynb: Acquiring ICD diagnosis codes for each individual from UKB
- 2_Lifestyle_factors.ipynb: Acquiring lifestyle factors for each individual from UKB
- 3_prepare_bmi_assoc_pheno_file.ipynb: Parse the phenotypes to create population specific files
- 4_create_step2_helper_files.ipynb: Prepare the genetic and phenotype data to use with REGENIE 

## src/assoc_test/ukb
- 0_qc_snp_data.sh: Quality control filter for genotype files
- 1_run_regenie_step1_quant.sh: REGENIE Step 1
- 2_run_regenie_step2_quant_helper.sh: REGENIE Step 2 helper script called by the main script
- 2_run_regenie_step2.sh: REGENIE Step 2 main script

## src/post_processing/ukb
- 00_create_meta_file.ipynb: REGENIE population specific association test result parsing file
- 01_icd_enrichment.ipynb: Contingency table for obesity related comorbidities in UKB
- 02_obesity_cat_enrichment.ipynb: Contingency table for obesity categories in UKB
- 03_proteomics.ipynb: UKB proteomics data analysis
- 04_bmi_and_pgs_interaction.ipynb: Multiplicative model for BMI PGS interaction in UKB
- 05_pgs_int_plots.ipynb: Interaction plots of genes showing non-additive interaction with PGS
- 06_lf_ancova.ipynb: ANCOVA model of discovered genes interacting with lifestyle to modulate BMI
- 07_bmi_variablity_plot.ipynb: BMI distribution of individuals carrying PTVs in discovered genes

Note: The AoU scripts are very similar to the UKB scripts except certain cohort specific difference in global variables. We intend to make our AoU workspace a community workspace soon to make the scripts accessible to all registered AoU researchers.

## notebooks/
- 0_meta_analysis.ipynb: Aggregating statistics from individual population association study using inverse variance weighted random effect model
- 1_known_obesity_genes.ipynb: Meta statistics of previously identified obesity genes with their effect across ancestries
- 2_icd_enrichment.ipynb: Enrichment of obesity related comorbidity in individuals with PTVs in discovered genes
- 3_bmi_cat_enrichment.ipynb: Enrichment of obesity clinical categories in individuals with PTVs in discovered genes
- 4_functional_enrichment.ipynb: Gene set enrichment analysis of the discovered genes
- 5_make_manual_lit_table.ipynb: Literature review of the discovered genes using APIs
- 6_proteomics_data.ipynb: Forest plot of the protein model coefficients
- 7_pgs_interaction.ipynb: Aggregating statistics of Gene-PGS interaction from UKB and AoU
- 8_drugbank.ipynb: Finding existing drug targets for discovered genes
- 9_cohort_stats.ipynb: Number of individuals per cohort 
- 10_known_genes_overlap.ipynb: Venn diagram of discovered genes found in complementary approaches
- 11_lf_bias.ipynb: Modulation of gene effects based on lifestyle in UKB and AoU
- 12_population_variance.ipynb: Interpopulation variance of the discovered and known obesity genes
- 13_all_associations.ipynb: Upset plot of discovered genes found in complementary approaches
- 14_pubmed_search.ipynb: Pubmed search of discovered genes and relation to obesity
- 15_supplementary_tables.ipynb: Supplementary table creation script
- 16_pgs_interaction_schematic.ipynb: Stick and ball schematic of PGS interaction modes

# Demo dataset
Input data for scripts in src/ folder can be accessed by registering to UK Biobank or All of Us portal. The scripts are parallelized by chromosomes and chrm21 can be used as a demo dataset to test the scripts. The exome annotation script for chmr 21 takes approximately 3 hours to run on a Jupyter Lab spark cluster instance in UK Biobank using 2 nodes with 48 cores each. For association testing using REGENIE, estimates for each step is available here: https://rgcgithub.github.io/regenie/recommendations/. All of Us time estimates are very similar to the UK Biobank.


# System requirements and installation guide
UK Biobank and All of Us both offer linux based cloud compute environment with most packages pre-installed. Python packages used in the source code such as statsmodels can be directly installed through pip. The first cell of each notebook contains a bash magic line to install packages that are not installed in the default jupyter environment provided by the respective cloud compute platform. It takes less than a minute to install these python packages. More details about the cloud compute platforms are available here:

1. UK Biobank: https://dnanexus.gitbook.io/uk-biobank-rap
2. All of Us: https://www.researchallofus.org/about-the-research-hub/



# Acknowledgements
We thank the participants and investigators in the UK Biobank and All of US research studies who made this work possible. This research has been conducted using the UK Biobank Resource under Application Number 45023. We also thank the National Institutes of Health’s All of Us Research Program for making available the participant data examined in this study.
