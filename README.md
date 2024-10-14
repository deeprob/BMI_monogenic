# Cross-ancestry analysis identifies genes associated with obesity risk and protection
Codebase for manuscript titled "Cross-ancestry analysis identifies genes associated with obesity risk and protection";

# Repository guidelines
This repo contains step-by-step guide to the analysis carried out for the project. A description of each sub dir is as follows:

1. src/: All preprocessing, variant annotation and association testing scripts are present here.
2. notebooks/: All downstream analysis scripts are present here.

## src/exome_annot

### ukb
- 0_MT.ipynb: Converting the pVCF files to Hail Matrix Table format
- 1_variant_annot_vep109.ipynb: Annotating the variants using VEP

### aou


## src/phenotype_prep

### ukb
- 0_bmi_covariates.ipynb: Acquiring BMI and other covariates from UKB
- 1_ICD10.ipynb: Acquiring ICD diagnosis codes for each individual from UKB
- 2_Lifestyle_factors.ipynb: Acquiring lifestyle factors for each individual from UKB
- 3_prepare_bmi_assoc_pheno_file.ipynb: Parse the phenotypes to create population specific files
- 4_create_step2_helper_files.ipynb: Prepare the genetic and phenotype data to use with REGENIE 

### aou

## src/assoc_test

### ukb
- 0_qc_snp_data.sh: Quality control filter for genotype files
- 1_run_regenie_step1_quant.sh: REGENIE Step 1
- 2_run_regenie_step2_quant_helper.sh: REGENIE Step 2 helper script called by the main script
- 2_run_regenie_step2.sh: REGENIE Step 2 main script


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