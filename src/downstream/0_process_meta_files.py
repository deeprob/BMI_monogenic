import os
import pandas as pd
import re
import numpy as np
from utils import meta

PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"


if __name__=="__main__":
    hgnc_file = os.path.join(PROJECT_DIR, "data/hgnc/protein-coding_gene.txt")
    hgnc_df = pd.read_csv(
        hgnc_file, sep="\t", 
        usecols=["hgnc_id", "symbol", "name", "alias_symbol", "alias_name", "prev_symbol", "prev_name"]
    )
    analyses = ["all_ancestry", "lovo", "meds", "shadow_effect", "male", "female"]
    meta_classes = [meta.Meta, meta.Lovo, meta.Meds, meta.Meds, meta.Meta, meta.Meta]
    phenotypes =[("bmi_rint", "bmi"), ("bmi_rint", ), ("bmi_rint", ), ("bmi_rint", ), ("bmi_rint", ), ("bmi_rint", )]
    reruns = [False, False, False, False, True, True]

    for analysis, meta_class, phenos, rerun in zip(analyses, meta_classes, phenotypes, reruns):
        if rerun:
            meta_class = meta_class()
            for biobank in ["aou", "ukb"]:
                for ancestry in ["afr", "amr", "eas", "eur", "sas", "mid"]:
                    for pheno in phenos:
                        filename = f"{pheno}_{ancestry}_{biobank}_meta_w_samples.tsv.gz"
                        filepath = os.path.join(PROJECT_DIR, "data/meta/raw", analysis, filename)
                        meta_df = pd.read_csv(filepath, sep="\t")
                        processed_meta_df = meta_class.harmonize_gene_symbols_and_filter(meta_df, hgnc_df)
                        save_filename = f"{pheno}_{ancestry}_{biobank}.tsv.gz"
                        save_filepath =  os.path.join(PROJECT_DIR, "data/meta/processed", analysis, save_filename)
                        os.makedirs(os.path.dirname(save_filepath), exist_ok=True)
                        processed_meta_df.to_csv(save_filepath, sep="\t", index=False)
