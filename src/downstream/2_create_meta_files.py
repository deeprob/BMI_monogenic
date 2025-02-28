import os
import pandas as pd
from utils import meta

PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"

def get_most_deleterious_idx(ser):
    all_masks = set(ser.unique())
    most_del = "Missense_lenient"
    if "pLoF" in all_masks:
        most_del = "pLoF"
    elif "Missense_strict" in all_masks:
        most_del = "Missense_strict"
    most_del_idx = ser.loc[ser==most_del].index[0]
    return most_del_idx

if __name__ == '__main__':
    analyses = ["all_ancestry", "meds", "shadow_effect", "lovo"]
    meta_classes = [meta.Meta, meta.Meds, meta.Meds, meta.Lovo]
    include_samples = [ True, False, False, False]

    for analysis, meta_class, isamples in zip(analyses, meta_classes, include_samples):
        print(analysis)
        filename = os.path.join(PROJECT_DIR, f"data/meta/results/{analysis}/ivw_fixed/bmi_rint_monogenic_meta.tsv")
        sig_meta_res_df = pd.read_csv(filename, sep="\t")
        meta_class = meta_class()
        if analysis == "lovo":
            meta_df = meta_class.create_supp_table(sig_meta_res_df)
        else:
            most_del_sig_meta_df = sig_meta_res_df.loc[sig_meta_res_df.groupby("gene")["gene_mask"].apply(get_most_deleterious_idx)].sort_values("p_value").reset_index(drop=True)
            anc_dir = os.path.join(PROJECT_DIR, f"data/meta/processed/{analysis}/")
            meta_df = meta_class.create_supp_table(most_del_sig_meta_df, anc_dir, isamples)
        meta_df_filename = os.path.join(PROJECT_DIR, f"data/meta/tables/{analysis}.xlsx")
        os.makedirs(os.path.dirname(meta_df_filename), exist_ok=True)
        meta_df.to_excel(meta_df_filename)
