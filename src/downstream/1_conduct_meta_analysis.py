import os
from utils import meta

BONF_P = 0.05/(20000*3)
PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"



if __name__=="__main__":
    analyses = ["all_ancestry", "lovo", "meds", "shadow_effect"]
    meta_classes = [meta.Meta, meta.Lovo, meta.Meds, meta.Meds]

    for analysis, meta_class  in zip(analyses, meta_classes):
        print(analysis)
        # create meta df
        proj_dir = os.path.join(PROJECT_DIR, "data/meta/processed", analysis)
        meta_class = meta_class()
        meta_res_df = meta_class.perform_meta_analysis(proj_dir)
        meta_res_filename = "meta_results.tsv.gz"
        meta_res_filepath = os.path.join(PROJECT_DIR, "data/meta/results", analysis, "ivw_fixed", meta_res_filename)
        os.makedirs(os.path.dirname(meta_res_filepath), exist_ok=True)
        meta_res_df.to_csv(meta_res_filepath, index=False, sep="\t")
        sig_meta_res_df = meta_res_df.loc[
            ((meta_res_df.ep_value<BONF_P)|(meta_res_df.nep_value<BONF_P))&
            (meta_res_df.p_value<BONF_P)
            ]
        sig_meta_res_filename = "monogenic_meta.tsv"
        sig_meta_res_filepath = os.path.join(PROJECT_DIR, "data/meta/results", analysis, "ivw_fixed", sig_meta_res_filename)
        sig_meta_res_df.to_csv(sig_meta_res_filepath, index=False, sep="\t")

