import os
from utils import meta

BONF_P = 0.05/(20000*3)
PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"



if __name__=="__main__":
    analyses = ["all_ancestry", "lovo", "meds", "shadow_effect", "male", "female"]
    meta_classes = [meta.Meta, meta.Lovo, meta.Meds, meta.Meds, meta.Meta, meta.Meta]
    phenotypes =[("bmi_rint", "bmi"), ("bmi_rint", ), ("bmi_rint", ), ("bmi_rint", ), ("bmi_rint", ), ("bmi_rint", )]
    reruns = [False, False, False, False, True, True]

    for analysis, meta_class, phenos, rerun in zip(analyses, meta_classes, phenotypes, reruns):
        if rerun:
            print(analysis)
            # create meta df
            proj_dir = os.path.join(PROJECT_DIR, "data/meta/processed", analysis)
            meta_class = meta_class()
            for pheno in phenos:
                meta_res_df = meta_class.perform_meta_analysis(proj_dir, pheno)
                meta_res_filename = f"{pheno}_meta_results.tsv.gz"
                meta_res_filepath = os.path.join(PROJECT_DIR, "data/meta/results", analysis, "ivw_fixed", meta_res_filename)
                os.makedirs(os.path.dirname(meta_res_filepath), exist_ok=True)
                meta_res_df.to_csv(meta_res_filepath, index=False, sep="\t")
                sig_meta_res_df = meta_res_df.loc[
                    ((meta_res_df.ep_value<BONF_P)|(meta_res_df.nep_value<BONF_P))&
                    (meta_res_df.p_value<BONF_P)
                    ]
                sig_meta_res_filename = f"{pheno}_monogenic_meta.tsv"
                sig_meta_res_filepath = os.path.join(PROJECT_DIR, "data/meta/results", analysis, "ivw_fixed", sig_meta_res_filename)
                sig_meta_res_df.to_csv(sig_meta_res_filepath, index=False, sep="\t")

