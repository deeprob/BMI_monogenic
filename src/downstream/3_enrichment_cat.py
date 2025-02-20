import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from scipy import stats


def get_table_icd_helper(ser, comorbidity):
    table = [
        [ser[f"gene_{comorbidity}"], ser[f"gene_non{comorbidity}"]],
        [ser[f"nongene_{comorbidity}"], ser[f"nongene_non{comorbidity}"]]
    ]
    res = fisher_exact(table)
    or_study = odds_ratio(table)
    cil, cih = or_study.confidence_interval(confidence_level=0.95)
    return pd.Series({"gene": ser.gene, "gene_mask": ser.gene_mask, "comorbidity": comorbidity, "OR": res.statistic, "p_value": res.pvalue, "ci_low": cil, "ci_high": cih})

def get_most_deleterious_idx(ser):
    all_masks = set(ser.unique())
    most_del = "Missense_lenient"
    if "pLoF" in all_masks:
        most_del = "pLoF"
    elif "Missense_strict" in all_masks:
        most_del = "Missense_strict"
    most_del_idx = ser.loc[ser==most_del].index[0]
    return most_del_idx

def conduct_enrichment(aou_file, ukb_file, comorbidities, most_del_sig_meta_df):
    aou_df = pd.read_csv(aou_file, index_col=[0, 1])
    ukb_df = pd.read_csv(ukb_file, index_col=[0, 1])
    meta_df = aou_df+ukb_df
    meta_df = meta_df.reset_index()
    meta_df = meta_df.rename(columns={"mask": "gene_mask"})
    meta_df = meta_df.merge(most_del_sig_meta_df.loc[:, ["gene", "gene_mask"]], on=["gene", "gene_mask"])
    df = pd.DataFrame()
    for comorbidity in comorbidities:
        cdf = meta_df.apply(get_table_icd_helper, axis=1, args=(comorbidity,))
        df = pd.concat((df, cdf))
    df = df.reset_index(drop=True)
    return df

def create_supp_table(enrich_df, most_del_sig_meta_df):
    stats_cols = ["OR", "ci_low", "ci_high", "p_value"]
    df_pivot = enrich_df.pivot(index=["gene", "gene_mask"], columns="comorbidity", values=stats_cols)
        # Flatten the MultiIndex columns and group by gene and gene mask
    df_pivot.columns = pd.MultiIndex.from_tuples(
        [(com, stat) for stat, com in df_pivot.columns],
        names=["Comorbidity", "Statistic"]
    )
    df_pivot = df_pivot.loc[
        [(g,m) for g,m in most_del_sig_meta_df.sort_values("beta", ascending=False).loc[:, ["gene", "gene_mask"]].values],
        [(com, stat) for com in comorbidities for stat in stats_cols ] 
    ]
    return df_pivot


if __name__=="__main__":
    PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"
    sig_meta_file = os.path.join(PROJECT_DIR, f"data/meta/results/all_ancestry/ivw_fixed/bmi_rint_monogenic_meta.tsv")
    sig_meta_res_df = pd.read_csv(sig_meta_file, sep="\t")
    most_del_sig_meta_df = sig_meta_res_df.loc[sig_meta_res_df.groupby("gene")["gene_mask"].apply(get_most_deleterious_idx)].sort_values("p_value").reset_index(drop=True)

    comorbid_list = [
        ["nu", "ovw", "ob", "sob"],
        ["cvd", "cad", "ht", "t2d", "hf", "af", "pe", "vt", "avs", "grd", "cls", "ccs", "cd", "nfld", "koa"]
    ]
    for analysis, comorbidities in zip(["bmi_cat", "comorbid"], comorbid_list):
        aou_file = os.path.join(PROJECT_DIR, f"data/enrichment/{analysis}/monogenic_aou_{analysis}.csv.gz")
        ukb_file = os.path.join(PROJECT_DIR, f"data/enrichment/{analysis}/monogenic_ukb_{analysis}.csv.gz")
        enrich_df = conduct_enrichment(aou_file, ukb_file, comorbidities, most_del_sig_meta_df)
        sig_enrich_file = os.path.join(PROJECT_DIR, f"data/enrichment/{analysis}/monogenic_enrich_{analysis}.csv")
        BONF_P = 0.05/len(enrich_df)
        enrich_df.loc[enrich_df.p_value<BONF_P].to_csv(sig_enrich_file, index=False)
        df_pivot = create_supp_table(enrich_df, most_del_sig_meta_df)
        supp_file = os.path.join(PROJECT_DIR, f"data/enrichment/{analysis}/monogenic_enrich_{analysis}.xlsx")
        df_pivot.to_excel(supp_file)
