import pandas as pd
import re
import os
from functools import reduce
import numpy as np
from scipy.stats import norm


def filter_df(df, gene):
    df = df.loc[df.snp=="gene"]
    df = df.rename(columns={"SE":  "se", "OR": "or", "p": "p_value"})
    df["gene"] = gene
    keep_cols = ["gene", "phenotype", "beta", "se", "or", "p_value", "n_total", "n_cases", "n_controls"]
    return df.loc[df.n_cases>50, keep_cols]

def create_meta_df(proj_dir, gene, ancestry=["afr", "amr", "eas", "eur", "sas", "mid"], biobank=["aou", "ukb"]):
    assoc_df = []
    merge_columns = ["gene", "phenotype"]
    stats_columns = ["beta", "se", "or", "p_value", "n_total", "n_cases", "n_controls"]
    for a in ancestry:
        for b in biobank:
            df = pd.read_csv(
                os.path.join(proj_dir, gene, f"phewas_{a}_{b}.csv.gz"), 
                dtype={"phenotype": str}
                )
            df = filter_df(df, gene)
            df.columns =[f"{c}_{a}_{b}" if c not in merge_columns else c for c in df.columns]
            assoc_df.append(df)

    meta_df = reduce(lambda x,y: x.merge(
        y, 
        on=merge_columns,
        how="outer"
        ), assoc_df
    )
    return meta_df

def get_meta_stats_helper(ser, ancestry, alpha=0.05):
    effect_sizes = np.array([ser[f"beta_{a}"] for a in ancestry if not pd.isnull(ser[f"beta_{a}"])])
    inverse_variances = np.array([(1/ser[f"se_{a}"]**2) for a in ancestry if not pd.isnull(ser[f"se_{a}"])])
    assert len(effect_sizes)==len(inverse_variances)

    if len(effect_sizes)==0:
        return pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA

    weighted_effect_size = np.sum(effect_sizes*inverse_variances)/np.sum(inverse_variances)
    weighted_variance = 1/np.sum(inverse_variances)
    weighted_se = np.sqrt(weighted_variance)

    # Calculate the Z-score
    z_score = weighted_effect_size / weighted_se
    # Calculate the two-tailed p-value
    p_value = 2 * norm.sf(abs(z_score)) #(1 - norm.cdf(abs(z_score)))

    # Calculate critical value for confidence interval
    z_critical = norm.ppf(1 - alpha / 2) 

    # Calculate confidence interval bounds
    lower_ci = weighted_effect_size - z_critical * weighted_se
    upper_ci = weighted_effect_size + z_critical * weighted_se

    n_cases = sum([ser[f"n_cases_{a}"] for a in ancestry if not  pd.isnull(ser[f"n_cases_{a}"])])
    n_controls = sum([ser[f"n_controls_{a}"] for a in ancestry if not  pd.isnull(ser[f"n_controls_{a}"])])
    n_total = sum([ser[f"n_total_{a}"] for a in ancestry if not  pd.isnull(ser[f"n_total_{a}"])])
    return weighted_effect_size, weighted_se, lower_ci, upper_ci, z_score, p_value, n_cases, n_controls, n_total

def get_meta_stats(ser):
    eur_ancestry = ["eur_aou", "eur_ukb"]
    noneur_ancestry =  ["afr_aou", "sas_aou", "eas_aou", "amr_aou",  "mid_aou", "afr_ukb", "sas_ukb", "eas_ukb", "amr_ukb",  "mid_ukb"]
    ees, ese, elci, ehci, ez_score, ep_value, ecases, econtrols, etotal = get_meta_stats_helper(ser, eur_ancestry)
    nees, nese, nelci, nehci, nez_score, nep_value, necases, necontrols, netotal = get_meta_stats_helper(ser, noneur_ancestry)
    es, se, lci, hci, z_score, p_value, cases, controls, total = get_meta_stats_helper(ser, eur_ancestry+noneur_ancestry)
    return pd.Series(
        {"gene": ser.gene, "phenotype": ser.phenotype, 
        "beta": es, "se": se, "ci_low": lci, "ci_high": hci, "z_score":z_score, "p_value": p_value, "cases": cases, "controls":  controls, "samples": total, 
        "ebeta": ees, "ese": ese, "eci_low": elci, "eci_high": ehci, "ez_score": ez_score, "ep_value": ep_value, "ecases": ecases, "econtrols":  econtrols, "esamples": etotal, 
        "nebeta": nees, "nese": nese, "neci_low": nelci, "neci_high": nehci, "nez_score": nez_score, "nep_value": nep_value, "necases": necases, "necontrols":  necontrols, "nesamples": netotal,
        })

def create_long_format_table(sig_meta_res_df):
    long_plot_df = pd.DataFrame()
    # Add meta statistics to long format
    plot_df = sig_meta_res_df.copy()
    plot_df = plot_df.set_index(["gene", "phenotype", "phenotype_def", "phenotype_cat"])
    plot_columns = [c for c in plot_df.columns if c!="gene"]
    meta_analysis = ["", "e", "ne"]
    meta_category = ["meta", "eur_meta", "non_eur_meta"]
    anc_category = ["All-ancestry", "European", "Non-european"]
    stats_cols = ["beta", "se", "ci_low", "ci_high", "p_value", "cases", "controls"]
    for m,c,ac in zip(meta_analysis, meta_category, anc_category):
        pdf = plot_df.loc[:, [f"{m}{s}" for s in stats_cols]]
        pdf.columns = stats_cols
        pdf["category"] = c
        pdf["anc_category"] = ac
        long_plot_df = pd.concat((long_plot_df, pdf))
    long_plot_df = long_plot_df.reset_index()
    long_plot_df = long_plot_df.rename(columns={"anc_category": "Group", "category": "Category"})
    return long_plot_df

def create_supp_table(long_plot_df):
    gene_columns = ["gene", "phenotype", "phenotype_def", "phenotype_cat"]
    columns_to_keep=["beta", "se", "ci_low", "ci_high", "p_value", "cases", "controls"]
    # Pivot the table
    df_pivot = long_plot_df.pivot(index=gene_columns, columns=["Group"], values=columns_to_keep)
    # Flatten the MultiIndex columns and group by gene and gene mask
    df_pivot.columns = pd.MultiIndex.from_tuples(
        [(g, stat) for stat, g in df_pivot.columns],
        names=["Group", "Statistic"]
    )
    df_pivot = df_pivot.loc[
        :, 
        [("European", s) for s in columns_to_keep] +
        [("Non-european", s) for s in columns_to_keep] +
        [("All-ancestry", s) for s in columns_to_keep]
    ]
    return df_pivot
    
if __name__=="__main__":
    PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"
    genes = ["YLPM1", "RIF1", "GIGYF1", "SLC5A3", "GRM7"]
    for gene in genes:
        gene_dir = os.path.join(PROJECT_DIR, f"data/phewas/")
        meta_df = create_meta_df(gene_dir, gene)
        BONF_P = 0.05/len(meta_df)
        phecode_def_df = pd.read_csv(os.path.join(PROJECT_DIR, "data/phewas/phecode_definitions1.2.csv"), dtype={"phecode": str})
        phecode_def_dict = phecode_def_df.loc[:, ["phecode", "phenotype"]].set_index("phecode").to_dict()["phenotype"]
        phecode_cat_dict = phecode_def_df.loc[:, ["phecode", "category"]].set_index("phecode").to_dict()["category"]
        meta_res_df = meta_df.apply(get_meta_stats, axis=1)
        meta_res_df["phenotype_def"] = meta_res_df.phenotype.map(phecode_def_dict)
        meta_res_df["phenotype_cat"] = meta_res_df.phenotype.map(phecode_cat_dict)
        phewas_file = os.path.join(PROJECT_DIR, f"data/phewas/{gene}/phewas_meta.csv.gz")
        meta_res_df.to_csv(phewas_file, index=False)

        sig_meta_res_df = meta_res_df.loc[
            (meta_res_df.p_value<BONF_P)
            ]
        long_plot_df = create_long_format_table(sig_meta_res_df)
        supp_table = create_supp_table(long_plot_df)
        supp_file = os.path.join(PROJECT_DIR, f"data/phewas/{gene}/phewas_meta.xlsx")
        supp_table.to_excel(supp_file)
