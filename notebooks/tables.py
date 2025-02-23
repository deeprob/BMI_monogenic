import os
import pandas as pd
from scipy.stats import norm, chi2


def calculate_ci(ser, alpha=0.05):
    # Calculate critical value for confidence interval
    z_critical = norm.ppf(1 - alpha / 2) 
    # Calculate confidence interval bounds
    lower_ci = ser.BETA - z_critical * ser.SE
    upper_ci = ser.BETA + z_critical * ser.SE
    return pd.Series({"ci_low": lower_ci, "ci_high": upper_ci})

def get_anc_stats(sig_df, anc_df, include_samples=False):
    sig_merge_cols = ["gene", "gene_mask"]
    anc_merge_cols = ["hgnc_gene", "gene_mask"]
    anc_stats_cols=["BETA", "SE", "p_value"]
    if include_samples: anc_stats_cols.append("nsamples")
    merged_df = sig_df.loc[:, sig_merge_cols].merge(anc_df.loc[:, anc_merge_cols + anc_stats_cols], left_on=sig_merge_cols, right_on=anc_merge_cols)
    if len(merged_df)>0:
        merged_df[["ci_low", "ci_high"]] = merged_df.apply(calculate_ci, axis=1)
    else:
        merged_df[["ci_low", "ci_high"]] = pd.NA, pd.NA
    return merged_df.rename(columns={"BETA": "beta", "SE": "se"}).drop(columns=["hgnc_gene"]).set_index(["gene", "gene_mask"])

def create_long_format_table(most_del_sig_meta_df, anc_dir, include_samples=False):
    long_plot_df = pd.DataFrame()
    # Add meta statistics to long format
    plot_df = most_del_sig_meta_df.copy()
    plot_df = plot_df.set_index(["gene", "gene_mask"])
    plot_df = plot_df.rename(columns={"esamples": "ensamples", "nesamples": "nensamples"})
    plot_columns = [c for c in plot_df.columns if c!="gene"]
    meta_analysis = ["", "e", "ne"]
    meta_category = ["meta", "eur_meta", "non_eur_meta"]
    anc_category = ["All-ancestry", "European", "Non-european"]
    stats_cols = ["beta", "se", "ci_low", "ci_high", "p_value"]
    if include_samples: stats_cols.append("nsamples")
    for m,c,ac in zip(meta_analysis, meta_category, anc_category):
        pdf = plot_df.loc[:, [f"{m}{s}" for s in stats_cols]]
        pdf.columns = stats_cols
        pdf["category"] = c
        pdf["anc_category"] = ac
        long_plot_df = pd.concat((long_plot_df, pdf))
    # Add ancestry wise statistics to long format
    ancestry = ["afr", "amr", "eas", "eur", "mid", "sas"]
    anc_category = ["Non-european", "Non-european", "Non-european", "European", "Non-european", "Non-european"]
    biobank = ["ukb", "aou"]
    # Add ancestry wise statistics to long format
    ancestry = ["afr", "amr", "eas", "eur", "mid", "sas"]
    anc_category = ["Non-european", "Non-european", "Non-european", "European", "Non-european", "Non-european"]
    biobank = ["ukb", "aou"]
    missing_cat = []
    for a,ac in zip(ancestry, anc_category):
        for b in biobank:
            anc_file = os.path.join(anc_dir, f"bmi_rint_{a}_{b}.tsv.gz")
            anc_df = pd.read_csv(anc_file, sep="\t")
            merged_anc_df = get_anc_stats(most_del_sig_meta_df, anc_df, include_samples)
            if len(merged_anc_df)==0:
                missing_cat.append(f"{a}_{b}")
                continue
            merged_anc_df["category"] = f"{a}_{b}"
            merged_anc_df["anc_category"] = ac
            long_plot_df = pd.concat((long_plot_df, merged_anc_df))
    long_plot_df = long_plot_df.reset_index()
    long_plot_df = long_plot_df.rename(columns={"anc_category": "Group", "category": "Category"})
    return long_plot_df, set(missing_cat)

def pivot_table(long_plot_df, include_samples=False):
    gene_columns = ["gene", "gene_mask"]
    columns_to_keep=["beta", "se", "ci_low", "ci_high", "p_value"]
    if include_samples: columns_to_keep.append("nsamples")
    # Pivot the table
    df_pivot = long_plot_df.pivot(columns=gene_columns, index=["Group", "Category"], values=columns_to_keep)
    # Flatten the MultiIndex columns and group by gene and gene mask
    df_pivot.columns = pd.MultiIndex.from_tuples(
        [(gene, mask, stat) for stat, gene, mask in df_pivot.columns],
        names=["Gene", "Gene Mask", "Statistic"]
    )
    df_pivot = df_pivot.sort_index(axis=1, level=["Gene", "Gene Mask"])
    return df_pivot


def create_supp_table(most_del_sig_meta_df, anc_dir, include_samples=False):
    columns_to_keep=["beta", "se", "ci_low", "ci_high", "p_value"]
    if include_samples: columns_to_keep.append("nsamples")
    long_plot_df, missing_cat = create_long_format_table(most_del_sig_meta_df, anc_dir, include_samples)
    df_pivot = pivot_table(long_plot_df, include_samples)
    df_pivot = df_pivot.loc[
        [("European", "eur_aou"), ("European", "eur_ukb"), ("European", "eur_meta")] + 
        [("Non-european", f"{a}_{b}") for a in ["afr", "amr", "eas", "mid", "sas"] for b in ["aou", "ukb"] if not f"{a}_{b}" in missing_cat] + [("Non-european", "non_eur_meta")] +
        [("All-ancestry", "meta")], [(g,m,s) for g,m in most_del_sig_meta_df.sort_values("beta", ascending=False).loc[:, ["gene", "gene_mask"]].values for s in columns_to_keep]
    ]
    return df_pivot
