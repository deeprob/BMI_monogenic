import pandas as pd
import re
import os
from functools import reduce
import numpy as np
from scipy.stats import norm

def get_meta_stats_helper(ser, ancestry, alpha=0.05):
    effect_sizes = np.array([ser[f"est_{a}"] for a in ancestry if not pd.isnull(ser[f"est_{a}"])])
    inverse_variances = np.array([(1/ser[f"se_{a}"]**2) for a in ancestry if not pd.isnull(ser[f"se_{a}"])])
    assert len(effect_sizes)==len(inverse_variances)

    if len(effect_sizes)==0:
        return pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA

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

    nobs = sum([ser[f"nobs_{a}"] for a in ancestry if not  pd.isnull(ser[f"nobs_{a}"])])
    return weighted_effect_size, weighted_se, lower_ci, upper_ci, z_score, p_value, nobs

def get_meta_stats(ser):
    biobank = ["aou", "ukb"]
    es, se, lci, hci, z_score, p_value, nsamples = get_meta_stats_helper(ser, biobank)
    return pd.Series(
        {"gene": ser.gene, "comorbidity": ser.comorbidity, 
        "lhs": ser.lhs, "op": ser.op, "rhs": ser.rhs, "exo": ser.exo, "label": ser.label,
        "beta": es, "se": se, "ci_low": lci, "ci_high": hci, "z_score":z_score, "p_value": p_value, "nobs": nsamples
        })

def create_meta_df(aou_file, ukb_file):
    aou_df = pd.read_csv(aou_file)
    ukb_df = pd.read_csv(ukb_file)
    meta_df = aou_df.merge(
        ukb_df, 
        on=["gene", "comorbidity", "label", "lhs", "rhs", "op", "exo"], 
        suffixes=("_aou", "_ukb")
        )
    meta_df = meta_df.loc[meta_df.label.notna()].reset_index(drop=True)
    return meta_df

if __name__=="__main__":
    PROJECT_DIR = "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"
    aou_file = os.path.join(PROJECT_DIR, "data/sem/sem_aou.csv.gz")
    ukb_file = os.path.join(PROJECT_DIR, "data/sem/sem_ukb.csv.gz")
    meta_df = create_meta_df(aou_file, ukb_file)
    meta_res_df = meta_df.apply(get_meta_stats, axis=1)
    meta_res_df = meta_res_df.set_index(["gene", "comorbidity"])
    meta_res_df.to_excel(os.path.join(PROJECT_DIR, "data/sem/sem_meta.xlsx"))
