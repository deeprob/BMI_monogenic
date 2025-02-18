import os
import re
from functools import reduce
from scipy.stats import norm
import pandas as pd
import numpy as np


class Meta:
    def get_alternate_df(self, df, alternate_column):
        alternate_df = df.loc[:, ["hgnc_id", "symbol", alternate_column]]
        alternate_df = alternate_df.dropna(subset=[alternate_column])
        alternate_df[alternate_column] = alternate_df[alternate_column].str.split("|")
        alternate_df = alternate_df.explode(alternate_column)
        return alternate_df

    def get_hgnc_harmony_dict(self, hgnc_df):
        current_symbols = dict(zip(hgnc_df.symbol, hgnc_df.symbol))
        prev_df = self.get_alternate_df(hgnc_df, "prev_symbol")
        prev_symbols = dict(zip(prev_df.prev_symbol, prev_df.symbol))
        current_symbols.update(prev_symbols)
        return current_symbols

    def get_gene_mask(self, regenie_id):
        pattern = re.compile("(.+)\.(pLoF|Missense_strict|Missense_lenient)\.0\.001")
        m = re.match(pattern, regenie_id)
        if not m:
            print(regenie_id)
        gene = m.group(1)
        mask = m.group(2)
        return pd.Series({"gene": gene, "gene_mask": mask})

    def harmonize_gene_symbols_and_filter(self, meta_df, hgnc_df):
        # Add gene name and mask
        if len(meta_df)==0:
            meta_df[["gene", "gene_mask"]] =  pd.Series({"gene": "", "gene_mask": ""})
        else:
            meta_df[["gene", "gene_mask"]] = meta_df.ID.apply(self.get_gene_mask)
        # Add hgnc current version annotation
        hgnc_dict = self.get_hgnc_harmony_dict(hgnc_df)
        meta_df["hgnc_gene"] = meta_df.gene.map(hgnc_dict)
        # Keep genes which have 
        # chrom and gene pos info, 
        # hgnc annotations, 
        # no duplicate hgnc annotations and gene mask
        # at least one sample should be present for the gene mask pair: 
        # There might be some variant-sample pair which are not filtered from the bim files 
        # because the variant genotype were filtered as missing due to low quality in hail 
        # but not in bim - something to check
        meta_df = meta_df.loc[
            (meta_df.CHROM.notna())& 
            (meta_df.GENPOS.notna())&
            (meta_df.hgnc_gene.notna())&
            (~meta_df.duplicated(subset=["hgnc_gene", "gene_mask"], keep=False))&
            (meta_df.nsamples>0)
        ].drop(columns=["TEST", "EXTRA"])
        # replace REGENIE IDs with hgnc gene and gene mask
        meta_df["ID"] = meta_df.hgnc_gene + "::" + meta_df.gene_mask 
        return meta_df

    def create_meta_df(self, proj_dir, ancestry=["afr", "amr", "eas", "eur", "sas", "mid"], biobank=["aou", "ukb"]):
        assoc_df = []
        merge_columns = ["hgnc_gene", "gene_mask", "ID", "ALLELE0", "ALLELE1"]
        stats_columns = ["N", "BETA", "SE", "CHISQ", "LOG10P", "p_value", "nsamples"]
        for a in ancestry:
            for b in biobank:
                df = pd.read_csv(
                    os.path.join(proj_dir, f"bmi_rint_{a}_{b}.tsv.gz"), 
                    sep="\t", usecols=merge_columns + stats_columns
                    )
                df.columns =[f"{c}_{a}_{b}" if c not in merge_columns else c for c in df.columns]
                assoc_df.append(df)

        meta_df = reduce(lambda x,y: x.merge(
            y, 
            on=merge_columns,
            how="outer"
            ), assoc_df
        )
        return meta_df

    def get_meta_stats_helper(self, ser, ancestry, alpha=0.05):
        effect_sizes = np.array([ser[f"BETA_{a}"] for a in ancestry if not pd.isnull(ser[f"BETA_{a}"])])
        inverse_variances = np.array([(1/ser[f"SE_{a}"]**2) for a in ancestry if not pd.isnull(ser[f"SE_{a}"])])
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

        nsamples = sum([ser[f"nsamples_{a}"] for a in ancestry if not  pd.isnull(ser[f"nsamples_{a}"])])
        return weighted_effect_size, weighted_se, lower_ci, upper_ci, z_score, p_value, nsamples

    def get_meta_stats(self, ser):
        eur_ancestry = ["eur_aou", "eur_ukb"]
        noneur_ancestry =  ["afr_aou", "sas_aou", "eas_aou", "amr_aou",  "mid_aou", "afr_ukb", "sas_ukb", "eas_ukb", "amr_ukb",  "mid_ukb"]
        ees, ese, elci, ehci, ez_score, ep_value, esamples =self.get_meta_stats_helper(ser, eur_ancestry)
        nees, nese, nelci, nehci, nez_score, nep_value, nesamples = self.get_meta_stats_helper(ser, noneur_ancestry)
        es, se, lci, hci, z_score, p_value, nsamples = self.get_meta_stats_helper(ser, eur_ancestry+noneur_ancestry)
        return pd.Series({
            "gene": ser.hgnc_gene, "gene_mask": ser.gene_mask, 
            "beta": es, "se": se, "ci_low": lci, "ci_high": hci, "z_score":z_score, "p_value": p_value, "nsamples": nsamples,
            "ebeta": ees, "ese": ese, "eci_low": elci, "eci_high": ehci, "ez_score": ez_score, "ep_value": ep_value, "esamples": esamples,
            "nebeta": nees, "nese": nese, "neci_low": nelci, "neci_high": nehci, "nez_score": nez_score, "nep_value": nep_value, "nesamples": nesamples,
        })

    def perform_meta_analysis(self, proj_dir, ancestry=["afr", "amr", "eas", "eur", "sas", "mid"], biobank=["aou", "ukb"]):
        meta_df = self.create_meta_df(proj_dir)
        meta_res_df = meta_df.apply(self.get_meta_stats, axis=1)
        return meta_res_df


class Lovo(Meta):
    def get_gene_mask_variant(self, regenie_id):
        pattern = re.compile("(.+)\.(pLoF|Missense_strict|Missense_lenient)\.0\.001_(.+)")
        m = re.match(pattern, regenie_id)
        if not m:
            pattern = re.compile("(.+)\.(pLoF|Missense_strict|Missense_lenient)\.0\.001")
            m = re.match(pattern, regenie_id)
            variant = "all"
        else:
            variant = m.group(3).lstrip("chr")
        gene = m.group(1)
        mask = m.group(2)
        return pd.Series({"gene": gene, "gene_mask": mask, "lovo_variant": variant})
    
    def harmonize_gene_symbols_and_filter(self, meta_df, hgnc_df):
        # Add gene name and mask
        if len(meta_df)==0:
            meta_df[["gene", "gene_mask", "lovo_variant"]] =  pd.Series({"gene": "", "gene_mask": "", "lovo_variant": ""})
        else:
            meta_df[["gene", "gene_mask", "lovo_variant"]] = meta_df.ID.apply(self.get_gene_mask_variant)
        # Add hgnc current version annotation
        hgnc_dict = self.get_hgnc_harmony_dict(hgnc_df)
        meta_df["hgnc_gene"] = meta_df.gene.map(hgnc_dict)
        # Keep genes which have 
        # chrom and gene pos info, 
        # hgnc annotations, 
        # no duplicate hgnc annotations and gene mask
        # at least one sample should be present for the gene mask pair: 
        # There might be some variant-sample pair which are not filtered from the bim files 
        # because the variant genotype were filtered as missing due to low quality in hail 
        # but not in bim - something to check
        meta_df = meta_df.loc[
            (meta_df.CHROM.notna())& 
            (meta_df.GENPOS.notna())&
            (meta_df.hgnc_gene.notna())&
            (~meta_df.duplicated(subset=["hgnc_gene", "gene_mask", "lovo_variant"], keep=False))
        ].drop(columns=["TEST", "EXTRA"])
        # replace REGENIE IDs with hgnc gene and gene mask
        meta_df["ID"] = meta_df.hgnc_gene + "::" + meta_df.gene_mask + "::" + meta_df.lovo_variant
        return meta_df
    
    def get_meta_stats_helper(self, ser, ancestry, alpha=0.05):
        effect_sizes = np.array([ser[f"BETA_{a}"] for a in ancestry if not pd.isnull(ser[f"BETA_{a}"])])
        inverse_variances = np.array([(1/ser[f"SE_{a}"]**2) for a in ancestry if not pd.isnull(ser[f"SE_{a}"])])
        assert len(effect_sizes)==len(inverse_variances)

        if len(effect_sizes)==0:
            return pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA

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

        return weighted_effect_size, weighted_se, lower_ci, upper_ci, z_score, p_value

    def get_meta_stats(self, ser):
        eur_ancestry = ["eur_aou", "eur_ukb"]
        noneur_ancestry =  ["afr_aou", "sas_aou", "eas_aou", "amr_aou",  "mid_aou", "afr_ukb", "sas_ukb", "eas_ukb", "amr_ukb",  "mid_ukb"]
        ees, ese, elci, ehci, ez_score, ep_value = self.get_meta_stats_helper(ser, eur_ancestry)
        nees, nese, nelci, nehci, nez_score, nep_value = self.get_meta_stats_helper(ser, noneur_ancestry)
        es, se, lci, hci, z_score, p_value = self.get_meta_stats_helper(ser, eur_ancestry+noneur_ancestry)
        return pd.Series(
            {"gene": ser.hgnc_gene, "gene_mask": ser.gene_mask, "lovo_variant": ser.lovo_variant,
            "beta": es, "se": se, "ci_low": lci, "ci_high": hci, "z_score":z_score, "p_value": p_value,
            "ebeta": ees, "ese": ese, "eci_low": elci, "eci_high": ehci, "ez_score": ez_score, "ep_value": ep_value,
            "nebeta": nees, "nese": nese, "neci_low": nelci, "neci_high": nehci, "nez_score": nez_score, "nep_value": nep_value
            })

    def create_meta_df(self, proj_dir, ancestry=["afr", "amr", "eas", "eur", "sas", "mid"], biobank=["aou", "ukb"]):
        assoc_df = []
        merge_columns = ["hgnc_gene", "gene_mask", "ID", "ALLELE0", "ALLELE1", "lovo_variant"]
        stats_columns = ["N", "BETA", "SE", "CHISQ", "LOG10P", "p_value"]
        for a in ancestry:
            for b in biobank:
                df = pd.read_csv(
                    os.path.join(proj_dir, f"bmi_rint_{a}_{b}.tsv.gz"), 
                    sep="\t", usecols=merge_columns + stats_columns
                    )
                df.columns =[f"{c}_{a}_{b}" if c not in merge_columns else c for c in df.columns]
                assoc_df.append(df)

        meta_df = reduce(lambda x,y: x.merge(
            y, 
            on=merge_columns,
            how="outer"
            ), assoc_df
        )
        return meta_df

    def perform_meta_analysis(self, proj_dir, ancestry=["afr", "amr", "eas", "eur", "sas", "mid"], biobank=["aou", "ukb"]):
        meta_df = self.create_meta_df(proj_dir)
        meta_df_filled = meta_df.groupby(["hgnc_gene", "gene_mask"], group_keys=False).apply(
            lambda g: g.fillna(g.loc[g.lovo_variant == "all", :].squeeze())
            )
        meta_res_df = meta_df_filled.apply(self.get_meta_stats, axis=1)
        return meta_res_df


class Meds(Meta):
    def harmonize_gene_symbols_and_filter(self, meta_df, hgnc_df):
        # Add gene name and mask
        if len(meta_df)==0:
            meta_df[["gene", "gene_mask"]] =  pd.Series({"gene": "", "gene_mask": ""})
        else:
            meta_df[["gene", "gene_mask"]] = meta_df.ID.apply(self.get_gene_mask)
        # Add hgnc current version annotation
        hgnc_dict = self.get_hgnc_harmony_dict(hgnc_df)
        meta_df["hgnc_gene"] = meta_df.gene.map(hgnc_dict)
        # Keep genes which have 
        # chrom and gene pos info, 
        # hgnc annotations, 
        # no duplicate hgnc annotations and gene mask
        # at least one sample should be present for the gene mask pair: 
        # There might be some variant-sample pair which are not filtered from the bim files 
        # because the variant genotype were filtered as missing due to low quality in hail 
        # but not in bim - something to check
        meta_df = meta_df.loc[
            (meta_df.CHROM.notna())& 
            (meta_df.GENPOS.notna())&
            (meta_df.hgnc_gene.notna())&
            (~meta_df.duplicated(subset=["hgnc_gene", "gene_mask"], keep=False))
        ].drop(columns=["TEST", "EXTRA"])
        # replace REGENIE IDs with hgnc gene and gene mask
        meta_df["ID"] = meta_df.hgnc_gene + "::" + meta_df.gene_mask
        return meta_df

    def create_meta_df(self, proj_dir, ancestry=["afr", "amr", "eas", "eur", "sas", "mid"], biobank=["aou", "ukb"]):
        assoc_df = []
        merge_columns = ["hgnc_gene", "gene_mask", "ID", "ALLELE0", "ALLELE1"]
        stats_columns = ["N", "BETA", "SE", "CHISQ", "LOG10P", "p_value"]
        for a in ancestry:
            for b in biobank:
                df = pd.read_csv(
                    os.path.join(proj_dir, f"bmi_rint_{a}_{b}.tsv.gz"), 
                    sep="\t", usecols=merge_columns + stats_columns
                    )
                df.columns =[f"{c}_{a}_{b}" if c not in merge_columns else c for c in df.columns]
                assoc_df.append(df)

        meta_df = reduce(lambda x,y: x.merge(
            y, 
            on=merge_columns,
            how="outer"
            ), assoc_df
        )
        return meta_df

    def get_meta_stats_helper(self, ser, ancestry, alpha=0.05):
        effect_sizes = np.array([ser[f"BETA_{a}"] for a in ancestry if not pd.isnull(ser[f"BETA_{a}"])])
        inverse_variances = np.array([(1/ser[f"SE_{a}"]**2) for a in ancestry if not pd.isnull(ser[f"SE_{a}"])])
        assert len(effect_sizes)==len(inverse_variances)
        if len(effect_sizes)==0:
            return pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA
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
        return weighted_effect_size, weighted_se, lower_ci, upper_ci, z_score, p_value

    def get_meta_stats(self, ser):
        eur_ancestry = ["eur_aou", "eur_ukb"]
        noneur_ancestry =  ["afr_aou", "sas_aou", "eas_aou", "amr_aou",  "mid_aou", "afr_ukb", "sas_ukb", "eas_ukb", "amr_ukb",  "mid_ukb"]
        ees, ese, elci, ehci, ez_score, ep_value = self.get_meta_stats_helper(ser, eur_ancestry)
        nees, nese, nelci, nehci, nez_score, nep_value = self.get_meta_stats_helper(ser, noneur_ancestry)
        es, se, lci, hci, z_score, p_value = self.get_meta_stats_helper(ser, eur_ancestry+noneur_ancestry)
        return pd.Series(
            {"gene": ser.hgnc_gene, "gene_mask": ser.gene_mask,
            "beta": es, "se": se, "ci_low": lci, "ci_high": hci, "z_score":z_score, "p_value": p_value,
            "ebeta": ees, "ese": ese, "eci_low": elci, "eci_high": ehci, "ez_score": ez_score, "ep_value": ep_value,
            "nebeta": nees, "nese": nese, "neci_low": nelci, "neci_high": nehci, "nez_score": nez_score, "nep_value": nep_value
            })

