import pandas as pd

if __name__ == "__main__":
    meta_file = "../data/meta/tables/all_ancestry_with_nvariants.xlsx"
    meta_df = pd.read_excel(meta_file, index_col=[0,1], header=[0,1,2])

    aou_nvars = pd.read_csv("../data/meta/tables/per_gene_variants_aou.csv", index_col=0)
    ukb_nvars = pd.read_csv("../data/meta/tables/per_gene_variants_ukb.csv", index_col=0).drop(columns=["oth"])

    aou_nvars.columns = [f"{col}_aou" for col in aou_nvars.columns]
    ukb_nvars.columns = [f"{col}_ukb" for col in ukb_nvars.columns]

    nvars_df = aou_nvars.merge(ukb_nvars, left_index=True, right_index=True, how="outer").T
    variant_type_map = dict(zip(meta_df.columns.get_level_values(0), meta_df.columns.get_level_values(1)))
    category_map = dict(zip(meta_df.index.get_level_values(1), meta_df.index.get_level_values(0)))

    nvars_df.columns = pd.MultiIndex.from_tuples([(g, variant_type_map[g], "nvariants") for g in nvars_df.columns], names=["Gene", "Gene Mask", "Statistic"])
    nvars_df.index = pd.MultiIndex.from_tuples([(category_map[c], c) for c in nvars_df.index], names = ["Group", "Category"])

    group_sums = nvars_df.groupby(level=0).sum()
    group_sums.index = pd.MultiIndex.from_tuples(
        [("European",    "eur_meta"),
        ("Non-european","non_eur_meta")],
        names=nvars_df.index.names
    )
    all_meta = group_sums.sum(axis=0).to_frame().T
    all_meta.index = pd.MultiIndex.from_tuples(
        [("All-ancestry","meta")],
        names=nvars_df.index.names
    )
    nvars_with_meta = pd.concat([nvars_df, group_sums, all_meta])

    meta_with_nvariants = meta_df.merge(nvars_with_meta, how="left", left_index=True, right_index=True)
    stat_order = ["beta","se","ci_low","ci_high","p_value","nsamples","nvariants"]
    genes = meta_df.columns.get_level_values(0).unique()
    masks_dict = dict(zip(meta_df.columns.get_level_values(0), meta_df.columns.get_level_values(1)))

    new_cols = []
    for gene in genes:
        for stat in stat_order:
            new_cols.append((gene, masks_dict[gene], stat))

    new_index = pd.MultiIndex.from_tuples(
        new_cols,
        names=meta_with_nvariants.columns.names  # should be ["Gene","Gene Mask","Statistic"]
    )
    meta_with_nvariants = meta_with_nvariants.reindex(columns=new_index)

    meta_with_nvariants.to_excel("../data/meta/tables/all_ancestry_with_nvariants.xlsx")

