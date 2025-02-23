import pandas as pd
import numpy as np
import itertools as it
import os

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial" # missing fonts:: https://alexanderlabwhoi.github.io/post/2021-03-missingfont/
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams.update({'font.size': 6, 'axes.linewidth': 1, 'xtick.major.width': 1, 'xtick.major.size': 5, 'ytick.major.width': 1, 'ytick.major.size': 5})
from matplotlib.backends.backend_pdf import PdfPages


def create_df_for_forestplot(most_del_known_df):
    plot_df = most_del_known_df.copy()
    plot_df = plot_df.set_index("gene")
    plot_columns = [c for c in plot_df.columns if c!="gene"]

    meta_analysis = ["", "e", "ne"]
    meta_category = ["all", "european", "non-european"]
    stats_cols = ["beta", "se", "ci_low", "ci_high", "z_score", "p_value"]

    long_plot_df = pd.DataFrame()

    for m,c in zip(meta_analysis, meta_category):
        pdf = plot_df.loc[:, [f"{m}{s}" for s in stats_cols]]
        pdf.columns = stats_cols
        pdf["category"] = c
        long_plot_df = pd.concat((long_plot_df, pdf))

    long_plot_df = long_plot_df.reset_index()
    long_plot_df["formatted_beta"] = long_plot_df.beta.apply(lambda x: f"{x:.2f}")
    long_plot_df["formatted_ci"] = "(" + long_plot_df.ci_low.apply(lambda x: f"{x:.2f}") + ", " + long_plot_df.ci_high.apply(lambda x: f"{x:.2f}") + ")"
    long_plot_df["formatted_se"] = long_plot_df.se.apply(lambda x: f"{x:.2f}")
    long_plot_df["formatted_z_score"] = long_plot_df.z_score.apply(lambda x: f"{x:.2f}")
    long_plot_df["formatted_p_value"] = long_plot_df.p_value.apply(p_value_formatter)
    long_plot_df["Beta (95% CI)"] = long_plot_df.apply(lambda ser: f"{ser.formatted_beta} {ser.formatted_ci}", axis=1)
    return long_plot_df

def p_value_formatter(pval):
    if pval==0:
        pval = "0"
    else:
        pval = f"{pval:.2e}"
        pval = pval.replace("e", " x 10")
    return pval

def save_pdf(save_file, fig):
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    pdf = PdfPages(save_file)
    pdf.savefig(fig, bbox_inches='tight',dpi=300)
    pdf.close()
    return

def create_forestplot(
    df, studies, categories, 
    labels_col, categories_col, effect_sizes_col, ci_low_col, ci_high_col,
    stats_cols, figsize=(3, 2.5)
):
    # Define markers for each category
    markers = dict(it.zip_longest(categories, "osD", fillvalue="o"))
    # Define colors for each category
    palette = ["royalblue", "indianred", "gray"] # ["lightgrey", "darkgrey", "black"]

    colors = dict(it.zip_longest(categories, palette, fillvalue="black"))

    # Create a figure with two axes
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [1.5, 2]},sharey=True)

    # Horizontal line at 0 for the null effect
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1)
    
    df["ci_low_error"] = df[effect_sizes_col] - df[ci_low_col]
    df["ci_high_error"] =  df[ci_high_col] - df[effect_sizes_col]
    study_offset=1.25
    cat_offset=1
    ### main plot ###
    # Plotting the effect sizes with different markers per category
    for i, study in enumerate(studies):
        offset = study_offset*i
        for j, category in enumerate(categories):
            # Get the index for this category
            idx = i * len(categories) + j*cat_offset
            effect_size = df.loc[(df[labels_col]==study)&(df[categories_col]==category), effect_sizes_col].values[0]
            ci_low_error = df.loc[(df[labels_col]==study)&(df[categories_col]==category), "ci_low_error"].values[0]
            ci_high_error = df.loc[(df[labels_col]==study)&(df[categories_col]==category), "ci_high_error"].values[0]
            errors = np.array([ci_low_error, ci_high_error]).reshape(2, 1)
            # Plot each category with a different marker
            ax.errorbar(
                effect_size, idx+offset, xerr=[[ci_low_error], [ci_high_error]], fmt=markers[category], 
                color=colors[category], capsize=2, elinewidth=1, ms=5
                )

        # Place the study name above the group of categories
        ax.text(min(df[ci_low_col])-0.5, len(categories)*i+1+offset, study, ha='center', va='center', fontweight='normal', style="italic") # * * len(categories)
    
    ax.text(min(df[ci_low_col])-0.5, -2, "Genes", ha='center', va='center', fontweight='normal')
    
    # Invert the y-axis so the studies are from top to bottom
    ax.invert_yaxis()

    # Customize the left axis: remove yticks but keep the ytick labels for study names
    ax.set_yticks([])
    ax.tick_params(axis='y', which='both', length=0) 


    # Labels for axes
    ax.set_xlabel('Effect Size', fontsize=7)
    # ax.set_title('Forest Plot')
    # Remove top, right, and left spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # ax.legend(labels=["vline"]+categories, frameon=False, loc="best", ncol=1)

    # for idx, stat in enumerate(statistics):
    WRITE_HEADER_COUNT=0

    for i, study in enumerate(studies):
        offset = study_offset*i
        for j, category in enumerate(categories):
            for s, st_col in enumerate(stats_cols):
                idx = i * len(categories) + j*cat_offset
                stat = df.loc[(df[labels_col]==study)&(df[categories_col]==category), st_col].values[0]
                ax2.text(0.25+0.5*s, idx+offset, f"{stat}", va='center', ha="center")
                if WRITE_HEADER_COUNT<len(stats_cols):
                    ax2.text(0.25+0.5*s, -2, f"{st_col.lstrip('formatted_')}", va='center', ha="center")
                    WRITE_HEADER_COUNT+=1

    # Remove the spines and ticks for ax2
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    # ax2.set_title('Statistics')

    plt.tight_layout()
    return fig
