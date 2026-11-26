# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: brainvar2_mbv
#     language: python
#     name: python3
# ---

# %% [markdown]
# # BrainVar2 MBV figures
#
# This notebook reads the lossless table produced by `collect_mbv_reports.py` and creates review figures. It does not assign RNA alignments to donors. Pairing decisions are made manually from these figures and recorded in the project PowerPoint/review record.

# %%
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from tqdm.auto import tqdm

DATA = Path('../results/mbv_all.tsv.gz')
FIGURE_DIR = Path('../results/figures')
TOP_LABELS_PER_METRIC = 3
MAKE_COHORT_HEATMAPS = True

FIGURE_DIR.mkdir(parents=True, exist_ok=True)
sns.set_theme(style='whitegrid', context='notebook')

# %%
mbv = pd.read_csv(DATA, sep='\t')
required = {
    'rna_alignment', 'candidate_donor',
    'perc_het_consistent', 'perc_hom_consistent',
    'n_het_covered', 'n_hom_covered',
}
missing = required.difference(mbv.columns)
if missing:
    raise ValueError(f'Missing required columns: {sorted(missing)}')

print(f"RNA alignments: {mbv['rna_alignment'].nunique():,}")
print(f"Candidate donors: {mbv['candidate_donor'].nunique():,}")
print(f"Candidate comparisons: {len(mbv):,}")


# %%
def plot_alignment(data, alignment):
    data = data.copy()
    top_het = data.nlargest(TOP_LABELS_PER_METRIC, 'perc_het_consistent')
    top_hom = data.nlargest(TOP_LABELS_PER_METRIC, 'perc_hom_consistent')
    labels = set(top_het['candidate_donor']) | set(top_hom['candidate_donor'])

    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    ax.scatter(
        data['perc_het_consistent'],
        data['perc_hom_consistent'],
        s=24, color='#7A7A7A', alpha=0.55, linewidth=0,
    )
    highlighted = data[data['candidate_donor'].isin(labels)]
    ax.scatter(
        highlighted['perc_het_consistent'],
        highlighted['perc_hom_consistent'],
        s=38, color='#B42318', alpha=0.9, linewidth=0,
    )
    texts = [
        ax.text(row.perc_het_consistent, row.perc_hom_consistent, row.candidate_donor, fontsize=8)
        for row in highlighted.itertuples()
    ]
    adjust_text(texts, ax=ax, arrowprops={'arrowstyle': '-', 'color': '#555555', 'lw': 0.6})
    ax.set(
        title=alignment,
        xlabel='Consistent covered heterozygous genotypes',
        ylabel='Consistent covered homozygous genotypes',
        xlim=(0, 1.01), ylim=(0, 1.01),
    )
    return fig

per_alignment_dir = FIGURE_DIR / 'per_alignment'
per_alignment_dir.mkdir(parents=True, exist_ok=True)
for alignment, data in tqdm(mbv.groupby('rna_alignment', sort=True)):
    fig = plot_alignment(data, alignment)
    safe_name = alignment.replace('/', '_')
    fig.savefig(per_alignment_dir / f'{safe_name}.svg')
    fig.savefig(per_alignment_dir / f'{safe_name}.png', dpi=200)
    plt.close(fig)

# %%
if MAKE_COHORT_HEATMAPS:
    for metric, title in [
        ('perc_het_consistent', 'Heterozygous concordance'),
        ('perc_hom_consistent', 'Homozygous concordance'),
    ]:
        matrix = mbv.pivot(index='rna_alignment', columns='candidate_donor', values=metric)
        width = max(10, matrix.shape[1] * 0.08)
        height = max(8, matrix.shape[0] * 0.08)
        fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
        sns.heatmap(matrix, cmap='viridis', vmin=0, vmax=1, xticklabels=False, yticklabels=False, ax=ax)
        ax.set(title=title, xlabel='Candidate donor', ylabel='RNA alignment')
        fig.savefig(FIGURE_DIR / f'{metric}.svg')
        fig.savefig(FIGURE_DIR / f'{metric}.png', dpi=200)
        plt.close(fig)

# %% [markdown]
# ## Manual review
#
# The generated figures are inputs to manual review. Organize the per-alignment figures in the controlled project PowerPoint, compare the strongest heterozygous and homozygous concordance signals, consult sample metadata where needed, and record the selected pairing or unresolved status outside this notebook. Do not infer pairings from filenames alone.
