# This script will generate a heatmap of the mutational frequency for different values of mutation and compensation. and a sample plot of the mutational frequency.
# %%
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib as mpl
import importlib
from matplotlib import colors
from collections import namedtuple
import common.dataAPI as dataAPI
import common.datatables as dataAPIExtensions
import heatmaps.mutantFrequency as mutantFrequencyMethods
import os

# some of the numpy erors get annoying
np.seterr(divide="ignore", invalid="ignore")

# reload modules (helpful for debugging)
importlib.reload(dataAPI)

workspace_paths = dataAPI.fetchWorkspaceEnv("../../")
ensemble_paths = dataAPI.fetchEnsemblePaths(workspace_paths["top_level_path"])

# pick one of the ensembles using 'example_ensemble_path'
example_ensemble_path = ensemble_paths[0]

# %%
# ensemble table generated from ensemble directory and options.csv
ensemble_table = dataAPI.ensembleTableInfo(example_ensemble_path.base)

# %%
importlib.reload(dataAPIExtensions)
dataAPIExtensions.addMetricsColumns(ensemble_table)

# %%
# let's get a random simulation and plot the heatmap, along with some summary statistics
# Create a namedtuple from the column names of the ensemble_table
EnsembleRecord = namedtuple("EnsembleRecord", ensemble_table.columns)

# Example usage
random_index = np.random.randint(0, len(ensemble_table))
example_record = EnsembleRecord(*ensemble_table.iloc[random_index])
print(example_record)
print(os.listdir(example_record.path))

# create side by side plots of heatmap and species snapshots
if os.path.exists(example_record.path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Plot the species snapshot
    img_path = os.path.join(example_record.path, "snapshots_ID_3.png")
    img = plt.imread(img_path)
    axes[0].imshow(img)
    axes[0].axis("off")

    # Plot the heatmap
    heatmap_grid = mutantFrequencyMethods.heatmap(
        example_record.path,
        example_record.width,
        example_record.height,
        example_record.numberTrials,
    )
    # plt.figure(figsize=(10, 8))
    im = sb.heatmap(
        heatmap_grid,
        cmap="turbo",
        vmin=0,
        vmax=0.5,
        cbar=True,
        xticklabels=False,
        yticklabels=False,
        ax=axes[1],
    )
    # plt.title("Heatmap from heatmap_ID3.arrow")
    axes[1].axis("off")
    plt.show()

# %%
# find the x-averaged mutational frequency for each selection and intensity
sb.set_theme(style="ticks")

df = ensemble_table
s_values = np.sort(df.selection.unique())
i_values = np.sort(df.intensity.unique())

_min, _max = i_values.min(), i_values.max()
norm = colors.TwoSlopeNorm(vmin=_min, vcenter=(_max + _min) / 2, vmax=_max)

cmap = mpl.cm.get_cmap("plasma")
fig, axes = plt.subplots(nrows=len(s_values), ncols=1, figsize=(8, 3.5*len(s_values)))

for i, sel in enumerate(s_values):
    _df = df[df.selection == sel]
    grouped = _df.groupby("intensity")

    for intensity in grouped.groups.keys():
        sub_group = grouped.get_group(intensity)

        lx, ly, trials = sub_group[["width", "height", "numberTrials"]].values[0]
        mean_hm = np.zeros((ly,))
        for seed in sub_group.rngSeed.values:
            sub_path = sub_group.path.values[0]
            hm = mutantFrequencyMethods.heatmap(sub_path, lx, ly, trials)
            tmp = hm.mean(axis=1)[::-1]
            mean_hm += tmp

        axes[i].plot(mean_hm / len(sub_group.rngSeed.values), color=cmap(norm(intensity)), linestyle="--", label=f"I: {intensity}")

    axes[i].set_ylim(-0.05, 1)
    axes[i].set_title(f"S: {sel}", fontsize=10)
plt.legend(title="Intensity", bbox_to_anchor=(1.05, 1), loc='upper left')
