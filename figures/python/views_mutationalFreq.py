# %%
import pandas as pd
import pyarrow.feather as feather
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.gridspec as gridspec
from matplotlib import colors
from collections import namedtuple
import matplotlib as mpl
import importlib
import common.dataAPI as dataAPI
import common.datatables as dataAPIExtensions
import metrics.mutantFrequency as mutantFrequencyMethods

# reload modules (helpful for debugging)
importlib.reload(dataAPI)

workspace_paths = dataAPI.fetchWorkspaceEnv("../../")
ensemble_paths = dataAPI.fetchEnsemblePaths(workspace_paths["top_level_path"])

# %%
# pick one of the ensembles using 'example_ensemble_path'
example_ensemble_path = ensemble_paths[0]

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

def fetchMutantFrequency(path, skip=0, file="mutationalFreq.arrow"):
    df = feather.read_feather(os.path.join(path, file))
    df["time"] = df.index
    df["mutationalFreq"] = 1 - df.mutationalFreq

    t = df.mutationalFreq.values
    first_zero_index = np.where(t == 0)[0]
    if len(first_zero_index) > 0:
        first_zero_value = first_zero_index[0]
    else:
        first_zero_value = len(t)

    df["f_m"] = np.sum(t[0:first_zero_value]) / first_zero_value
    return df

def parameterIter(df, parameter, record_opts=["mutation", "compensation"]):
    for i, record in df.iterrows():
        _df = fetchMutantFrequency(record.path)
        f_m = _df[parameter].values[0]
        m = record[record_opts[0]]
        c = record[record_opts[1]]
        yield (m, c, f_m)

# this will generate a sample plot of the mutational frequency
sns.set_theme(style="white")

df = ensemble_table

gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[:, 2:])
ax2 = plt.subplot(gs[0:2, 0:2])
ax3 = plt.subplot(gs[2::, 0:2])

samplePath = np.random.choice(df.path.values)
# samplePath = nprandom.choice(df[df.mutation == 0.05].path.values)

# image of species configuration in the simulation
image_path = os.path.join(example_record.path, "snapshots_ID_3.png")
image = plt.imread(image_path)
ax1.imshow(image)
ax1.axis("off")

axes = (ax2, ax3)
sampe_df = df[df["mutation"] == 0.05]

cmap = mpl.colormaps["coolwarm"]

record_opts = ["intensity", "selection"]
vmin = sampe_df[record_opts[0]].min()
vmax = sampe_df[record_opts[1]].max()
vcenter = (vmax + vmin) / 2
custom_norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

for i, row in sampe_df.iterrows():
    _df = fetchMutantFrequency(row.path)
    _second_parameter = row[record_opts[1]]
    x, y, c = _df.time, _df.mutationalFreq, cmap(custom_norm(_second_parameter))
    sns.lineplot(x=x, y=y, ax=axes[0], c=c, linewidth=2.5)

axes[0].set_xscale("log")

# heatmap grid of f_m
x, y, z = zip(*parameterIter(df, "f_m", record_opts=record_opts))
heatmap_data = pd.DataFrame({"x": x, "y": y, "z": z})
ps = heatmap_data.pivot_table(index="y", columns="x", values="z", aggfunc="mean")

sns.heatmap(
    ps,
    fmt=".2f",
    annot=False,
    ax=axes[1],
    cbar_kws={"label": "$f_m$"},
    cmap="coolwarm",
)

axes[1].set_xlabel("intensity")
axes[1].set_ylabel("selection")
axes[1].invert_yaxis()

# Add rectangular outline to heatmap at x-value 0.05
rect = plt.Rectangle((10, 0), 1, 10, edgecolor="green", facecolor="none", linewidth=2)
axes[1].add_patch(rect)

plt.tight_layout(h_pad=0.5, w_pad=0.15)
fig = plt.gcf()

imgs_path = f"{example_ensemble_path.img}/mutational_freq"
if not os.path.exists(imgs_path):
    os.makedirs(imgs_path, exist_ok=True)

fig_path = os.path.join(imgs_path, os.path.basename(os.path.normpath(example_record.path))) + ".png"
if os.path.exists(imgs_path):
    fig.savefig(fig_path, transparent=True)
plt.show()

# %%
