# %%
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sb
from matplotlib import colors
import pyarrow.feather as feather
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from functools import partial

basePath = input("base path")

imgPath = os.path.join(basePath, "images")
if not os.path.exists(imgPath):
    os.makedirs(imgPath)


# %%
def phaseSpaceTable(basePath):
    dictList = []
    for root, dirs, _ in os.walk(basePath):
        for directory in dirs:
            if not "env" in directory:
                continue
            subpath = os.path.join(root, directory)

            opts = pd.read_csv(
                os.path.join(subpath, "inputOpts.csv"), sep="\t", header=0
            ).to_dict(orient="index")[0]

            opts["path"] = subpath

            dictList.append(opts)
    return pd.DataFrame.from_dict(dictList)


def partition(x, threshold=0.002):
    _x = np.sort(np.array(x))
    _indx = np.argsort(np.array(x))
    cuts = np.where(np.diff(_x) > threshold)
    x_partition = np.split(_x, cuts[0] + 1)
    indx_partition = np.split(_indx, cuts[0] + 1)
    return x_partition, indx_partition


# doc genreate colorbar by itself and save to file
def colorbar_plot(cmap, custom_norm):
    fig, ax = plt.subplots(figsize=(6, 4))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=custom_norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, orientation="horizontal", label="", location="top")
    fig.savefig(f"{imgPath}/colorbar_{cmap}.svg", transparent=True)


# %%
df = phaseSpaceTable(basePath)
x_axis = "selection"
y_axis = "intensity"

sb.set_theme()

# %%
# > animate along second axis (intensity) / each landscape is a unique rng seed
def lineageHeatmap(i, trials=1):
    f = feather.read_feather(
        f"{paths[i]}/lineage_tracks_heatmap.arrow"
    ).lineage_tracks.values
    hm = f.reshape((lx, ly))
    axes[0].imshow(hm, vmin=0.001, vmax=0.3, origin="lower", cmap="Reds")

    f = feather.read_feather(f"{paths[i]}/heatmap_ID3.arrow").heatmap_ID3.values
    hm = (f.reshape((lx, ly)) / trials) - 1
    axes[1].imshow(hm, vmin=0, vmax=1, origin="lower", cmap="turbo")

landscape = df[df.rngSeed == 1]
selections = np.sort(landscape.selection.unique())

for i, _ in enumerate(selections):
    section = landscape[landscape.selection == selections[i]]
    lx, ly, trials = section[["width", "height", "numberTrials"]].values[0]

    paths = section.path.values[np.argsort(section.intensity.values)]

    fig, axes = plt.subplots(nrows=2, figsize=(8, 8))
    fig.tight_layout()

    for ax in axes:
        ax.axis("off")
    plt.subplots_adjust(wspace=0.01, hspace=0.01)

    ani = animation.FuncAnimation(
        fig,
        partial(lineageHeatmap, trials=trials),
        interval=50,
        blit=False,
        frames=range(len(paths)),
        repeat_delay=100,
    )

    if not os.path.exists(imgPath + "/heatmaps"):
        os.makedirs(imgPath + "/heatmaps")

    ani.save(f"{imgPath}/heatmaps/lineage_heatmap_{selections[i]}.gif", writer="imagemagick", fps=8)

#%%

