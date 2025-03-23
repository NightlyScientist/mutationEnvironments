# %%
import matplotlib.pyplot as plt
import os
from collections import namedtuple
import importlib
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import common.dataAPI as dataAPI
import common.datatables as dataAPIExtensions
import metrics.mutantFrequency as mutantFrequencyMethods
import argparse

parser = argparse.ArgumentParser(
    description="Generate 2D heatmap images for the mutant frequency"
)
parser.add_argument("--path", type=str, help="path to simulation directory")
args = parser.parse_args()


# reload modules (helpful for debugging)
importlib.reload(dataAPI)

ensemble_table = dataAPI.ensembleTableInfo(args.path)

# %%
importlib.reload(dataAPIExtensions)
dataAPIExtensions.addMetricsColumns(ensemble_table)

# let's get a random simulation and plot the heatmap, along with some summary statistics
# Create a namedtuple from the column names of the ensemble_table
EnsembleRecord = namedtuple("EnsembleRecord", ensemble_table.columns)

# iterate over the ensemble_table and print the first record
for i, _ in ensemble_table.iterrows():
    example_record = EnsembleRecord(*ensemble_table.iloc[i])
    print(example_record.path)

    hm = mutantFrequencyMethods.heatmap(
        example_record.path,
        example_record.width,
        example_record.height,
        example_record.numberTrials,
    )

    #fig, (ax,mid, right) = plt.subplots(1, 3, figsize=(20, 8))
    fig, ax = plt.subplots(1, 1, figsize=(4, 8))
    hm = hm[::-1]
    ax.imshow(hm, cmap="jet", aspect="auto")
    #cbar = fig.colorbar(
    #    ax.imshow(hm, cmap="jet", aspect="auto"), ax=ax, orientation="vertical"
    #)
    #cbar.set_label("Mutant Probability")
    #cbar.ax.tick_params(labelsize=10)
    #cbar.mappable.set_clim(0, 1)
    #fig.subplots_adjust(right=0.9)

    ax.invert_yaxis()

    ax.tick_params(axis='both', which='major', labelsize=26)
    ax.set_xlabel("x", fontsize=30)
    ax.set_ylabel("y", fontsize=30)
    #ax.set_title("Mutant Probability Heatmap")
    plt.tight_layout(pad=0)

    paths = dataAPI.setOutputPaths(example_record.path)
    plt.savefig(os.path.join(paths.img, "mutantFreqImg.png"))
    plt.close()


    fig, ax = plt.subplots(1, 1, figsize=(9, 8))
    ax.plot(hm.mean(axis=1))
    ax.set_xlabel("y", fontsize=20)
    ax.set_ylabel(r"$f_m(y)$", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.savefig(os.path.join(paths.img, "f_m_y.png"))
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    img_path = os.path.join(example_record.path, "snapshots_ID_4.png")
    img = plt.imread(img_path)
    ax.imshow(img)
    ax.axis('off')
    plt.tight_layout(pad=0)
    plt.savefig(os.path.join(paths.img, "mutations_snapshots.png"))
    plt.close()