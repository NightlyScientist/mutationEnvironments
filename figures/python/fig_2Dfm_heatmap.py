# %%
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sb
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections import namedtuple
import importlib
from scipy.ndimage import zoom
from scipy.interpolate import griddata
import common.dataAPI as dataAPI
import common.datatables as dataAPIExtensions
import metrics.mutantFrequency as mutantFrequencyMethods

# reload modules (helpful for debugging)
importlib.reload(dataAPI)

workspace_paths = dataAPI.fetchWorkspaceEnv("../../")
ensemble_paths = dataAPI.fetchEnsemblePaths(workspace_paths["top_level_path"])

def showFilteredPaths(ensemble_paths, filter=None):
    for i, path in enumerate(ensemble_paths):
        if filter is not None and filter not in path.base:
            continue
        print(f"[{i}] ", path.base)

selector = 9

# %%
# pick one of the ensembles using 'example_ensemble_path'
example_ensemble_path = ensemble_paths[selector]

# ensemble table generated from ensemble directory and options.csv
ensemble_table = dataAPI.ensembleTableInfo(example_ensemble_path.base)

importlib.reload(dataAPIExtensions)
dataAPIExtensions.addMetricsColumns(ensemble_table)

# let's get a random simulation and plot the heatmap, along with some summary statistics
# Create a namedtuple from the column names of the ensemble_table
EnsembleRecord = namedtuple("EnsembleRecord", ensemble_table.columns)

# Example usage
random_index = np.random.randint(0, len(ensemble_table))
example_record = EnsembleRecord(*ensemble_table.iloc[random_index])
print(example_record)
print(os.listdir(example_record.path))

labels = {
    "time_extinction": r"$t_e / t_f$",
    "v_fraction": r"f$_{MT}$",
    "n_fraction": r"$N_0 / N_{total}$",
    "selection": r"$Selection$",
    "intensity": r"$Strength$",
    "mutation": r"$Mutation Rate$",
    "density": r"$Density$",
    "radius": r"$Radius$",
    "xi_m": r"$<\xi>$",
    "xi_var": r"$\xi_\sigma$",
}

#%%
class HeatmapSeries:
    def __init__(self):
        self.data = dict()

    def add_dataset(self, heatmap, selection, intensity):
        if selection not in self.data:
            self.data[selection] = dict()
        self.data[selection][intensity] = heatmap.mean(axis=1)[::-1]

    def draw(self, imgPath, rng_seed):
        for row in self.data:
            fig, ax = plt.subplots(figsize=(8, 6))
            for column in self.data[row]:
                ax.plot(self.data[row][column], label=column, marker="o")
            ax.legend()
            ax.set_ylabel(r"$f_m(y)$", fontsize=16)
            ax.set_xlabel("y", fontsize=16)
            fig.savefig(
                f"{imgPath}/hm_series_seed_{rng_seed}_row_{row}.png",
                transparent=True,
            )
            plt.close(fig)


def mutantHeatmapGridView(dataframe, paths, append_name="", all_intensity=False):
    df = dataframe
    s_values = np.sort(df.selection.unique())
    i_values = np.sort(df.intensity.unique())
    _groupby = ["selection", "intensity"]

    # select every third element in the list
    s_values = s_values[::3]
    if not all_intensity:
        i_values = i_values[::4]

    points = [(s, i) for s in s_values for i in i_values]
    grouped = df.groupby(_groupby)

    nrows, ncols = len(s_values), len(i_values)

    for ensemble in [True, False]:
        for rng_seed in np.sort(df.rngSeed.unique()):
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
            plt.subplots_adjust(wspace=0.01, hspace=0.01)

            heatmap_series = HeatmapSeries()

            for i, k in enumerate(points):
                sub_group = grouped.get_group(k)
                sub_group = sub_group[sub_group.rngSeed == rng_seed]
                _path = sub_group.path.values[0]
                lx, ly, trials = sub_group[["width", "height", "numberTrials"]].values[
                    0
                ]

                row = ncols - (i // ncols) - 1
                col = i % ncols
                if ensemble:
                    hm = mutantFrequencyMethods.heatmap(_path, lx, ly, trials)
                    im = axes[row, col].imshow(
                        hm, cmap="turbo", interpolation="nearest", vmin=0, vmax=1
                    )
                    axes[row, col].axis("off")
                    heatmap_series.add_dataset(hm, row, col)
                else:
                    figPath = os.path.join(_path, "snapshots_ID_1.png")
                    img = plt.imread(figPath)
                    axes[row, col].imshow(img)
                    axes[row, col].axis("off")
            if not all_intensity:
                plt.savefig(
                    f"{paths.img}/grid_ensemble_{ensemble}_seed_{rng_seed}_{append_name}.png",
                    transparent=True,
                )

            if ensemble:
                fig.subplots_adjust(right=0.85)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
                fig.colorbar(im, cax=cbar_ax)
                heatmap_series.draw(imgPath=paths.img, rng_seed=rng_seed)
            plt.close(fig)


def htsptSep(r: int, p: float | list | np.ndarray):
    return r * np.sqrt(-np.pi / np.log(1 - p))


def sBalance_alt(
    phi: float, L: int, radius: int, l: float, I: float, g: list[float] | float = 1
):
    r = radius
    _l = l * np.sqrt(2)
    k = 2 * np.sqrt(2) / np.sqrt(3)
    S = (4 * r * I) / (2 * I * r + _l + _l * I)
    return S**2 / k / 2


def sBalance(
    phi: float, L: int, radius: int, l: float, I: float, g: list[float] | float = 1
):
    r = radius 
    r = radius / 1.1547
    if phi > 0.1:
        _l = l * 2 * np.sqrt(2) / np.sqrt(3) * np.sqrt(2)
        #_l = l * 2 * np.sqrt(2) / np.sqrt(3)
    else:
        _l = l * np.sqrt(2)
    k =  1 * 1.1547
    S = (2 * r / _l) * (I / (1 + I))
    return S**2 / k


def partition(x, threshold=0.002):
    _x = np.sort(np.array(x))
    _indx = np.argsort(np.array(x))
    cuts = np.where(np.diff(_x) > threshold)
    x_partition = np.split(_x, cuts[0] + 1)
    indx_partition = np.split(_indx, cuts[0] + 1)
    return x_partition, indx_partition


def phaseDiagramView(df, save=False, paths=None):
    """from a dataframe, draw a heatmap for fm(s,v)"""
    sb.set_theme(style="white")
    lower_bound = 0.245
    upper_bound = 0.255
    selected_df = df[(df.v_fraction >= lower_bound) & (df.v_fraction <= upper_bound)]

    ps = df.pivot_table(index=x_axis, columns=y_axis, values=z_axis, aggfunc="mean")

    custom_norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5)

    fig, ax = plt.subplots(figsize=(10, 8), nrows=1)
    # fig, (ax, sub) = plt.subplots(figsize=(8, 9), nrows=2)

    hm_2Darray = ps.to_numpy()
    # uniform_env = hm_2Darray[:, 0]
    # print(uniform_env)

    # plt.imshow(hm_2Darray, cmap="coolwarm", norm=custom_norm, interpolation="nearest")
    sb.heatmap(
        hm_2Darray,
        fmt=".2f",
        annot=False,
        ax=ax,
        cbar_kws={"label": ""},
        cmap="PuOr",
        norm=custom_norm,
    )

    # Interpolate the 2D array using bilinear spline
    zoom_factor = 1  # Adjust the zoom factor as needed
    hm_2Darray_interpolated = zoom(hm_2Darray, zoom_factor, order=1)

    # Create a grid for the interpolated data
    x = np.linspace(0, hm_2Darray.shape[1] + 0.5, hm_2Darray_interpolated.shape[1])
    y = np.linspace(0, hm_2Darray.shape[0] + 0.5, hm_2Darray_interpolated.shape[0])
    X, Y = np.meshgrid(x, y)

    # Find the 0.25 contour
    contour = ax.contour(X, Y, hm_2Darray_interpolated, levels=[0.25], colors="green", linewidths=5)

    radius = np.mean(df.radius)
    density = np.mean(df.density)
    Lx = df.width.values[0]

    num_1 = 1 + 6 / 0.5  # range / step
    intensity = num_1 * np.linspace(np.min(df.intensity), np.max(df.intensity), 200)

    num_2 = 1 + 0.1 / 0.01  # range / step
    for f, c in [(sBalance, "black"), (sBalance_alt, "green")]:
        selection = num_2 * f(density, Lx, radius, htsptSep(radius, density), intensity)
        ax.plot(0.5 + intensity, 0.5 + 10 * selection, c=c, linewidth=10, alpha=0.5)
        break

    # sb.scatterplot(
    #    data=selected_df,
    #    x="intensity",
    #    y="selection",
    #    hue="v_fraction",
    #    ax=sub,
    # )

    # intensity = intensity / num_1
    # for f, c in [(sBalance, "black"), (sBalance_alt, "green")]:
    #    selection = f(density, Lx, radius, htsptSep(radius, density), intensity)
    #    sub.scatter(intensity, selection, c=c, s=10)

    # sub.set_ylim(0, 0.1)

    # ax.set_xlim(0, 12)

    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel("", fontsize=20)
    cbar.ax.tick_params(labelsize=22)

    # . Reduce number of tick labels by 2
    ax.set_xticklabels([f"{x / 2 - 0.25:.1f}" for x in ax.get_xticks()])
    ax.set_yticklabels([f"{x * 0.01 - 0.005:.2f}" for x in ax.get_yticks()])
    ax.set_xticks(ax.get_xticks()[::2])
    ax.set_yticks(ax.get_yticks()[::2])

    # . increase tick label size
    ax.tick_params(labelsize=18)

    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_ylabel("")

    if save:
        fig.savefig(f"{paths.img}/phase_diagram_revised.pdf", transparent=True)
        plt.close(fig)
    else:
        plt.show()

df = ensemble_table
x_axis, y_axis, z_axis = "selection", "intensity", "v_fraction"
label = labels[z_axis]
phaseDiagramView(df, save=False, paths=example_ensemble_path)
phaseDiagramView(df, save=True, paths=example_ensemble_path)

# %%
# > sequence of images
sb.set_theme()

s_values = np.sort(df.selection.unique())
i_values = np.sort(df.intensity.unique())
_groupby = ["selection", "intensity"]

# select every third element in the list
s_values = s_values[::3]
i_values = i_values[::4]

points = [(s, i) for s in s_values for i in i_values]
grouped = df.groupby(_groupby)

nrows, ncols = len(s_values), len(i_values)

# draw the entire grid view for the heatmps, snapshots, and draw x-averaged f_m
mutantHeatmapGridView(dataframe=df, paths=example_ensemble_path, all_intensity=False)

# %%
# > x-averaged f_m for each rng seed and selection
s_values = np.sort(df.selection.unique())
s_values = s_values[::3]
i_values = np.sort(df.intensity.unique())

# set the color map and normalizer
_min, _max = i_values.min(), i_values.max()
norm = colors.TwoSlopeNorm(vmin=_min, vcenter=(_max + _min) / 2, vmax=_max)
cmap = mpl.cm.get_cmap("plasma")

sb.set_theme(style="ticks")

for rng_seed in np.sort(df.rngSeed.unique()):
    for row, s_value in enumerate(s_values):
        heatmap_series = HeatmapSeries()
        sub_dataframe = df[(df.rngSeed == rng_seed) & (df.selection == s_value)]

        for i, record in sub_dataframe.iterrows():
            hm = mutantFrequencyMethods.heatmap(
                record.path, record.width, record.height, record.numberTrials
            )
            heatmap_series.add_dataset(hm, record.selection, record.intensity)

        fig, ax = plt.subplots(figsize=(8, 6))

        for selection in np.sort(list(heatmap_series.data.keys())):
            for intensity in np.sort(list(heatmap_series.data[selection])):
                label = f"{intensity}"
                ax.plot(
                    heatmap_series.data[selection][intensity],
                    label=label,
                    marker="o",
                    color=cmap(norm(intensity)),
                )
        # Add inset colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.046, pad=0.04)
        cbar.set_label("Intensity", fontsize=16)

        # ax.legend(title="Intensity")
        ax.set_ylabel(r"$f_m(y)$", fontsize=16)
        ax.set_xlabel("y", fontsize=16)
        ax.set_ylim(0, 0.75)
        fig.suptitle(f"selection: {selection}")
        fig.savefig(
            f"{example_ensemble_path.img}/hm_series_seed_{rng_seed}_row_{row}.png",
            transparent=True,
        )
        plt.close(fig)
    # break
# %%
df = ensemble_table
# > x-averaged f_m for each selection and averaged across each rng seed
s_values = np.sort(df.selection.unique())
s_values = s_values[::3]
i_values = np.sort(df.intensity.unique())
i_values = i_values[::4]

# set the color map and normalizer
_min, _max = i_values.min(), i_values.max()
norm = colors.TwoSlopeNorm(vmin=_min, vcenter=(_max + _min) / 2, vmax=_max)
cmap = mpl.cm.get_cmap("winter")

sb.set_theme(style="ticks")
matrix_size = (len(df.rngSeed.unique()), df.height.values[0])

for row, s_value in enumerate(s_values):
    if row != 1:
        continue

    fig, ax = plt.subplots(figsize=(8, 6))

    for i, i_value in enumerate(i_values):
        sub_dataframe = df[(df.selection == s_value) & (df.intensity == i_value)]
        heatmap_series = np.zeros(matrix_size)

        # iterate through each record associated with different rng seeds
        indices_skipped = []
        for j, record in sub_dataframe.iterrows():
            if not os.path.exists(f"{record.path}/heatmap_ID3.arrow"):
                indices_skipped.append(record.rngSeed)
                continue
            hm = mutantFrequencyMethods.heatmap(
                record.path, record.width, record.height, record.numberTrials
            )
            # heatmap_series.add_dataset(hm, record.selection, record.intensity)
            heatmap_series[record.rngSeed - 1, :] = hm.mean(axis=1)[::-1]
        
        if len(indices_skipped) > 0:
            print(f"Skipped indices: {len(indices_skipped)}")

        # average across each rng seed and create uncertainty band plot
        # Skip indices in mean and std from indices_skipped
        valid_indices = [i for i in range(matrix_size[0]) if i not in indices_skipped]
        mean_heatmap = heatmap_series[valid_indices, :].mean(axis=0)
        std_heatmap = heatmap_series[valid_indices, :].std(axis=0)
        #mean_heatmap = heatmap_series.mean(axis=0)
        #std_heatmap = heatmap_series.std(axis=0)

        ax.plot(
            mean_heatmap,
            label=f"{i_value}",
            marker="o",
            color=cmap(norm(i_value)),
        )
        ax.fill_between(
            range(len(mean_heatmap)),
            mean_heatmap - std_heatmap,
            mean_heatmap + std_heatmap,
            color=cmap(norm(i_value)),
            alpha=0.3,
        )
    
    ax.axvline(x=20, color='red', linestyle='--', linewidth=2)

    # Add inset colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.046, pad=0.04)
    cbar.set_label("Intensity", fontsize=16)

    ax.legend(title="Intensity")
    ax.set_ylabel(r"$f_m(y)$", fontsize=16)
    ax.set_xlabel("y", fontsize=16)
    #ax.set_ylim(0, 0.75)
    fig.suptitle(f"selection: {s_value}")
    #fig.savefig(
    #    f"{example_ensemble_path.img}/hm_series_row_{row}.png",
    #    transparent=True,
    #)
    #plt.close(fig)
    print(s_value)
    plt.show()
    # break


# %%
for ensemble in [True, False]:
    for rng_seed in df.rngSeed.unique():
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
        plt.subplots_adjust(wspace=0.01, hspace=0.01)

        for i, k in enumerate(points):
            sub_group = grouped.get_group(k)
            sub_group = sub_group[sub_group.rngSeed == rng_seed]
            _path = sub_group.path.values[0]
            lx, ly, trials = sub_group[["width", "height", "numberTrials"]].values[0]

            row = ncols - (i // ncols) - 1
            col = i % ncols
            if ensemble:
                hm = mutantFrequencyMethods.heatmap(_path, lx, ly, trials)
                im = axes[row, col].imshow(
                    hm, cmap="turbo", interpolation="nearest", vmin=0, vmax=1
                )
                axes[row, col].axis("off")
            else:
                figPath = os.path.join(_path, "snapshots_ID_1.png")
                img = plt.imread(figPath)
                axes[row, col].imshow(img)
                axes[row, col].axis("off")
        plt.savefig(
            f"{imgPath}/grid_ensemble_{ensemble}_seed_{rng_seed}.png", transparent=True
        )

        if ensemble:
            fig.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
            fig.colorbar(im, cax=cbar_ax)
        plt.close(fig)

# %%
# >try to estimate the rate of change of Pm(h)
random_path = np.random.choice(
    df.path[(df.selection > 0.01) & (df.selection < 0.5) & (df.intensity > 4)].values, 1
)[0]
hm = mutantFrequencyMethods.heatmap(random_path, 2000, 1000, 100)

fig, ax = plt.subplots(figsize=(12, 6), ncols=2)
ax[0].imshow(hm, cmap="turbo", interpolation="nearest", vmin=0, vmax=0.5)

# mean value along the x-axis
mean_hm = hm.mean(axis=1)[::-1]
ax[1].plot(mean_hm)
ax[1].set_ylim(0, 0.55)
# ax[1].set_yscale("log")

# >sequence of images
# this section will produce a sequence of images for a given range of mutation rates
sb.set_theme()

s_values = np.sort(df.compensation.unique())
i_values = np.sort(df.mutation.unique())
_groupby = ["compensation", "mutation"]

# select every third element in the list
# s_values = s_values[3:6]
s_values = s_values[::3]
i_values = i_values[::4]
# s_values = [0.3, 0.2, 0.15, 0.05]
# i_values = [0, 0.03, 0.05, 0.1]

points = [(s, i) for s in s_values for i in i_values]
grouped = df.groupby(_groupby)

nrows, ncols = len(s_values), len(i_values)

for ensemble in [True, False]:
    for rng_seed in df.rngSeed.unique():
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
        plt.subplots_adjust(wspace=0.01, hspace=0.01)

        for i, k in enumerate(points):
            sub_group = grouped.get_group(k)
            sub_group = sub_group[sub_group.rngSeed == rng_seed]
            _path = sub_group.path.values[0]
            lx, ly, trials = sub_group[["width", "height", "numberTrials"]].values[0]

            row = ncols - (i // ncols) - 1
            col = i % ncols
            if ensemble:
                hm = mutantFrequencyMethods.heatmap(_path, lx, ly, trials)
                im = axes[row, col].imshow(
                    hm, cmap="turbo", interpolation="nearest", vmin=0, vmax=1
                )
                axes[row, col].axis("off")
            else:
                figPath = os.path.join(_path, "snapshots_ID_3.png")
                img = plt.imread(figPath)
                axes[row, col].imshow(img)
                axes[row, col].axis("off")
        plt.savefig(
            f"{imgPath}/grid_ensemble_{ensemble}_seed_{rng_seed}.png", transparent=True
        )

        if ensemble:
            fig.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
            fig.colorbar(im, cax=cbar_ax)
        break

# %%
x_axis = "mutation"
y_axis = "compensation"
# z_axis = "n_fraction"

# > phase diagram (heatmap)
# this generats a heatmap for fm(s,v)
sb.set_theme(style="white")

for z_axis in ["n_fraction", "v_fraction"]:
    label = labels[z_axis]

    ps = df.pivot_table(index=y_axis, columns=x_axis, values=z_axis, aggfunc="mean")

    custom_norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax = sb.heatmap(
        ps,
        fmt=".2f",
        annot=False,
        ax=ax,
        cbar_kws={"label": ""},
        cmap="vlag",
        # norm=custom_norm,
    )

    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel(label, fontsize=20)

    # for i, s in zip(i_values, s_values):
    #    ax.scatter((i + 0.005) * (10 / 0.1), (s + 0.025) * (6 / 0.3), color='black', alpha = 0.9, marker='x', s=400)

    ax.invert_yaxis()
    ax.tick_params(labelsize=12)
    ax.set_xlabel("$\mu$")
    ax.set_ylabel("C")
    fig.savefig(
        f"{imgPath}/phase_diagram_{x_axis}_{y_axis}_{z_axis}.png", transparent=True
    )
