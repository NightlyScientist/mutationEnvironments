# %%
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sb
from matplotlib import colors
import pyarrow.feather as feather
from mpl_toolkits.axes_grid1 import make_axes_locatable

basePath = input("base path")

imgPath = os.path.join(basePath, "images")
if not os.path.exists(imgPath):
    os.makedirs(imgPath)

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

sb.set_theme(style="white")


# %%
def speciesFraction(df, num=2, mN=0, mLx=0):
    n = [f"n_{i}" for i in range(1, num + 1)]
    v = [f"v_{i}" for i in range(1, num + 1)]
    n_fraction = ((df[n[-1]] - mLx) / (df[n].sum(axis=1))).mean()
    v_fraction = ((df[v[-1]] - mN) / (df[v].sum(axis=1))).mean()
    return n_fraction, v_fraction


def phaseSpaceTable(path: str):
    dictList = []
    for root, dirs, _ in os.walk(path):
        for directory in dirs:
            if not "env" in directory:
                continue
            subpath = os.path.join(root, directory)

            opts = pd.read_csv(
                os.path.join(subpath, "inputOpts.csv"), sep="\t", header=0
            ).to_dict(orient="index")[0]

            # .if there is mutations, get background noise as mu * N
            lx, ly = opts["width"], opts["height"]
            mN = opts["mutation"] * lx * ly
            mLx = opts["mutation"] * lx

            data = pd.read_csv(os.path.join(subpath, "table.csv"), header=0)
            n_fraction, v_fraction = speciesFraction(data, num=2, mN=mN, mLx=mLx)

            opts["path"] = subpath
            opts["time_extinction"] = (data.time_extinction / data.time)[
                data.time_extinction > 0
            ].mean()
            opts["xi_m"] = data.xi_m.mean()
            opts["xi_var"] = np.sqrt(data.xi_var.mean())
            opts["n_fraction"] = n_fraction
            opts["v_fraction"] = v_fraction

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
    fig, ax = plt.subplots(figsize=(6, 6))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=custom_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", label="", location="right")
    cbar.ax.tick_params(labelsize=18)

    fig.savefig(f"{imgPath}/colorbar_{cmap}.svg", transparent=True)


# %%
colorbar_plot("coolwarm", colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5))
colorbar_plot("turbo", colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5))

# %%
df = phaseSpaceTable(basePath)
x_axis = "selection"
y_axis = "intensity"
z_axis = "v_fraction"
label = labels[z_axis]


# %%
# > phase diagram (heatmap)
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
    if phi > 0.1:
        _l = l * 2 * np.sqrt(2) / np.sqrt(3)
    else:
        _l = l * np.sqrt(2)
    k = 2 * 1.1547
    S = (2 * r / _l) * (I / (1 + I))
    return S**2 / k


def phase_diagram(save=False, base_path=None):
    lower_bound = 0.245
    upper_bound = 0.255
    selected_df = df[(df.v_fraction >= lower_bound) & (df.v_fraction <= upper_bound)]

    ps = df.pivot_table(index=x_axis, columns=y_axis, values=z_axis, aggfunc="mean")

    custom_norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5)

    fig, ax = plt.subplots(figsize=(10, 8), nrows=1)
    # fig, (ax, sub) = plt.subplots(figsize=(8, 9), nrows=2)

    sb.heatmap(
        data=ps,
        fmt=".2f",
        annot=False,
        ax=ax,
        cbar_kws={"label": ""},
        cmap="coolwarm",
        norm=custom_norm,
    )

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

    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel("", fontsize=20)
    cbar.ax.tick_params(labelsize=22)

    # . Reduce number of tick labels by 2
    ax.set_xticks(ax.get_xticks()[::2])
    ax.set_yticks(ax.get_yticks()[::2])

    # . increase tick label size
    ax.tick_params(labelsize=18)

    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_ylabel("")

    if save:
        imgPath = os.path.join(base_path, "images")
        if not os.path.exists(imgPath):
            os.makedirs(imgPath)

        fig.savefig(f"{imgPath}/phase_diagram_alt.pdf", transparent=True)
        plt.close(fig)
    else:
        plt.show()


#%%
phase_diagram(save=True, base_path=basePath)

# %%
top_directory = (
    "workspace/experiments"
)

directories = [
    d
    for d in os.listdir(top_directory)
    if os.path.isdir(os.path.join(top_directory, d))
]
directories = [d for d in directories if "phase_diagram_no_mutations" in d]

for dir in directories:
    _path = os.path.join(top_directory, dir)
    df = phaseSpaceTable(_path)
    x_axis = "selection"
    y_axis = "intensity"
    z_axis = "v_fraction"
    label = labels[z_axis]
    phase_diagram(save=True, base_path=_path)


# %%
# >sequence of images
def mutantHeatmap(path: str, lx, ly, trials):
    f = feather.read_feather(f"{path}/heatmap_ID3.arrow").heatmap_ID3.values
    hm = (f.reshape((lx, ly)) / trials) - 1
    return np.flip(hm, axis=0)


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
                hm = mutantHeatmap(_path, lx, ly, trials)
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
# >  maximal mutatant survival
max_surv_path = ""
df = phaseSpaceTable(max_surv_path)

imgPath = os.path.join(basePath, "images")
if not os.path.exists(imgPath):
    os.makedirs(imgPath)
#%%
limiter = 0.09
observable = "n_fraction"
cmap = mpl.colormaps.get_cmap("plasma")

fig, ax = plt.subplots(ncols=1, figsize=(6, 6))

for name, group in df.groupby("selection"):
    S = float(name)

    indx_partition = partition(group.density)[1]
    radius = group.radius.values[0]
    diameter = 2 * radius

    # mean phi for each partition
    density = np.array([group.density.values[p].mean() for p in indx_partition])
    cs = cmap((density - density.min()) / (density.max() - density.min()))

    # convert phi to lambda
    hl = np.array([htsptSep(radius, d) for d in density])

    # fetch observable values
    obs = np.array([group[observable].values[p].mean() for p in indx_partition])
    obs_std = [
        group[observable].values[p].std() / np.sqrt(p.size) for p in indx_partition
    ]
    ax.errorbar(
        np.sqrt(2) * hl / diameter,
        #density,
        obs,
        yerr=obs_std,
        fmt="o",
        color=cmap(S / 0.10),
    )
    # ax.plot(density / diameter, v_fraction, linestyle="-", color=cmap(S / 0.10))

    r = group.radius.values[0]
    #r = group.radius.values[0] * np.sqrt(3)
    #intensity = group.intensity.values[0]
    #s = group.selection.values[0]
    #z_c = 2 * r * intensity / (1 + intensity) / s

    # .percolation of overlapping circles
    for z in [htsptSep(r, 0.5), htsptSep(r, 0.68)]:
        ax.axvline(x= np.sqrt(2) * z / diameter, color="black", linestyle="--")
        break

    # add phi = 0.68 line
    #ax.axvline(x= 0.68, color="black", linestyle="--")
    #ax.axvline(x= 0.54, color="black", linestyle="--")
    #ax.axvline(x= 0.4, color="black", linestyle="--")

v_min, v_max = [0, 0.1]
norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
divider = make_axes_locatable(ax)
colorbar_axes = divider.append_axes("right", size="3%", pad=0.05)
cb = plt.colorbar(sm, cax=colorbar_axes, orientation="vertical", label="")
cb.set_ticks([v_min, v_max])

ax.set_xlim(0.5, 6)
ax.tick_params(labelsize=12)
fig.savefig(f"{imgPath}/fmt_vs_lambda_{observable}.png", transparent=True)

# %%
indx_partition = partition(df.density)[1]
df["mean_density"] = np.full(df.shape[0], 0)

for p in indx_partition:
    df["mean_density"][p] = np.full(p.size, df.density.values[p].mean())

# %%
cmap = mpl.cm.get_cmap("plasma")
normalized = mpl.colors.Normalize(
    vmin=df.mean_density.min(), vmax=df.mean_density.max()
)
fig, ax = plt.subplots(figsize=(6, 6))

for name, group in df.groupby("mean_density"):
    if name > 0.3:
        continue

    _, indices = partition(group.selection.values, threshold=0.002)
    selection = np.array([group["selection"].values[p].mean() for p in indices])
    v_fraction = np.array([group["v_fraction"].values[p].mean() for p in indices])
    v_fraction_std = np.array(
        [df["v_fraction"].values[p].std() / np.sqrt(p.size) for p in indices]
    )

    cs = cmap(normalized(group.mean_density.values[0]))
    ax.plot(selection, v_fraction, color=cs)
    ax.errorbar(selection, v_fraction, yerr=v_fraction_std, fmt="o", color=cs)

    # .v -> inf limit of z
    r = group.radius.values[0]
    nu = group.intensity.values[0]
    l = htsptSep(r * np.sqrt(3), name)
    S_crit = r / (np.sqrt(2) * l) * 2
    # ax.axvline(x=S_crit, color="black", linestyle="--")

v_min, v_max = df.mean_density.min(), df.mean_density.max()
norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
divider = make_axes_locatable(ax)
colorbar_axes = divider.append_axes("right", size="3%", pad=0.05)
cb = plt.colorbar(sm, cax=colorbar_axes, orientation="vertical")
cb.set_label(label="$\phi$", weight="bold")
cb.set_ticks(np.linspace(v_min, v_max, 5))

ax.set_xlabel(labels["selection"], fontsize=18)
ax.set_ylabel(labels["v_fraction"], fontsize=18)
ax.tick_params(labelsize=12)

# %%
