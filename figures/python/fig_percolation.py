# %%
import importlib
import pandas as pd
import pyarrow.feather as feather
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import common.dataAPI as dataAPI
import common.datatables as dataAPIExtensions

# reload modules (helpful for debugging)
importlib.reload(dataAPI)

workspace_paths = dataAPI.fetchWorkspaceEnv("../../")
print(workspace_paths.keys())
print(workspace_paths["top_level_path"])

ensemble_paths = dataAPI.fetchEnsemblePaths(workspace_paths["top_level_path"])

for i, path in enumerate(ensemble_paths):
    print(f"[{i}] ", path.base)

# %%
# pick one of the ensembles using 'example_ensemble_path'
selector = 0
#example_ensemble_path = ensemble_paths[selector]
_path = "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/experiments/percolation/percolation_limits/"
example_ensemble_path = dataAPI.setOutputPaths(_path)
print(example_ensemble_path.base)

# ensemble table generated from ensemble directory and options.csv
ensemble_table = dataAPI.ensembleTableInfo(example_ensemble_path.base)

importlib.reload(dataAPIExtensions)
dataAPIExtensions.addMetricsColumns(ensemble_table)


# %%
def partition(x, threshold=0.002):
    _x = np.sort(np.array(x))
    _indx = np.argsort(np.array(x))
    cuts = np.where(np.diff(_x) > threshold)
    x_partition = np.split(_x, cuts[0] + 1)
    indx_partition = np.split(_indx, cuts[0] + 1)
    return x_partition, indx_partition


def htsptSep(r: int, p: float | list | np.ndarray):
    return r * np.sqrt(-np.pi / np.log(1 - p))


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

# %%
df = ensemble_table
limiter = 0.09
observable = "n_fraction"
cmap = mpl.colormaps.get_cmap("plasma")

fig, ax = plt.subplots(ncols=1, figsize=(6, 6))

counter = 0
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
    obs_std = np.array([
        group[observable].values[p].std() / np.sqrt(p.size) for p in indx_partition
    ])

    ax.fill_between(
        np.sqrt(2) * hl / diameter,
        obs - obs_std,
        obs + obs_std,
        color=cmap(S / 0.10),
        alpha=0.3,
    )
    ax.errorbar(
        np.sqrt(2) * hl / diameter,
        # density,
        obs,
        yerr=obs_std,
        fmt="o",
        color=cmap(S / 0.10),
    )
    #ax.plot(density / diameter, v_fraction, linestyle="-", color=cmap(S / 0.10))

    r = group.radius.values[0]
    # r = group.radius.values[0] * np.sqrt(3)
    # intensity = group.intensity.values[0]
    # s = group.selection.values[0]
    # z_c = 2 * r * intensity / (1 + intensity) / s

    # .percolation of overlapping circles
    for z in [htsptSep(r, 0.5), htsptSep(r, 0.68)]:
        ax.axvline(x=np.sqrt(2) * z / diameter, color="black", linestyle="--")
        break

    # add phi = 0.68 line
    # ax.axvline(x= 0.68, color="black", linestyle="--")
    # ax.axvline(x= 0.54, color="black", linestyle="--")
    # ax.axvline(x= 0.4, color="black", linestyle="--")
    #if counter == 2:
    #    break
    #counter += 1

v_min, v_max = [0, 0.1]
norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
divider = make_axes_locatable(ax)
colorbar_axes = divider.append_axes("right", size="3%", pad=0.05)
cb = plt.colorbar(sm, cax=colorbar_axes, orientation="vertical", label="")
cb.set_ticks([v_min, v_max])

ax.set_xlim(0.5, 6)
ax.tick_params(labelsize=12)
fig.savefig(f"{example_ensemble_path.img}/fmt_vs_lambda_{observable}.png", transparent=True)
plt.show()


#%%
indx_partition = partition(df.density)[1]
df["mean_density"] = np.full(df.shape[0], 0)

for p in indx_partition:
    #df["mean_density"][p] = np.full(p.size, df.density.values[p].mean())
    df.iloc[p, "mean_density"] = np.full(p.size, df.density.values[p].mean())

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
