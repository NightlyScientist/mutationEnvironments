# %%
import pandas as pd
import pyarrow.feather as feather
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os
import common.dataAPI as dataAPI
import common.datatables as dataAPIExtensions
import importlib
from collections import namedtuple

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

# %%
sns.set_theme(style="white")
np.seterr(divide="ignore", invalid="ignore")

labels = {
    "selection": r"$Selection$",
    "intensity": r"$Strength$",
    "mutation": r"$Mutation Rate$",
    "density": r"$Density$",
    "radius": r"$Radius$",
    "sectorPos": r"$\Delta x$",
    "sectorVar": r"$\sigma(y)$",
    "roughness": r"$\sigma$",
    "time": r"$y$",
}

def findTerminalxPos(x, L, threshold=10.0):
    terminal_index = np.where((np.abs(x) > (L - threshold)))
    if len(terminal_index) == 0:
        return len(x)
    return terminal_index[0][0]


def sectorAngle(df, x, y, opts):
    _x, _y = df[x].values, df[y].values
    terminal_point = findTerminalxPos(_x, opts["width"] // 2)
    theta = np.arctan2(_y[terminal_point], _x[terminal_point])
    return 180 / np.pi * theta


def getValues(df, x, y, skip=0):
    _x, _y = df[x].values[skip::], df[y].values[skip::]
    return _x, _y


def extractFit(x, y, scale=lambda t: t, _mask=None):
    X, Y = scale(x), scale(y)
    if _mask is None:
        _mask = np.fill(True, len(X))
    mask = np.where(np.isfinite(X) & np.isfinite(Y) & _mask)
    return stats.linregress(X[mask], Y[mask])


def drawLinearFit(ax, res, x):
    Y = res.intercept + res.slope * x
    label = f"$m = {res.slope:.2f} +/- {np.sqrt(res.stderr):.2f}$"
    ax.plot(x, Y, "r", label=label)


def drawPowerLawFit(ax, res, x):
    Y = np.exp(res.intercept) * x**res.slope
    label = f"$m = {res.slope:.2f} +/- {np.sqrt(res.stderr):.2f}$"
    ax.plot(x, Y, "r", label=label)


def retrieveData(path, skip=0):
    opts = pd.read_csv(os.path.join(path, "inputOpts.csv"), sep="\t", header=0).to_dict(
        orient="index"
    )[0]

    df = feather.read_feather(os.path.join(path, "boundary.arrow"))
    df["time"] = np.arange(0, len(df))
    df["roughness"] = np.sqrt(df["sectorVar"])
    df["sectorPos"] = df["sectorPos"] - df["sectorPos"].values[0]

    # _df = feather.read_feather(os.path.join(path, "mutationalFreq.arrow"))
    # df["mutationalFreq"] = _df["mutationalFreq"]
    return df, opts

def addMetricsColumns(df: pd.DataFrame):
    df[["xi_m", "xi_var", "v_fraction", "n_fraction", "time_extinction"]] = np.nan

    error_counter = 0 
    for index, path in enumerate(df.path):
        abs_path = os.path.join(path, "table.csv")
        if os.stat(abs_path).st_size == 0:
            error_counter += 1
            continue
        data = pd.read_csv(abs_path, header=0)


        # if there is a non-zero mutation rate, get background noise as mu * N
        width = df["width"][index]
        height = df["height"][index]
        mutation_probability = df["mutation"][index]
        mN = mutation_probability * width * height
        mLx = mutation_probability * width


        df.at[index, "time_extinction"] = (data.time_extinction / data.time)[
            data.time_extinction > 0
        ].mean()
        df.at[index, "xi_m"] = data.xi_m.mean()
        df.at[index, "xi_var"] = np.sqrt(data.xi_var.mean())
        df.at[index, "n_fraction"] = n_fraction
        df.at[index, "v_fraction"] = v_fraction
    print(f"Error counter: {error_counter}")

#%%

def phaseSpaceTable(basePath):
    dictList = []
    faulty = []
    for root, dirs, _ in os.walk(basePath):
        for directory in dirs:
            if not "env" in directory:
                continue
            subpath = os.path.join(root, directory)

            df, opts = retrieveData(subpath)

            # . fit the sector position
            data_pairs = [
                ("sectorPos", "sector_slope", lambda x: x),
                ("sectorVar", "roughness_slope", lambda x: np.log(x)),
            ]
            for k, n, f in data_pairs:
                x, y = getValues(df, "time", k, skip=10)
                t_p = findTerminalxPos(x, opts["width"] // 2, threshold=200)
                mask = (x > 10) & (x < t_p)
                fit = extractFit(x, y, scale=f, _mask=mask)
                opts[n] = fit.slope
                opts[n + "_stderr"] = np.sqrt(fit.stderr)

            opts["sector_angle"] = sectorAngle(df, "time", "sectorPos", opts)

            opts["path"] = subpath

            dictList.append(opts)
    return pd.DataFrame.from_dict(dictList), faulty


# %%
df, faulty = phaseSpaceTable(basePath)

# %%
# > phase diagram of the sector angle and roughness
from matplotlib import colors

x_axis = "mutation"
y_axis = "compensation"
z_axis = "sector_slope"
ps = df.pivot_table(index=y_axis, columns=x_axis, values=z_axis, aggfunc="mean")

vmin = np.min(df[z_axis].values)
vmax = np.max(df[z_axis].values)
# custom_norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

fig, ax = plt.subplots(figsize=(8, 6))
ax = sns.heatmap(
    ps,
    fmt=".2f",
    annot=False,
    ax=ax,
    cbar_kws={"label": ""},
    cmap="coolwarm",
    # norm=custom_norm,
)

radius = np.mean(df.radius) * 1000
density = np.mean(df.density)
num = 1 + 6 / 0.5  # range / step
intensity = num * np.linspace(np.min(df.intensity), np.max(df.intensity), 200)

cbar = ax.collections[0].colorbar
cbar.ax.set_ylabel("", fontsize=20)

ax.invert_yaxis()
ax.tick_params(labelsize=12)
ax.set_xlabel("")
ax.set_ylabel("")
fig.savefig(f"{imgPath}/phase_diagram_{x_axis}_{y_axis}_{z_axis}.png", transparent=True)

# %%
# > sample plot of the sector angle and roughness
import matplotlib.gridspec as gridspec

gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[:, 2:])
ax2 = plt.subplot(gs[0:2, 0:2])
ax3 = plt.subplot(gs[2::, 0:2])

if len(faulty) > 0:
    samplePath = faulty[0]
else:
    samplePath = np.random.choice(df.path.values)
    # samplePath = np.random.choice(df[(df.compensation==0) & (df.mutation == 0.0)].path.values)

image_path = os.path.join(samplePath, "snapshots_ID_3.png")
image = plt.imread(image_path)
ax1.imshow(image)
ax1.axis("off")

sample, opts = retrieveData(samplePath)
axes = (ax2, ax3)

axes[0].set(xscale="linear", yscale="linear")
axes[1].set(xscale="log", yscale="log")

# x = sample.time.values
# y = sample.sectorPos.values
x, y = getValues(sample, "time", "sectorPos", skip=10)
terminal_point = findTerminalxPos(x, opts["width"] // 2, threshold=25)
mask = (x > 10) & (x < terminal_point)
fit = extractFit(x, y, _mask=mask)

drawLinearFit(axes[0], fit, sample.time.values[0:terminal_point])
sns.lineplot(x=x[mask], y=y[mask], alpha=0.5, ax=axes[0])

f = lambda t: np.log(t)
# y = sample.sectorVar.values
x, y = getValues(sample, "time", "sectorVar", skip=10)
mask = (x > 10) & (x < terminal_point)
fit = extractFit(x, y, _mask=mask, scale=f)

drawPowerLawFit(axes[1], fit, sample.time.values[0:terminal_point])
sns.lineplot(x=x[mask], y=y[mask], alpha=0.5, ax=axes[1])

axes[0].set_xlabel(labels["time"], fontsize=14)
axes[0].set_ylabel(labels["sectorPos"], fontsize=14)

axes[1].set_xlabel(labels["time"], fontsize=14)
axes[1].set_ylabel(labels["sectorVar"], fontsize=14)
axes[1].legend(loc="upper left")

plt.tight_layout(h_pad=0.05, w_pad=0.05)
plt.show()

# %%
# > sample plot of the mutational frequency
ax = sns.lineplot(x=sample.time, y=sample.mutationalFreq)
