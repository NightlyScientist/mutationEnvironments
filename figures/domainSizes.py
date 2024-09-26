# %%
basePath = input()

# %%
import pandas as pd
import pyarrow.feather as feather
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os
import warnings
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d

#sns.set_theme(style="white")
imgPath = os.path.join(basePath, "images")
if not os.path.exists(imgPath):
    os.makedirs(imgPath)

# %%
def phaseSpaceTable(basePath):
    dictList = []
    for root, dirs, _ in os.walk(basePath):
        for dir in dirs:
            if not "env" in dir:
                continue
            subpath = os.path.join(root, dir)

            opts = pd.read_csv(
                os.path.join(subpath, "inputOpts.csv"), sep="\t", header=0
            ).to_dict(orient="index")[0]
            opts["path"] = subpath
            dictList.append(opts)
    return pd.DataFrame.from_dict(dictList)


def crossover_freq(N=10**6, alpha=2 / 5, beta=4):
    return N ** ((-1 + alpha) / (beta - alpha))


def crossover_P(N=10**6, alpha=2 / 5, beta=4):
    return N ** ((-alpha * beta + alpha) / (beta - alpha))


def pdf(inputPath):
    df = feather.read_feather(f"{inputPath}/domain_sizes.arrow")
    binSize = 1 / df.bin_counts.size
    df["sizes"] = [(i - 0.5) * binSize for i in range(1, 1 + df.bin_counts.size)]
    return df


def get_pdf(gdf, cdf=False):
    curveSets = dict()
    for key in gdf.groups.keys():
        group_df = gdf.get_group(key)
        ggdf = group_df.groupby(["selection", "mutation"])

        curveSet = dict()
        for i, second_key in enumerate(ggdf.groups):
            sizes = dict()
            bin_counts = dict()

            sims = ggdf.get_group(second_key)
            for j, path in enumerate(sims.path):
                sim = pdf(path)
                sizes[j] = sim.sizes.values
                bin_counts[j] = sim.bin_counts.values

            x = np.vstack(list(sizes.values()))
            y = np.vstack(list(bin_counts.values()))
            bin_counts = np.sum(y, axis=0)

            x = np.mean(x[:, bin_counts > 1], axis=0)
            y = np.sum(y[:, bin_counts > 1], axis=0)

            curveSet[second_key] = (x, y / np.sum(y))
        curveSets[key] = curveSet
    return curveSets


# %%
# .group entries
df = phaseSpaceTable(basePath)
gdf = df.groupby("intensity")
curveSets = get_pdf(gdf, cdf=True)

# %%
# >cdfs for various parameters
cmap = mpl.colormaps.get_cmap("plasma")

counter = 0
for intensity, curveSet in curveSets.items():
    counter += 1
    #if intensity < 6.0:
        #continue

    fig, ax = plt.subplots(nrows=1, figsize=(8, 8))
    ax.set(xscale="log", yscale="log")
    ax.set_ylim(10 ** (-6), 1)

    for (s, m), curve in curveSet.items():
        x = curve[0]
        y = 1 - np.cumsum(curve[1])
        ax.scatter(x, y, alpha=0.5, color=cmap(float(s) / 0.10), s=20)

    # .add scaling for bubbles
    N = 1000 * 1500
    px = np.linspace(x.min(), crossover_freq(N=N), 1000)
    py = (N * px) ** (-2 / 5)
    ax.plot(px, py, color="black", linestyle="--")
    ax.text(px[10] * 2, py[10] * 1.1, "$x^{-2/5}$", fontsize=22)

    # .add scaling for sectors
    if counter == 1:
        px = np.linspace(crossover_freq(N=N), 1, 1000)
        py = (px ** (-5) / N) **2
        ax.plot(px, py, color="black", linestyle="--")
        ax.text(0.2, 0.0001, "$x^{-10}$", fontsize=22)

    #. add crossover frequency
    ax.text(crossover_freq(N=N) + 0.01, crossover_P(N=N), "$x_c$", fontsize=22)

    # .add line with slope x^{-1}
    # px = np.linspace(x.min(), 1, 1000)
    # py = (10**6 * px) ** (-1)
    # ax.plot(px, py, color="black", linestyle="--")

    # .add crossover length
    if counter > 1:
        ax.axvline(x=crossover_freq(N=N, beta=4), color="black", linestyle="--")

    plt.tick_params(axis="both", which="major", labelsize=16)
    plt.tick_params(axis="both", which="minor", labelsize=16)

    # .add colorbar
    v_min, v_max = [0, 0.1]
    norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    divider = make_axes_locatable(ax)
    colorbar_axes = divider.append_axes("right", size="3%", pad=0.05)
    cb = plt.colorbar(sm, cax=colorbar_axes, orientation="vertical")
    cb.set_ticks([v_min, v_max])

    fig.savefig(f"{imgPath}/clone_sizes_{intensity}.png", transparent=True)
    plt.show()
    plt.close(fig)
    #break

#%%
#> compare neutral evolution in both landscape scenarios
fig, ax = plt.subplots(nrows=1, figsize=(8, 8))
ax.set(xscale="log", yscale="log")
ax.set_ylim(10 ** (-6), 1)

counter = 0 
colors = ["black", "green"]
legend = ["Uniform Environment", "Disordered Hotspot Landscape"]

for intensity, curveSet in curveSets.items():
    if 0 < intensity < 6.0:
        continue

    print(intensity)

    for (s, m), curve in curveSet.items():
        x = curve[0]
        y = 1 - np.cumsum(curve[1])
        ax.scatter(x, y, alpha=0.5, c=colors[counter], s=20, label=legend[counter])
        break
    counter += 1

    # .add scaling for bubbles
    N = 1000 * 1500
    px = np.linspace(x.min(), crossover_freq(N=N), 1000)
    py = (N * px) ** (-2 / 5)
    ax.plot(px, py, color="black", linestyle="--")
    ax.text(px[10] * 2, py[10] * 1.1, "$x^{-2/5}$", fontsize=22)

    # .add scaling for sectors
    px = np.linspace(crossover_freq(N=N), 1, 1000)
    py = (px ** (-5) / N) **2
    ax.plot(px, py, color="black", linestyle="--")
    ax.text(0.3, 0.0001, "$x^{-10}$", fontsize=22)

    #. add crossover frequency
    ax.text(crossover_freq(N=N) + 0.01, crossover_P(N=N), "$x_c$", fontsize=22)

    # .add line with slope x^{-1}
    # px = np.linspace(x.min(), 1, 1000)
    # py = (10**6 * px) ** (-1)
    # ax.plot(px, py, color="black", linestyle="--")

    # .add crossover length
    #ax.axvline(x=crossover_freq(N=N, beta=4), color="black", linestyle="--")

plt.tick_params(axis="both", which="major", labelsize=16)
plt.tick_params(axis="both", which="minor", labelsize=16)
plt.legend(fontsize=13)
fig.savefig(f"{imgPath}/clone_sizes_comparison.png", transparent=True)
plt.show()


# %%
fig, ax = plt.subplots(ncols=1, figsize=(6, 6))
# ax.set(xscale="log", yscale="log")

interpolated_curves = dict()
x_range = np.linspace(10 ** (-4), 10 ** (-1), 1000)

for intensity, curveSet in curveSets.items():
    interpolated_curves_set = dict()
    for (s, m, nu), curve in curveSet.items():
        x = curve[0]
        y = 1 - curve[1]

        # Perform linear spline interpolation
        f = interp1d(x, y, kind="linear")
        x_interp = np.linspace(x.min(), x.max(), 1000)
        y_interp = f(x_interp)

        interpolated_curves_set[(s, m, nu)] = (f, x.min(), x.max())
    interpolated_curves[intensity] = interpolated_curves_set

data = dict()
for intensity, curveSet in interpolated_curves.items():
    ref_key = list(curveSet.keys())[0]
    ref_f, ref_x_min, ref_x_max = curveSet[ref_key]
    for (s, m, nu), (f, x_min, x_max) in curveSet.items():
        if ref_key == (s, m, nu):
            continue
        _mask = (x_range >= x_min) & (x_range <= x_max)
        y_interp = f(x_range[_mask])
        y_ref = ref_f(x_range[_mask])
        result = np.sum((y_interp - y_ref) ** 2)

        if s not in data:
            data[s] = [(intensity, result)]
        data[s].append((intensity, result))

for s in data.keys():
    intensity, result = zip(*data[s])
    ax.plot(intensity, result, color=cmap(float(s) / 0.10))

# ax.errorbar(
#     xaxis,
#     np.array(avg_res),
#     yerr=np.std(np.array(avg_res)) / np.sqrt(len(avg_res)),
#     fmt="o",
#     color="black",
# )
# ax.set_ylabel(r"$\langle \sum_j (y_i(x_j) - y_0(x_j)^2) \rangle_i$")
ax.set_ylabel("$\sum \Delta(P(X > x))$", fontsize=20)
ax.set_xlabel("Intensity", fontsize=20)
ax.tick_params(axis="both", which="major", labelsize=20)
ax.tick_params(axis="both", which="minor", labelsize=20)

#%%
def htsptSep(r: int, p: float | list | np.ndarray):
    return np.sqrt(-np.pi * r**2 / np.log(1 - p))
htsptSep(10 * 1.15, 0.1)

X_c = 0.1
A = X_c * df.width.values[0] * df.height.values[0]
lperp = A ** (2/5)
lpar = lperp ** (3/2)
print(lpar)

# %%
