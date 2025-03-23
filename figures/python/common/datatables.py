import numpy as np
import os
import pandas as pd


def speciesFraction(df, num: int = 2, mN: float = 0, mLx: float = 0):
    n = [f"n_{i}" for i in range(1, num + 1)]
    v = [f"v_{i}" for i in range(1, num + 1)]

    n_fraction = ((df[n[-1]] - mLx) / (df[n].sum(axis=1))).mean()
    v_fraction = ((df[v[-1]] - mN) / (df[v].sum(axis=1))).mean()
    n_fraction = n_fraction if not np.isnan(n_fraction) else 0
    return n_fraction, v_fraction


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

        n_fraction, v_fraction = speciesFraction(data, num=2, mN=mN, mLx=mLx)

        df.at[index, "time_extinction"] = (data.time_extinction / data.time)[
            data.time_extinction > 0
        ].mean()
        df.at[index, "xi_m"] = data.xi_m.mean()
        df.at[index, "xi_var"] = np.sqrt(data.xi_var.mean())
        df.at[index, "n_fraction"] = n_fraction
        df.at[index, "v_fraction"] = v_fraction
    print(f"Error counter: {error_counter}")
