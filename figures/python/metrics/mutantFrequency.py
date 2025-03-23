from pyarrow import feather
import numpy as np
import os
import pandas as pd

def heatmap(path: str, lx, ly, trials):
    """fetch hatmap of mutant frequency"""
    f = feather.read_feather(f"{path}/heatmap_ID3.arrow").heatmap_ID3.values
    hm = (f.reshape((ly, lx)) / trials) - 1
    return np.flip(hm, axis=0)

def fetchMutantFrequency_y(path, skip=0, file="mutationalFreq.arrow") -> pd.DataFrame:
    """fetch x-averaged mutant frequency as a function of expansion distance"""
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