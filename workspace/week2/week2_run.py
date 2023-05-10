import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def get_histogram_data(file):
    df = pd.read_csv(file, sep="\t", header=None)

    df.columns = [
        "5pchr",
        "5pstart",
        "5pend",
        "depth",
        "exonchr",
        "exonstart",
        "exonend",
        "gene",
        "start_codon",
        "strand",
    ]

    df["pos"] = df["5pstart"] - df["start_codon"]

    hist_list = []
    for n, i in enumerate(list(df["pos"])):
        hist_list += [i] * df["depth"][n]

    return hist_list


# print(list(df["pos"]))

plt.figure(figsize=(10, 4))
fig, ax = plt.subplots(2, 1)

ax[0].set_ylabel("siLuc - Raw read count")
ax[0].set_xticks(np.arange(-50, 50, 10))
ax[0].hist(
    get_histogram_data("fivepcounts-filtered-RPF-siLuc.txt"),
    bins=range(-50, 50),
    color="black",
)


ax[1].set_ylabel("siLin28a - Raw read count")
ax[1].set_xticks(np.arange(-50, 50, 10))
ax[1].hist(
    get_histogram_data("fivepcounts-filtered-RPF-siLin28a.txt"),
    bins=range(-50, 50),
    color="black",
)

plt.savefig("Histogram.svg")
