import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


df = pd.read_csv("fivepcounts-filtered-RPF-siLuc.txt", sep="\t", header=None)


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

# print(list(df["pos"]))

plt.figure(figsize=(10, 2))
plt.ylabel("siLuc - Raw read count")
plt.xticks(np.arrange(-50, 50, 10))
plt.hist(hist_list, bins=range(-50, 50), color="black")

plt.savefig("Histogram.svg")
