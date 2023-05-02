import pandas as pd

cnts = pd.read_csv("read-counts.txt", sep="\t", comment="#", index_col=0)


# Rename columns (read-counts file generated from relative [path])
column_renamer = {}
for col in cnts.columns:
    if col.startswith("../../"):
        column_renamer[col] = col.split("/")[-1]

print(column_renamer)

cnts = cnts.rename(columns=column_renamer)
print(cnts.columns)
print(cnts)

"""
cnts["clip_enrichment"] = cnts["CLIP-35L33G.bam"] / cnts["RNA-control.bam"]
cnts["rden_change"] = (cnts["RPF-siLin28a.bam"] / cnts["RNA-siLin28a.bam"]) / (
    cnts["RPF-siLuc.bam"] / cnts["RNA-siLuc.bam"]
)
"""
# Remove rows if any of the values are 0
cnts["Check_0"] = (
    cnts["CLIP-35L33G.bam"]
    * cnts["RNA-control.bam"]
    * cnts["RPF-siLin28a.bam"]
    * cnts["RNA-siLin28a.bam"]
    * cnts["RPF-siLuc.bam"]
    * cnts["RNA-siLuc.bam"]
)

# Read ount cutoff
cnts = cnts[cnts["CLIP-35L33G.bam"] >= 30]
cnts = cnts[cnts["RNA-control.bam"] >= 30]
cnts = cnts[cnts["RNA-siLin28a.bam"] >= 30]
cnts = cnts[cnts["RNA-siLuc.bam"] >= 30]
cnts = cnts[cnts["RPF-siLuc.bam"] >= 80]


# Get RPKM
cnts["CLIP-35L33G.bam"] = (
    cnts["CLIP-35L33G.bam"].astype(float)
    / cnts["Length"].astype(float)
    * 1000
    / sum(cnts["CLIP-35L33G.bam"])
    * 10**6
)
cnts["RNA-control.bam"] = (
    cnts["RNA-control.bam"].astype(float)
    / cnts["Length"].astype(float)
    * 1000
    / sum(cnts["RNA-control.bam"])
    * 10**6
)
cnts["RPF-siLin28a.bam"] = (
    cnts["RPF-siLin28a.bam"].astype(float)
    / cnts["Length"].astype(float)
    * 1000
    / sum(cnts["RPF-siLin28a.bam"])
    * 10**6
)
cnts["RNA-siLin28a.bam"] = (
    cnts["RNA-siLin28a.bam"].astype(float)
    / cnts["Length"].astype(float)
    * 1000
    / sum(cnts["RNA-siLin28a.bam"])
    * 10**6
)
cnts["RPF-siLuc.bam"] = (
    cnts["RPF-siLuc.bam"].astype(float)
    / cnts["Length"].astype(float)
    * 1000
    / sum(cnts["RPF-siLuc.bam"])
    * 10**6
)
cnts["RNA-siLuc.bam"] = (
    cnts["RNA-siLuc.bam"].astype(float)
    / cnts["Length"].astype(float)
    * 1000
    / sum(cnts["RNA-siLuc.bam"])
    * 10**6
)


cnts["clip_enrichment"] = cnts["CLIP-35L33G.bam"] / cnts["RNA-control.bam"]

cnts["rden_change"] = (cnts["RPF-siLin28a.bam"] / cnts["RNA-siLin28a.bam"]) / (
    cnts["RPF-siLuc.bam"] / cnts["RNA-siLuc.bam"]
)


from matplotlib import pyplot as plt
import numpy as np

cnts["x"] = np.log2(cnts["clip_enrichment"])
cnts["y"] = np.log2(cnts["rden_change"])

cnts.to_csv("Debug.csv")


print(cnts)


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.scatter(
    np.log2(cnts["clip_enrichment"]),
    np.log2(cnts["rden_change"]),
    s=1,
    c="#000000",
    alpha=0.3,
)

# plt.savefig("week1_1.svg")

import ssl

ssl._create_default_https_context = ssl._create_unverified_context
mouselocal = pd.read_csv(
    "https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt",
    sep="\t",
    index_col=0,
)

print(mouselocal.head())

# mouselocal.to_csv("local.csv")

cnts.index = cnts.index.str.split(".").str[0]
cnts = pd.merge(cnts, mouselocal, left_index=True, right_index=True)


colors = {
    "nucleus": "tab:blue",
    "cytoplasm": "tab:green",
    "integral membrane": "tab:pink",
}


cnts["color"] = cnts["type"].apply(lambda x: colors[x])

plt.cla()
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
plt.xlim([-6, 4])
plt.ylim([-2, 3])
ax.scatter(
    np.log2(cnts["clip_enrichment"]),
    np.log2(cnts["rden_change"]),
    s=1.5,
    c=cnts["color"],
    cmap="jet",
)

plt.savefig("week1_2.svg")
