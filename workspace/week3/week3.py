import pandas as pd
import re

pileup = pd.read_csv(
    "CLIP-let7g-gene.pileup",
    sep="\t",
    names=["chrom", "pos", "_ref", "count", "basereads", "quals"],
)
pileup.tail()

toremove = re.compile("[<>$*#^]")
pileup["matches"] = pileup["basereads"].apply(lambda x: toremove.sub("", x))

print(pileup[["chrom", "pos", "matches"]])
print(pileup[pileup["pos"] == 106056094].iloc[0]["matches"])

print(pileup)
