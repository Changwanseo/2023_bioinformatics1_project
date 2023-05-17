import pandas as pd
import re
import math


def write_bedgraph(gene, position):
    pileup = pd.read_csv(
        f"CLIP-{gene}-gene.pileup",
        sep="\t",
        names=["chrom", "pos", "_ref", "count", "basereads", "quals"],
    )
    pileup.tail()

    toremove = re.compile("[<>$*#^]")
    pileup["matches"] = pileup["basereads"].apply(lambda x: toremove.sub("", x))

    # print(pileup[["chrom", "pos", "matches"]])
    # print(pileup[pileup["pos"] == 106056094].iloc[0]["matches"])

    # Assisted with chatgpt
    def calculate_shannon_entropy(string):
        base_counts = {}
        total_bases = 0

        for base in string:
            if base in base_counts:
                base_counts[base] += 1
            else:
                base_counts[base] = 1
            total_bases += 1

        entropy = 0.0
        for count in base_counts.values():
            frequency = count / total_bases
            entropy -= frequency * math.log2(frequency)
        return entropy

    pileup["shannon"] = pileup["matches"].apply(calculate_shannon_entropy)

    pileup["end"] = pileup["pos"] + 1

    bedgraph = pileup[["chrom", "pos", "end", "shannon"]]
    bedgraph.to_csv("tmp_csv", sep="\t", header=False, index=False)

    header = f'browser position {position}\nbrowser hide all\nbrowser pack refGene encodeRegions\nbrowser full altGraph\ntrack type=bedGraph name="{gene}" description="{gene}" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n'

    with open("tmp_csv", "r") as f:
        content = f.read()
        with open(f"{gene}.bedgraph", "a") as fw:
            fw.write(header)
            fw.write(content)


write_bedgraph(gene="let7g", position="chr9:106056039-106056126")
write_bedgraph(gene="let7f-1", position="chr13:48691305-48691393")
write_bedgraph(gene="let7d", position="chr13:48689488-48689590")
