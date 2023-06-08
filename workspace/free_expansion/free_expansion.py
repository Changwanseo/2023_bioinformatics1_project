import subprocess
import pandas as pd
import re
import math


## Grep all mirlet7 from gencode.gtf and save to mirlet.gtf
cmd = "grep -i mirlet ../../data/gencode.vM27.annotation.gtf > mirlet.gtf"
run = subprocess.call(cmd, shell=True)

## Read mirlet.gtf and process available regions
df = pd.read_csv(
    "mirlet.gtf",
    sep="\t",
    names=[
        "chr",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ],
    header=None,
)

## To remove duplicates, only get transcript rows
# Manually checked gene/transcript/exon showed same starting and end point
df = df[df["feature"] == "transcript"]


## Shorten attribute dataframe to show gene name
# A function to simplify attribute column
def get_genename(string):
    return string.split('gene_name "')[1].split('";')[0]


df["attribute"] = df["attribute"].apply(get_genename)
# Reset index of dataframe to enumerate
df.reset_index(inplace=True)

df.to_excel("Mirlet_list.xlsx")


###########################################################################################


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


## Loop week3 works
for n, gene in enumerate(df["attribute"]):
    cmds = [
        f"samtools view -b -o CLIP-{gene}.bam ../../data/project/CLIP-35L33G.bam {df['chr'][n]}:{df['start'][n]}-{df['end'][n]}",
        f"samtools view CLIP-{gene}.bam | wc -l",
        f"samtools mpileup CLIP-{gene}.bam > CLIP-{gene}.pileup",
        f"wc -l CLIP-{gene}.pileup",
        f"awk '$2 >= {df['start'][n]} && $2 <= {df['end'][n]} {{ print $0; }}' CLIP-{gene}.pileup > CLIP-{gene}-gene.pileup",
    ]

    for cmd in cmds:
        run = subprocess.call(cmd, shell=True)

    write_bedgraph(gene, f"{df['chr'][n]}:{df['start'][n]}-{df['end'][n]}")

    ## Make track files
    cmds = [
        f"make_tracks_file --trackFiles {gene}.bedgraph -o {gene}.ini",
        f"pyGenomeTracks --tracks {gene}.ini --region {df['chr'][n]}:{df['start'][n]}-{df['end'][n]} --trackLabelFraction 0.2 --width 38 --dpi 300 --fontSize 20 -o {gene}.svg",
    ]

    for cmd in cmds:
        run = subprocess.call(cmd, shell=True)

    """
    samfile = pysam.AlignmentFile("CLIP-Mirlet7g.bam", "rb")
    for read in 
    """
