import subprocess
import pandas as pd


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
