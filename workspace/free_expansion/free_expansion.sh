# Should consider more about which sites to grep



# Grep from gencode.gtf. Expansion this to several informative sites with iteration and paraellization
grep -i mirlet7g gencode.gtf

# Loop multiple bam
samtools view -b -o CLIP-let7g.bam CLIP-35L33G.bam chr9:106056039-106056126
samtools view CLIP-let7g.bam | wc -l
samtools mpileup CLIP-let7g.bam > CLIP-let7g.pileup
wc -l CLIP-let7g.pileup

awk '$2 >= 106056039 && $2 <= 106056126 { print $0; }' CLIP-let7g.pileup > CLIP-let7g-gene.pileup
tail CLIP-let7g-gene.pileup

# Pandas works, make each of the pileup files and merge


# Genomewide visualization
