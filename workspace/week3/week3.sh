grep -i mirlet7g ../../data/gencode.vM27.annotation.gtf

samtools view -b -o CLIP-let7g.bam ../../data/project/CLIP-35L33G.bam chr9:106056039-106056126
samtools view CLIP-let7g.bam | wc -l
samtools mpileup CLIP-let7g.bam > CLIP-let7g.pileup
wc -l CLIP-let7g.pileup
head CLIP-let7g.pileup
awk '$2 >= 106056039 && $2 <= 106056126 { print $0; }' CLIP-let7g.pileup > CLIP-let7g-gene.pileup
tail CLIP-let7g-gene.pileup

#chr13:48691305-48691393
grep -i mirlet7f-1 ../../data/gencode.vM27.annotation.gtf
samtools view -b -o CLIP-let7f-1.bam ../../data/project/CLIP-35L33G.bam chr13:48691305-48691393
samtools view CLIP-let7f-1.bam | wc -l
samtools mpileup CLIP-let7f-1.bam > CLIP-let7f-1.pileup
wc -l CLIP-let7f-1.pileup
head CLIP-let7f-1.pileup
awk '$2 >= 48691305 && $2 <= 48691393 { print $0; }' CLIP-let7f-1.pileup > CLIP-let7f-1-gene.pileup
tail CLIP-let7f-1-gene.pileup



#chr13:48689488-48689590
grep -i mirlet7d ../../data/gencode.vM27.annotation.gtf
samtools view -b -o CLIP-let7d.bam ../../data/project/CLIP-35L33G.bam chr13:48689488-48689590
samtools view CLIP-let7d.bam | wc -l
samtools mpileup CLIP-let7d.bam > CLIP-let7d.pileup
wc -l CLIP-let7d.pileup
head CLIP-let7d.pileup
awk '$2 >= 48689488 && $2 <= 48689590 { print $0; }' CLIP-let7d.pileup > CLIP-let7d-gene.pileup
tail CLIP-let7d-gene.pileup