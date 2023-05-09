grep '	start_codon	.*	+	.*transcript_support_level "1"' ../../data/gencode.vM27.annotation.gtf | \
 sed -e 's/	[^	]*transcript_id "\([^"]*\)".*$/	\1/g' > gencode-start.gtf

head gencode-start.gtf

# Because of autoindentation of sublime text, this command has been edited and excuted manually
#grep ' exon .*  + ' ../../data/gencode.vM27.annotation.gtf | \
# sed -e 's/ [^  ]*transcript_id "\([^"]*\)".*$/ \1/g' > gencode-plusexon.gtf

head gencode-plusexon.gtf

# Because of autoindentation of sublime text, this command has been edited and excuted manually
#bedtools intersect -a gencode-start.gtf -b gencode-plusexon.gtf -wa -wb | \
# awk -F'    ' -v OFS='  ' '$9 == $18 { print $10, $13-1, $14, $18, $4-1, $16; }' | \
# sort -k1,1 -k2,3n -k4,4 > gencode-exons-containing-startcodon.bed

head gencode-exons-containing-startcodon.bed; tail gencode-exons-containing-startcodon.bed

(samtools view -H ../../data/project/RPF-siLuc.bam; \
  samtools view -F20 ../../data/project/RPF-siLuc.bam | \
  bioawk -c sam '{ if (length($seq) >= 25) print $0; }') | \
 samtools view -b -o filtered-RPF-siLuc.bam

ls -alh *RPF-siLuc.bam

bedtools genomecov -ibam filtered-RPF-siLuc.bam -bg -5 > fivepcounts-RPF-siLuc.bed

head fivepcounts-RPF-siLuc.bed

bedtools intersect -a fivepcounts-RPF-siLuc.bed -b gencode-exons-containing-startcodon.bed -wa -wb -nonamecheck > fivepcounts-filtered-RPF-siLuc.txt

head fivepcounts-filtered-RPF-siLuc.txt