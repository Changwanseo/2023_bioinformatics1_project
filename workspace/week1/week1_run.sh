# Check md5sum
md5sum ../../data/*

# Feature counts
featureCounts -a ../../data/gencode.vM27.annotation.gtf -o read-counts.txt ../../data/project/*.bam