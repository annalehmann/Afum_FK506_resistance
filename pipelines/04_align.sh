for f in *R1*.fastq.gz
do
n="${f%_R1*}"
minimap2 \
-ax sr \
CEA10_genome.fa \
"$n"_R1_001.fastq.gz \
"$n"_R2_001.fastq.gz \|
samtools view -q 30 \|
samtools sort -o "$n"_filteredq30.bam
samtools index "$n"_filteredq30.bam
done
