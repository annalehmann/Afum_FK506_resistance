for f in *R1*.fastq.gz
do
n="${f%_R1*}"
snippy --reference Aspergillus_fumigatusa1163.ASM15014v1.dna.nonchromosomal.fa \
--R1 "$n"_R1_001.fastq.gz \
--R2 "$n"_R2_001.fastq.gz \
--outdir "$n"_snippy_results
done
