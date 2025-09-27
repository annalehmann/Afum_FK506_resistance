mkdir filteredVariants

# get list of snps (all) in progenitor A1163 strain
grep -o 'AFUB_[[:alnum:]]\{6\}' snpEff_output/A1163-YPD_S79_L008_snpEff_all.vcf | \
uniq > filteredVariants/A1163_parent_all_snps_geneIDs.txt

# get lists of all and unique gene IDs in each vcf file
for f in snpEff_output/*all.vcf
do
n=$(basename "$f" _snpEff_all.vcf)
    grep -o 'AFUB_[[:alnum:]]\{6\}' "$f" | \
    uniq > filteredVariants/"$n"_all_snps_geneIDs.txt
    diff filteredVariants/A1163_parent_all_snps_geneIDs.txt filteredVariants/"$n"_all_snps_geneIDs.txt | \
    grep '> AFUB_[[:alnum:]]\{6\}' | \
    cut -c 3-13 > filteredVariants/"$n"_all_unique_snps.txt
done

# get list of snps (coding only) in progenitor A1163 strain
grep -o 'AFUB_[[:alnum:]]\{6\}' snpEff_output/A1163-YPD_S79_L008_snpEff_coding.vcf | \
uniq > filteredVariants/A1163_parent_coding_snps_geneIDs.txt

# get lists of all and unique gene IDs in each vcf file
for f in snpEff_output/*coding.vcf
do
n=$(basename "$f" _snpEff_coding.vcf)
    grep -o 'AFUB_[[:alnum:]]\{6\}' "$f" | \
    uniq > filteredVariants/"$n"_coding_snps_geneIDs.txt
    diff filteredVariants/A1163_parent_coding_snps_geneIDs.txt filteredVariants/"$n"_coding_snps_geneIDs.txt | \
    grep '> AFUB_[[:alnum:]]\{6\}' | \
    cut -c 3-13 > filteredVariants/"$n"_coding_unique_snps.txt
done

cd filteredVariants

for f in ../snpEff_output/*all.vcf
do
n=$(basename "$f" _snpEff_all.vcf)
    while read -r line; do
        gene_ID="$line"
        awk -v gene_ID="$gene_ID" '$0 ~ gene_ID' "$f"
    done < "$n"_all_unique_snps.txt > "$n"_all_unique_snps.vcf
done

for f in ../snpEff_output/*coding.vcf
do
n=$(basename "$f" _snpEff_coding.vcf)
    while read -r line; do
        gene_ID="$line"
        awk -v gene_ID="$gene_ID" '$0 ~ gene_ID' "$f"
    done < "$n"_coding_unique_snps.txt > "$n"_coding_unique_snps.vcf
done

# get number of snps in parent
grep -o 'AFUB_[[:alnum:]]\{6\}' snpEff_output/A1163-YPD_S79_L008_snpEff_all.vcf | \
sort -u > A1163_parent_all_snps_geneIDs.txt

grep -o 'AFUB_[[:alnum:]]\{6\}' snpEff_output/A1163-YPD_S79_L008_snpEff_coding.vcf | \
sort -u > A1163_parent_coding_snps_geneIDs.txt
