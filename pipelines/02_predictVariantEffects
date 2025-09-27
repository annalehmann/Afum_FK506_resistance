# import A1163 database
java -jar ~/snpEff/snpEff.jar download -v Aspergillus_fumigatusa1163

mkdir snpEff_output

for f in *_snippy_results
do
n="${f%_snippy_results}"
java -jar ~/snpEff/snpEff.jar ann \
-no-downstream \
-no-upstream \
-no-intergenic \
-c ~/snpEff/snpEff.config \
-o vcf \
Aspergillus_fumigatusa1163 \
"$f"/snps.vcf > snpEff_output/"$n"_snpEff_coding.vcf

java -jar ~/snpEff/snpEff.jar ann \
-c ~/snpEff/snpEff.config \
-o vcf \
Aspergillus_fumigatusa1163 \
"$f"/snps.vcf > snpEff_output/"$n"_snpEff_all.vcf
done
