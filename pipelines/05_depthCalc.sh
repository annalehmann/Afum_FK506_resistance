for f in *_filteredq30.bam
do
n="${f%_filteredq30.bam}"
mosdepth "$n" "$f" \
--by 500 \
--threads 2
done
