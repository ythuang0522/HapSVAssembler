curDir=`pwd`

echo -e "\n$1" >> summary
$curDir/contigLength.pl $1 length
sort -nr length > length_sort
$curDir/N50 length_sort 0
rm length_sort length
echo -n "contig num:" >> summary
grep '>' $1 | wc -l >> summary
echo -n "#N:" >> summary
grep -o 'N' $1 | wc -l >> summary
