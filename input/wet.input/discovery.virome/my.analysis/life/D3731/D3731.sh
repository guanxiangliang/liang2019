source activate lifecyc
cd /home/guanxian/new/D3731
for i in *.viral.contig.txt
do
cat ${i} | while read p
do
base=$(echo ${p}| awk -F. '{print $1}')
awk "/^>/{x = /${p}_/;}(x)" /home/guanxian/new/orfs/${base}/${base}.contig.3k.orf.fa > ${p}.fa
perl ~/software/PHACTS/phacts.pl -f ${p}.fa --classes ~/software/PHACTS/classes_lifestyle -r 100 > ${p}.lifestyle
rm ${p}.fa
done
done
