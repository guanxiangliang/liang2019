for i in D*
do
hmmscan --domtblout ${i}/hmmout --domE 1e-5 /media/lorax/users/guanxiang/2_veo/liang.veo/input/database/Pfam-A.hmm ${i}/${i}.contig.3k.orf.fa
awk 'BEGIN{print "target""\t""acc""\t""tlen""\t""query""\t""qlen""\t""evalue"}{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' ${i}/hmmout | sed '/\#/d' > ${i}/${i}.hmm.out
done


for i in D*
do
grep ">" ${i}/${i}.contig.3k.orf.fa | awk -F\# '{print $1"\t"$1"\t"$2"\t"$3}' | sed -e 's/\ //g' | sed -e 's/>//g' |awk '{gsub("_\\S+","",$2); print}' | awk '$5="+"' | awk 'BEGIN{print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${i}/${i}.orf.saf
featureCounts -a ${i}/${i}.orf.saf -F SAF -o ${i}/${i}.orf.count ../contigs/${i}/${i}.bam -T 1
done
