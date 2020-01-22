for i in /media/lorax/users/guanxiang/1_igram/1_virome/analysis/3_vlp_vliadtion/sunbeam/sunbeam_output/qc/decontam/D*_1*
do
base=$(basename ${i} | awk -F_ '{print $1}')
mkdir ${base}
megahit -1 /media/lorax/users/guanxiang/1_igram/1_virome/analysis/3_vlp_vliadtion/sunbeam/sunbeam_output/qc/decontam/${base}_1.fastq.gz -2 /media/lorax/users/guanxiang/1_igram/1_virome/analysis/3_vlp_vliadtion/sunbeam/sunbeam_output/qc/decontam/${base}_2.fastq.gz -o ${base}/contig
done


for i in D*
do
vsearch --sortbylength ${i}/contig/final.contigs.fa  --relabel ${i}.contig   --outpu  ${i}/${i}.contig.fa
vsearch --sortbylength ${i}/${i}.contig.fa  --relabel ${i}.contig  --minseqlength 3000 --maxseqlength -1 --outpu  ${i}/${i}.contig.3k.fa
prodigal -i ${i}/${i}.contig.3k.fa  -a ${i}/${i}/${i}.contig.3k.orf.fa -p meta
hmmscan --domtblout ${i}/hmmout /media/lorax/users/guanxiang/2_veo/liang.veo/input/database/Pfam-A.hmm ${i}/${i}.contig.3k.orf.fa
awk 'BEGIN{print "target""\t""acc""\t""tlen""\t""query""\t""qlen""\t""evalue"}{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' ${i}/hmmout | sed '/\#/d' > ${i}/${i}.hmm.out
done



for i in D*
do
bowtie2-build ${i}/${i}.contig.fa ${i}/${i}.contig
bowtie2 -t -q  -x ${i}/${i}.contig -1 /media/lorax/users/guanxiang/1_igram/1_virome/analysis/3_vlp_vliadtion/sunbeam/sunbeam_output/qc/decontam/${i}_1.fastq.gz -2 /media/lorax/users/guanxiang/1_igram/1_virome/analysis/3_vlp_vliadtion/sunbeam/sunbeam_output/qc/decontam/${i}_2.fastq.gz  -S ${i}/dna.sam -p 4
samtools view -bS ${i}/dna.sam > ${i}/nonsorted.bam
samtools sort ${i}/nonsorted.bam -o ${i}/${i}.bam
rm ${i}/dna.sam
rm ${i}/nonsorted.bam
samtools index ${i}/${i}.bam
grep ">" ${i}/${i}.contig.3k.orf.fa | awk -F\# '{print $1"\t"$1"\t"$2"\t"$3}' | sed -e 's/\ //g' | sed -e 's/>//g' |awk '{gsub("_\\S+","",$2); print}' | awk '$5="+"' | awk 'BEGIN{print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${i}/${i}.orf.saf
featureCounts -a ${i}/${i}.orf.saf -F SAF -o ${i}/${i}.orf.count ${i}/${i}.bam -T 1
done





