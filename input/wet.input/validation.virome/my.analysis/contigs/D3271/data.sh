for j in *-contigs.fa
do
i=$(echo $j|awk -F- '{print $1}' )
vsearch --sortbylength ${i}-contigs.fa  --relabel ${i}.contig   --output  ${i}.contig.fa
vsearch --sortbylength ${i}.contig.fa  --minseqlength 3000 --maxseqlength -1 --output  ${i}.contig.3k.fa
prodigal -i ${i}.contig.3k.fa  -a ${i}.contig.3k.orf.fa -p meta
hmmscan --domtblout hmmout /home/guanxian/pfam/Pfam-A.hmm ${i}.contig.3k.orf.fa
awk 'BEGIN{print "target""\t""acc""\t""tlen""\t""query""\t""qlen""\t""evalue"}{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' hmmout | sed '/\#/d' > ${i}.hmm.out
rm hmmout
bowtie2-build ${i}.contig.fa ${i}.contig
bowtie2 -t -q  -x ${i}.contig -1 /home/guanxian/magic3/sunbeam_output_remove/qc/decontam/${i}_1.fastq.gz -2 /home/guanxian/magic3/sunbeam_output_remove/qc/decontam/${i}_2.fastq.gz  -S dna.sam -p 4
samtools view -bS dna.sam > nonsorted.bam
samtools sort nonsorted.bam -o ${i}.bam
rm dna.sam
rm nonsorted.bam
samtools index ${i}.bam
grep ">" ${i}.contig.3k.orf.fa | awk -F\# '{print $1"\t"$1"\t"$2"\t"$3}' | sed -e 's/\ //g' | sed -e 's/>//g' |awk '{gsub("_\\S+","",$2); print}' | awk '$5="+"' | awk 'BEGIN{print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${i}.orf.saf
featureCounts -a ${i}.orf.saf -F SAF -o ${i}.orf.count ${i}.bam -T 4
rm ${i}.bam*
rm ${i}*bt2
done




