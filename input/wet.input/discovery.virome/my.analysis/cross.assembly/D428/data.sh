for i in *.contig.fa
do
base=$(echo ${i}|awk -F. '{print $1}')
#vsearch --sortbylength ${i} --relabel ${base}.contig  --minseqlength 0 --maxseqlength -1 --output ${base}.contig.fa
#vsearch --sortbylength ${base}.contig.fa  --minseqlength 500 --maxseqlength -1 --output ${base}.contig.3k.fa
#rm ${base}-contigs.fa
bowtie2-build ${base}.contig.fa ${base}.contig
for j in /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/D*_1.fastq.gz
do
k=$(basename ${j}|awk -F_ '{print $1}')
bowtie2 -x ${base}.contig -1 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/${k}_1.fastq.gz -2 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/${k}_2.fastq.gz  -S good.sam -p 4
samtools view -bS good.sam > nonsorted.bam
samtools sort nonsorted.bam -o sorted.bam
samtools index sorted.bam
samtools idxstats sorted.bam > ${base}.${k}.read.txt
rm good.sam
rm *.bam
done

#bowtie2 -x ${base}.contig -1 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/dnc_1.fastq.gz -2 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/dnc_2.fastq.gz  -S good.sam -p 4
#samtools view -bS good.sam > nonsorted.bam
#samtools sort nonsorted.bam -o ${base}.nc.bam
#samtools index ${base}.nc.bam
#samtools idxstats ${base}.nc.bam > ${base}.nc.read.txt
#rm good.sam
#rm *.bam
#bowtie2 -x ${base}.contig -1 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/${base}1_1.fastq.gz -2 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/${base}1_2.fastq.gz  -S good.sam -p 4
#samtools view -bS good.sam > nonsorted.bam
#samtools sort nonsorted.bam -o ${base}.1.bam
#samtools index ${base}.1.bam
#samtools idxstats ${base}.1.bam > ${base}.1.read.txt
#rm good.sam
#rm *.bam
#bowtie2 -x ${base}.contig -1 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/${base}4_1.fastq.gz -2 /home/guanxian/4_phage.bac.infant/phage/discovery_decontam/${base}4_2.fastq.gz  -S good.sam -p 4
#samtools view -bS good.sam > nonsorted.bam
#samtools sort nonsorted.bam -o ${base}.4.bam
#samtools index ${base}.4.bam
#samtools idxstats ${base}.4.bam > ${base}.4.read.txt
#rm good.sam
#rm *.bam
#prodigal -i ${base}.contig.3k.fa  -a ${base}.contig.3k.orf.fa -p meta
#blastp -num_threads 4 -evalue 1e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" -query ${base}.contig.3k.orf.fa -db /home/guanxian/liang2019/input/database/uniprot.virome/uniprot.virome.db -out ${base}.3k.uniprot
#grep ">" ${base}.contig.3k.orf.fa | sed 's/>//g' | awk '{print $1}' > ${base}.orf.number
done




