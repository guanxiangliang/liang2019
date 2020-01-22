for i in *_1.fastq.gz
do
genome="adeno"
j=$(echo ${i}|awk -F_ '{print $1}')
#bwa mem  -t 4 /home/guanxian/liang2019/input/database/anmimal.virus.genome/${genome}.fasta  ${j}.1  ${j}.2 -o /home/guanxian/magic/sunbeam_output/qc/decontam/mapping/${genome}/${j}.good.sam
bowtie2 -x /home/guanxian/liang2019/input/database/anmimal.virus.genome/${genome}.fasta -1 ${j}_1.fastq.gz  -2 ${j}_2.fastq.gz  -S /home/guanxian/coverage/mapping/${genome}/${j}.good.sam -p 4
samtools view -Sb /home/guanxian/coverage/mapping/${genome}/${j}.good.sam > /home/guanxian/coverage/mapping/${genome}/${j}.nonSorted.bam
samtools sort /home/guanxian/coverage/mapping/${genome}/${j}.nonSorted.bam -o /home/guanxian/coverage/mapping/${genome}/${j}.Sorted.bam
samtools rmdup -s /home/guanxian/coverage/mapping/${genome}/${j}.Sorted.bam /home/guanxian/coverage/mapping/${genome}/${j}.bam
samtools index /home/guanxian/coverage/mapping/${genome}/${j}.bam
bedtools genomecov -ibam /home/guanxian/coverage/mapping/${genome}/${j}.bam |grep -v  "genome" > /home/guanxian/coverage/mapping/${genome}/${j}.coverage
python python/code/bam.count.py /home/guanxian/coverage/mapping/${genome}/${j}.bam > /home/guanxian/coverage/mapping/${genome}/${j}.read
#samtool idxstats /home/guanxian/coverage/mapping/${genome}/${j}.bam > /home/guanxian/coverage/mapping/${genome}/${j}.read
rm /home/guanxian/coverage/mapping/${genome}/${j}.good.sam
rm /home/guanxian/coverage/mapping/${genome}/${j}.nonSorted.bam
rm /home/guanxian/coverage/mapping/${genome}/${j}.Sorted.bam
rm /home/guanxian/coverage/mapping/${genome}/${j}.bam*
done
