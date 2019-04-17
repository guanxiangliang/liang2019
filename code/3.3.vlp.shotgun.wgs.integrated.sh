# To make mapping files
bac="wgs25";ind="DNM3951";vir="D3951";ind="DM3951";shot="s3951"
cd ${genome}
mkdir map
cd map
bowtie2 -t -q  -x ${genome}/${bac}/contigs/${bac}.sorted.contig -1 ${iread}/${ind}_1.fastq.gz -2 $iread/${ind}_2.fastq.gz -S ${ind}.sam -p 8
samtools view -bS ${ind}.sam > ${ind}.bam
samtools sort ${ind}.bam -o ${ind}.sorted.bam
rm ${ind}.sam
rm ${ind}.bam
mv ${ind}.sorted.bam ${ind}.bam
samtools index ${ind}.bam
bamCoverage -b ${ind}.bam --normalizeUsing RPKM --binSize 100 -e=200 --outFileFormat bedgraph -o ${ind}.bg

bowtie2 -t -q  -x ${genome}/${bac}/contigs/${bac}.sorted.contig -1 ${vread}/${vir}_1.fastq.gz -2 $vread/${vir}_2.fastq.gz -S ${vir}.sam -p 8
samtools view -bS ${vir}.sam > ${vir}.bam
samtools sort ${vir}.bam -o ${vir}.sorted.bam
rm ${vir}.sam
rm ${vir}.bam
mv ${vir}.sorted.bam ${vir}.bam
samtools index ${vir}.bam
bamCoverage -b ${vir}.bam --normalizeUsing RPKM --binSize 100 -e=200 --outFileFormat bedgraph -o ${vir}.bg

## To prepare Fig2G,F
cd ~/project/igram/code_result
grep "Mito" tmp.csv | awk -F, '{print $4"\t"$5}' | sed 's/"//g' > /media/lorax/users/guanxiang/1_igram/1_virome/analysis/1_vlp_analysis/my.analysis/induced.vlp.to.stool.vlp/data.txt
grep "inducer" tmp.csv | awk -F, '{print $4"\t"$5}' | sed 's/"//g' > /media/lorax/users/guanxiang/1_igram/1_virome/analysis/1_vlp_analysis/my.analysis/induced.vlp.to.stool.vlp/data.inducer.txt

cd /media/lorax/users/guanxiang/1_igram/1_virome/analysis/1_vlp_analysis/my.analysis/induced.vlp.to.stool.vlp/

# 1. grep the sequence of induced VLP reads that can map to identified viral contigs
cat data.txt | while read  p
do
ind=$(echo ${p}|awk '{print $2}')
base=$(echo ${ind}|awk -F. '{print $1}')
if [ ! -d "$base" ]; then
mkdir ${base}
fi
echo ${ind}| python ~/software/python-tools/pick.fasta.pipe.py /media/lorax/users/guanxiang/1_igram/1_virome/analysis/2_induction_analysis/my.analysis/contigs/${base}/${base}.contig.sorted.fa >> ${base}/${base}.viral.contig.fa
bowtie2-build ${base}/${base}.viral.contig.fa ${base}/${base}
done

vlp="/media/lorax/users/guanxiang/1_igram/1_virome/analysis/1_vlp_analysis/sunbeam/sunbeam_output/qc/decontam"
for i in *
do
for j in *
do
base=$(echo ${j}| sed 's/^..//')
bowtie2 -t -q  -x ${i}/${i} -1 ${vlp}/"D"${base}_1.fastq.gz -2 ${vlp}/"D"${base}_2.fastq.gz  -S ${i}/dna.sam -p 4
samtools view -bS ${i}/dna.sam > ${i}/nonsorted.bam
samtools sort ${i}/nonsorted.bam -o ${i}/"D"${base}.${i}.bam
rm ${i}/dna.sam
rm ${i}/nonsorted.bam
samtools index ${i}/"D"${base}.${i}.bam
samtools idxstats ${i}/"D"${base}.${i}.bam > ${i}/"D"${base}.${i}.reads.txt
rm ${i}/*.bam*
done
done


