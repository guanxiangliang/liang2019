## To qc wgs reads
mkdir ${input}/wet.input/wgs
cd ${input}/wet.input/wgs
mkdir data_fp
# save or download here
source activate sunbeam
sunbeam init --data_fp data_fp ./
sed -i 's/fwd_adapters:.*/fwd_adapters: []/' ./sunbeam_config.yml
sed -i 's/rev_adapters:.*/rev_adapters: []/' ./sunbeam_config.yml
sed -i "s#host_fp.*#host_fp: \x27${input}/databae/sunbeam.host.genome\x27#" ./sunbeam_config.yml
sed -i "s#kraken_db_fp.*#kraken_db_fp: \x27${input}/database/kraken.db\x27#" ./sunbeam_config.yml
sed -i "s#genomes_fp:.*#genomes_fp: \x27${input}/database/anmimal.virus.genome\x27#" ./sunbeam_config.yml
sunbeam run -- --configfile ./sunbeam_config.yml all_dedup
sunbeam run -- --configfile ./sunbeam_config.yml all_decontam
rm -r ./sunbeam_output/qc/0*
rm -r ./sunbeam_output/qc/*clean*


## Assembly
pa="/media/lorax/users/guanxiang/1_igram/3_wgs/sunbeam_output_no_mask/assembly/spades"
reads="/media/lorax/users/guanxiang/1_igram/3_wgs/sunbeam_output_no_mask/qc/decontam"
ssp="/media/lorax/users/guanxiang/1_igram/3_wgs/sunbeam_output_no_mask/assembly/sspace"
for i in ${reads}/*_1.fastq.gz
do
base=$(basename ${i} | awk -F_ '{print $1}')
cd ${spa}
mkdir ${base}
# Spades
/home/guanxiang/software/SPAdes-3.12.0-Linux/bin/spades.py --careful -o ${base} -1 ${reads}/${base}_1.fastq.gz -2 ${reads}/${base}_2.fastq.gz
cd ${base}
find . ! -name 'scaffolds.fasta' -type f -exec rm -f {} + # remove other files
find ./*  -type d -exec rm -r {} + # remove other folders
# SSPACE
cd ${ssp}
mkdir ${base}
cd ${base}
echo -e "${base}\t${reads}/${base}_1.fastq.gz\t${reads}/${base}_2.fastq.gz\t550\t1\tFR" >> lib.txt
perl /home/guanxiang/software/sspace_basic/SSPACE_Basic_v2.0.pl -l lib.txt -s ${spa}/${base}/scaffolds.fasta -b ${base} -T 8
find . ! -name '*scaffolds.fasta' -type f -exec rm -f {} + # remove other files
find ./*  -type d -exec rm -r {} + # remove other folders
done


## CheckM
cd sspace
for i in wgs*
do
mkdir ${i}/checkm
mkdir ${i}/checkm/bin
cp ${i}/*.fasta ${i}/checkm/bin
mv ${i}/checkm/bin/*.fasta ${i}/checkm/bin/${i}.fa
checkm lineage_wf -t 8 -x fa ${i}/checkm/bin ${i}/checkm/out
id+="${i}\n"  #add each id to this $id with \n
echo -e ${id} > tmp1
com+="$(awk -FCompleteness '{print $2}' ${i}/checkm/out/storage/bin_stats_ext.tsv | awk -FGCN '{print $1}' | awk -F\  '{print $2}' | cut -c1-5)\n"
echo -e ${com} > tmp2
con+="$(awk -FContamination '{print $2}' ${i}/checkm/out/storage/bin_stats_ext.tsv | awk -F\# '{print $1}' | awk -F\  '{print $2}'| cut -c1-3)\n"
echo -e ${con} > tmp3
done
sed -i "1i id" tmp1
sed -i "1i completeness" tmp2
sed -i "1i contamination" tmp3
paste tmp1 tmp2 tmp3 > checkm.summary
rm tmp*
unset id
unset com
unset con


