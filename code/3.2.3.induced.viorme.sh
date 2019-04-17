## This step follow default sunbeam to process induced vlp sequence
mkdir ${input}/wet.input/induced.virome
cd ${input}/wet.input/induced.virome
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
##


wd="/media/lorax/users/guanxiang/1_igram/1_igram_virome_hiseq/analysis/2_induction_analysis/my.analysis"
reads="/media/lorax/users/guanxiang/1_igram/1_igram_virome_hiseq/analysis/2_induction_analysis/sunbeam_output/qc/decontam"
for i in ${wd}/contigs/*; do base=$(basename ${i});
mkdir ${wd}/orfs/${base}
vsearch --sortbylength ${wd}/contigs/${base}/${base}.contig.fa  --relabel ${base}.contig  --minseqlength 3000 --maxseqlength -1 --output ${wd}/orfs/${base}/${base}.contig.3k.fa
prodigal -i ${wd}/orfs/${base}/${base}.contig.3k.fa  -a ${wd}/orfs/${base}/${base}.contig.3k.orf.fa -p meta
blastp -num_threads 8 -evalue 1e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" -query ${wd}/orfs/${base}/${base}.contig.3k.orf.fa -db ~/database/virome_database/uniprot_virome/uniprot_virome -out ${wd}/orfs/${base}/${base}.3k.uniprot
grep ">" ${wd}/orfs/${base}/${base}.contig.3k.orf.fa | sed 's/>//g' | awk '{print $1}' > ${wd}/orfs/${base}/${base}.orf.number 
