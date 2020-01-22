#! /bin/bash


# 2.2.2 This step is for validation virome
mkdir ${bao}/input/wet.input/validation.virome
cd ${bao}/input/wet.input/validation.virome
mkdir data_fp
# save or download here
source activate sunbeam
sunbeam init --data_fp data_fp ./
sed -i 's/fwd_adapters:.*/fwd_adapters: []/' ./sunbeam_config.yml
sed -i 's/rev_adapters:.*/rev_adapters: []/' ./sunbeam_config.yml
sed -i "s#host_fp.*#host_fp: \x27${input}/databae/sunbeam.host.genome\x27#" ./sunbeam_config.yml
sed -i "s#kraken_db_fp.*#kraken_db_fp: \x27${input}/database/kraken.db\x27#" ./sunbeam_config.yml
sed -i "s#genomes_fp:.*#genomes_fp: \x27${input}/database/anmimal.virus.genome\x27#" ./sunbeam_config.yml
cat $SUNBEAM_DIR/extensions/sbx_dedup/config.yml >> ./sunbeam_config.yml
sunbeam run -- --configfile ./sunbeam_config.yml all_dedup
sunbeam run -- --configfile ./sunbeam_config.yml all_decontam
sunbeam run -- --configfile ./sunbeam_config.yml all_mapping
rm -r ./sunbeam_output/qc/0*
rm -r ./sunbeam_output/qc/*clean*
cd ${bao}
# To process mapping data to obtain the animal virus coverage data
cd ${input}/input/wet.input/validation.virome/sunbeam_output/mapping
rm -r intermediates
rm adeno/R*
rm anello/R*
rm parvo/R*
rm polymavi/R*
rm astro/D*
rm calici/D*
rm picona/D*
for j in *
do
        cd ${j}
        for i in *.bam
        do
        base=$(echo ${i} |awk -F. '{print $1}')
        samtools index ${i}
        bedtools genomecov -ibam ${base}.bam |grep -v  "genome" > ${base}.coverage
        python code/python.code/couht.bam.py ${base}.bam > ${base}.read
	#samtools idxstats ${base}.bam > ${base}.read
        rm *bam*
        cd ..
        done
done
cd ${bao}
mkdir ${input}/wet.input/coverage
mkdir ${input}/wet.input/coverage/vali
cp -r ${input}/wet.input/validation.virome/sunbeam_output/mapping/* ${input}/wet.input/coverage/vali








