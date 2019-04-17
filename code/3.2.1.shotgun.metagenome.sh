mkdir ${input}/wet.input/shotgun.metagenome
cd ${input}/wet.input/shotgun.metagenome
mkdir data_fp
# save or download here
source activate sunbeam
sunbeam init --data_fp data_fp ./
sed -i 's/fwd_adapters:.*/fwd_adapters: []/' ./sunbeam_config.yml
sed -i 's/rev_adapters:.*/rev_adapters: []/' ./sunbeam_config.yml
sed -i "s#host_fp.*#host_fp: \x27${input}/databae/sunbeam.host.genome\x27#" ./sunbeam_config.yml
sed -i "s#kraken_db_fp.*#kraken_db_fp: \x27${input}/database/kraken.db\x27#" ./sunbeam_config.yml
cat $SUNBEAM_DIR/extensions/sbx_dedup/config.yml >> ./sunbeam_config.yml
sunbeam run -- --configfile ./sunbeam_config.yml all_dedup
sunbeam run -- --configfile ./sunbeam_config.yml all_decontam
sunbeam run -- --configfile ./sunbeam_config.yml all_classify
rm -r ./sunbeam_output/qc/0*
rm -r ./sunbeam_output/qc/*clean*

cd ${bao}
