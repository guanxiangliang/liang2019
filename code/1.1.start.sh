git clone -b stable https://github.com/sunbeam-labs/sunbeam ./code/sunbeam-stable
cd ./code/sunbeam-stable
./install.sh


git clone https://github.com/guanzhidao/sbx_dedup extensions/sbx_dedup
conda activate sunbeam
conda install --file extensions/sbx_dedup/requirements.txt
cd ../../
conda deactivate

conda activate liang2019

mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

echo "export bao=$PWD" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "export input=$PWD/input">> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "export code=$PWD/code" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "unset bao" >>$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo "unset input" >>$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo "unset code" >>$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh




