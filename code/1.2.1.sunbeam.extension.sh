#! /bin/bash

cd ${bao}
git clone -b stable https://github.com/eclarke/sunbeam ${code}/sunbeam-stable
cd ${code}/sunbeam-stable
bash install.sh
git clone https://github.com/guanzhidao/sbx_dedup extensions/sbx_dedup
conda activate sunbeam
conda install --file extensions/sbx_dedup/requirements.txt
conda deactivate 
cd ${bao}



