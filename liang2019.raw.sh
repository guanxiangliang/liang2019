#!/bin/bash
#set -e
#set -x
#set -u

echo "###################################################"
echo "#This is the code for Liang et al 2019 manuscript #"
echo "#                                                 #"
echo "#By Guanxiang Liang, 2019                         #"
echo "#                                                 #"
echo "#guanxiang.liang@pennmedicine.upenn.edu           #"
echo "###################################################"


# 1. Preparation 

# 1.1 Setting up working directory
bash ./code/1.1.start.sh

# 1.2 Preparing database

## 1.2.1 Kraken custom database. Please follow instruction "https://github.com/zhaoc1/sunbeam_databases/blob/master/build_krakendb.sh". Save Krakendb in 'input/kraken.db'  
bash ${code}/1.2.1.sunbeam.extension.sh

## 1.2.2 Uniprot viral protein database. This is inspired by publication PMID:26489866.
bash ${code}/1.2.2.uniprot.data.sh

## 1.2.3 Target viral genomes
# Genomes from both RefSeq and Neighbor nucleotides represent 7 families, these genomes were downloaded by NCBI and saved into ${bao}/input/databbase/animal.virus.genome
#Adenoviridae (n = 464): https://www.ncbi.nlm.nih.gov/nuccore?term=Adenoviridae[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+AC_000001:AC_999999[pacc]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+(%22vhost+human%22[Filter])
#Polymaviridae (n =1224): https://www.ncbi.nlm.nih.gov/nuccore?term=Polyomaviridae[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+AC_000001:AC_999999[pacc]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+(%22vhost+human%22[Filter])
#Anelloviridae (n = 776): https://www.ncbi.nlm.nih.gov/nuccore?term=Anelloviridae[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+AC_000001:AC_999999[pacc]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+(%22vhost+human%22[Filter])
#Parvovirirdae (n = 60): https://www.ncbi.nlm.nih.gov/nuccore?term=Parvovirinae[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+AC_000001:AC_999999[pacc]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+(%22vhost+human%22[Filter])
#Caliciviridae (n = 1409): https://www.ncbi.nlm.nih.gov/nuccore?term=Caliciviridae[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+AC_000001:AC_999999[pacc]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+(%22vhost+human%22[Filter])
#Picornaviridae (n = 3728): https://www.ncbi.nlm.nih.gov/nuccore?term=Picornaviridae[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+AC_000001:AC_999999[pacc]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+(%22vhost+human%22[Filter])

## 1.2.4 Host genomes used for Sunbeam pipeline.
# dowload the lastest asembly from RefSeq
bash ${code}/1.2.4.download.refseq.sh 

##########################################################################################################################################################################################################################################################################

# 2. Virome data process and make plot for Fig1.B,C,D,E; Fig3A,B,C,D

# 2.1 Parepare wet.lab data and make Fig1B, 1C, 3D 
# Save metadata for both discovery and validation cohort in `input/meta.data`
# Save VLP particles, 16s qPCR, and adenovirus qPCR data in `input/wet.input`
# Save adeno qPCR data to `input/wet/input`
# Make sure all the libraries that required has been installed in terminal R
Rscript ${code}/2.1.vlp.qpcr.r

# 2.2 Process discovery raw reads and do kraken classification to make Fig 1D, E and do coverage analysis to plot Fig 3A,B,C
## 2.2.1 This step is for disvoery virome 
bash ${code}/2.2.1.discovery.virome.sh

## 2.2.2 This step is for validation virome
bash ${code}/2.2.2.validation.virome.sh

## 2.3 Use R to make Fig.1D,E; Fig.3A,B,C,D
Rscript  ${code}/2.3.virome.r


#######################################################################################################################################################################################################################################################################

# 3. Integrated analysis for virome, shotgun metagenome, wgs, and induced virome. To make Figure 2.

cd ${bao}

# 3.1 Save induced phage number in `/wet.data`

# 3.2 Process shotgun metagenome, wgs, and induced virome sequence by sunbeam. 
## 3.2.1 This step is for shotgun metagenome
bash ${code}/3.2.1.shotgun.metagenome.sh

## 3.2.2 This step is for wgs qc, and assembly
bash ${code}/3.2.2.wgs.sh

## 3.2.3 This step is for induced VLP sequencing
bash ${code}/3.2.3.induced.virome.sh

# 3.3 Integrated analysis to map stool vlp reads, induced VLP reads to wgs genome, and compare stool vlp and induced vlp
bash ${code}/3.3.vlp.shotgun.wgs.integrated.sh

# 3.4 To plot Fig.2B-G
Rscript ${code}/3.4.integrated.analysis.r


###########################################################################################################################################################################################################################################################################


# 4.1 
Rscript ${code}/4.1.supplementary.r




