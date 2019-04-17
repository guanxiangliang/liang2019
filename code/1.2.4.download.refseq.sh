mkdir ./input/database/sunbeam.host.genome
cd ./input/database/sunbeam.host.genome
rm *
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/assembly_summary.txt 
tail -n 1 assembly_summary.txt | awk '{print $20"/*genomic.fna.gz"}' > download.list
rm assembly_summary.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/assembly_summary.txt
tail -n 1 assembly_summary.txt | awk '{print $20"/*genomic.fna.gz"}' >> download.list
rm assembly_summary.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Escherichia_virus_phiX174/assembly_summary.txt
tail -n 1 assembly_summary.txt | awk '{print $18"/*genomic.fna.gz"}' >> download.list
rm assembly_summary.txt
wget -R "*_from_genomic.fna.gz" -i download.list
mv *GRCh38* human.fasta.gz
mv *GRCm38* mouse.fasta.gz
mv *Viral* phix174.fasta.gz
gunzip *.gz
rm download.list
cd ${bao}



