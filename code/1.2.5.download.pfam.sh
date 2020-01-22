cd input/database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/database_files/pfamA_ncbi.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/database_files/pfamA_tax_depth.txt.gz
gunzip pfamA_ncbi.txt.gz
gunzip pfamA_tax_depth.txt.gz
cd ${bao}
