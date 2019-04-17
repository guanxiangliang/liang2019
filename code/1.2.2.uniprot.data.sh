#! bash

mkdir ./input/database/uniprot.virome
cd ./input/database/uniprot.virome
rm *
wget  ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz # download uniprot tremble protein fasta database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz # downlaod uniprot sprot data
cat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > uniprot.fasta.gz # combine tremble and sprot 
gunzip uniprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.dat.gz # download uniprot tremble virus information
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz # download uniprot sprot virus information
cat uniprot_trembl_viruses.dat.gz uniprot_sprot_viruses.dat.gz > uniprot.viruses.dat.gz # combine
gunzip uniprot.viruses.dat.gz
grep '^ID' uniprot.viruses.dat | awk '{print $2}' > virus.list.txt # grep the accession id of virus, and get ready to pull out viral protein seq from uniprot protein fasta database
awk '{print $1}' uniprot.fasta  |  awk -F\| '{if($1==">tr"|| $1==">sp") print ">"$3; else print $1 }' > tmp.fa # process the header of fasta and ready for seqtk to use for subset
seqtk subseq tmp.fa virus.list.txt > uniprot.virome.fasta # seqtk is super fast to subset million reads
makeblastdb -dbtype prot -in uniprot.virome.fasta  -out uniprot.virome.db # to make the blast database
egrep "^(OX|ID)" ./uniprot.viruses.dat | sed -e :a -e '$!N;s/\n\(OX\)/ /;ta' -e 'P;D'| awk '{print $6}' | awk -F= '{print $2}' > tmp.tax # to get the taxon id
egrep "^(DE|ID)" ./uniprot.viruses.dat | sed -e :a -e '$!N;s/\n\(DE\)/ /;ta' -e 'P;D' | awk '{print $2}' > tmp.id # to get the uniprot id
egrep "^(DE|ID)" ./uniprot.viruses.dat | sed -e :a -e '$!N;s/\n\(DE\)/ /;ta' -e 'P;D' | awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$10"\t"$12}' | awk -F= '{print $2}' | awk -F{ '{print $1}' | sed s/"\t"$//g | sed s/"\t"/_/g > tmp.name # to get the uniprot protein name
paste tmp.id tmp.tax tmp.name > uniprot.virome.meta # to merge them to each other
rm tmp*
rm *.gz
rm *.dat
rm virus.list.txt
rm uniprot.fasta





