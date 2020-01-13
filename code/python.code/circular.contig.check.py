#! python3
'''
$ python circular.contig.py  original.fasta
'''

from Bio import SeqIO
import sys


print("contig.id"+"\t"+"length"+"\t"+"circularity")
for contig in  SeqIO.parse(open(sys.argv[1]), 'fasta'):
	cir = False
	length = len(contig.seq)
	#if length < 500: 
	#	print(contig.id +"\t"+ str(length)+"\t"+"short")
#	for k in range(10,1000):
	for k in range(10,150):
		if contig.seq[0:k] == contig.seq[length-k:]:
			cir = True
			break
	if cir == True:
		print(contig.id+"\t"+str(length)+"\t"+"circular")
	else:
		print(contig.id+"\t"+str(length)+"\t"+"non-circular")
