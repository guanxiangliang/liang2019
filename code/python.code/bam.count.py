#! /usr/bin/python3

import pysam
import numpy as np
import Bio 
from Bio import SeqIO
import sys

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
#for read in samfile.fetch('D333.contig1'):
#	print(read)
#for x in  samfile.pileup('D333.contig1'):
#	print(x)
contigs = {}
mapped_set = set()
total_set = set()

def read_good(bam):
	#min_mapid = 0
	min_mapq = 10
	#min_reafq = 0
	#min_aln_cov = 0
	#align_len = len(bam.query_alignment_sequence)
	#query_len = bam.query_length
	#if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < min_mapid:
	#    return False
	#if np.mean(bam.query_qualities) < min_reafq:
	#    return False
	if bam.mapping_quality < min_mapq:
		return False
	#elif align_len/float(query_len) < min_aln_cov:
	#    return False
	else:
		return True

def read_bad(bam):
	max_mapq = 10
	if bam.mapping_quality < max_mapq:
		return True
	else:
		return False

#count=samfile.count("D333.contig2",read_callback="nofilter")

for read in samfile.fetch(until_eof=True):
	total_set.add(read.query_name)

for read in samfile.fetch(until_eof=False):
        mapped_set.add(read.query_name)

for read in Bio.SeqIO.parse(open(sys.argv[2]), 'fasta'):
	contigs[read.id] = str(read.seq)

print('\t'.join(["contig.id","contig.length","mapped","mapped.low.qc","total.mapped","total"]))

for contig_id in list(contigs.keys()):
	contig_length = len(contigs[contig_id])
	count=samfile.count(contig_id,read_callback=read_good)
	countbad=samfile.count(contig_id,read_callback=read_bad)
	print('\t'.join([contig_id, str(contig_length), str(count),str(countbad),str(len(mapped_set)*2),str(len(total_set)*2)]))




