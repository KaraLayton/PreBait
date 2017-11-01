#!/usr/bin/python
from Bio import SeqIO
import re
import os,sys

usage="""
RBlast_Parser.py --btxt Chr_Chr_blast.txt

parses output from command:
blastn -task blastn -query Chr_mag_best.fasta -db blastdb/Chr_wes_trinity.fasta -evalue 1e-10 -outfmt '10 qseqid qseq sseqid sseq pident' -max_target_seqs 1 -max_hsps 1 -num_threads $threads >> Chr_Chr_blast.txt

This script will take the results from the reciprical blast and make an exon alignment for each exon targeted.
Blast results that were identical were filtered because they are not informative between Chr spp.
Additionally if the blast results were <92_percent identical they were not included in the resulting alignments because they are probably paralogs.
"""
args = sys.argv[1:]
#Print usage if no input is put in by user
if not args:
	print usage
	sys.exit(1)
if args[0] == '--btxt':
	file_path = os.path.abspath(args[1])
	base_path = os.path.dirname(file_path)
	seq_path = "%s/exon_seqs"%(base_path)
	aln_path = "%s/exon_aln"%(base_path)
#make out-directory
if not os.path.exists(seq_path):
	os.mkdir(seq_path)
if not os.path.exists(aln_path):
	os.mkdir(aln_path)

#Load the reciprical blast file into memory
with open(file_path) as file_handle:
	blines=[x.strip('\n') for x in file_handle.readlines()]

filtercount=0
exon_name_list=[]
#parse each output
for exon in blines:
	raw_blast=exon.split(',')
	exonerate_string,qseq,sseqid,sseq,p_id=raw_blast
	#just need the exon and the seq name for this step.
	exon_name,qseqid=exonerate_string.split(':')[0:2]
	#Filter if too diverget or identical between Chrom spp.
	if float(p_id)>90 and float(p_id)!=100:
		exon_name_list.append(exon_name)
		pass
	else:
		filtercount+=1

with open("exon_names.txt",'w') as exon_name_file:
	for exon_name in exon_name_list:
		exon_name_file.write("%s\n"%(exon_name))
print filtercount
print len(blines)-filtercount
