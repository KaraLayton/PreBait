#!/usr/bin/python
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
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
	chr_exon_path ="%s/Chr_exon_names.txt"%(base_path)
#make out-directory
if not os.path.exists(seq_path):
	os.mkdir(seq_path)
if not os.path.exists(aln_path):
	os.mkdir(aln_path)

#Load the reciprical blast file into memory
with open(file_path) as file_handle:
	blines=[x.strip('\n') for x in file_handle.readlines()]
#Load list of all Chr exons to see if there was a Rblast hit
with open(chr_exon_path) as chr_exon_handle:
	all_chr_exons=[x.strip('\n') for x in chr_exon_handle.readlines()]



Fil_Rblast_hits=[]
throw_out=[]
#parse each output
for exon in blines:
	raw_blast=exon.split(',')
	exonerate_string,qseq,sseqid,sseq,p_id=raw_blast
	#just need the exon and the seq name for this step.
	exon_name,qseqid=exonerate_string.split(':')[0:2]
	#Filter if too diverget or identical between Chrom spp.
	if float(p_id)>90 and float(p_id)!=100:
		Fil_Rblast_hits.append(exon)
		with open("%s/%s.fasta"%(seq_path,exon_name),"w")as seq_handle:
			seq_handle.write(">%s\n%s\n>%s\n%s\n"%(qseqid,qseq,sseqid,sseq))
	else:
		throw_out.append(exon_name)

#Align with mafft
outname = "%s%s/%s.fasta" % (path, outpathname, exon)
        SeqIO.write(exon_seq_list, outname, "fasta")
        mafft_cline = MafftCommandline(input=outname)
        stdout, stderr = mafft_cline()
        with open("%s%s/aligned/%s_aligned.fasta" % (path, outpathname, exon),
                  "w") as handle:
            handle.write(stdout)
        align = AlignIO.read("%s%s/aligned/%s_aligned.fasta" %
                             (path, outpathname, exon), "fasta")



#Which seqs were not in the blast results? Pull out the singletons and write to txt file
singleton_count=0
with open("singletons.txt",'w') as singleton_file:
	for exon_name in all_chr_exons:
		if exon_name not in Fil_Rblast_hits and exon_name not in throw_out:
			singleton_file.write("%s\n"%(exon_name))
			singleton_count+=1



#Print results
print "Number of All exons %d"%len(all_chr_exons)
print "Number of filtered rblast hits %d"%len(Fil_Rblast_hits)
print "Number of exons filtered %d"%len(throw_out)
print "Number of singletons %d"%singleton_count
