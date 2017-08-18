#!/usr/bin/python
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import re

target_set = "CleanLongAplyExons_"
path = "/Users/josec/Desktop/NudiSilicoTest/Exonerate/Outputs/"
exon_list_path = "/Users/josec/Desktop/NudiSilicoTest/Exonerate/Outputs/LongExons/Exon_list.txt"
# inputs=['Doris','Chr_wes','Chr_mag','Bathi','Act_QLD']
# inputs=['Doris','Chr_wes','Bathi','Act_QLD']
inputs = ['Chr_wes', 'Chr_mag']
outpathname = "Clean_exons"
L_threshold = int("115")
P_threshold = int("65")

with open(exon_list_path, "rU") as listhandle:
    Aply_exons = [x.strip(">\n") for x in listhandle.readlines()]

for exon in Aply_exons:
    exon_seq_list = []
    for input_string in inputs:
        inx = "%s%s" % (target_set, input_string)
        handle = open("%s%s.fasta" % (path, inx))
        input_seq_iterator = SeqIO.parse(handle, "fasta")
        for record in input_seq_iterator:
            if exon == record.name.split(":")[0]:
                match = re.search(r"PP(\d*)", record.name)
                score = int(match.group(1))
                lmatch = re.search(r"LL(\d*)", record.name)
                length = int(lmatch.group(1))
                if length > L_threshold and score > P_threshold:
                    exon_seq_list.append(record)
                    break

        handle.close()
    if len(exon_seq_list) > 1:
        outname = "%s%s/%s.fasta" % (path, outpathname, exon)
        SeqIO.write(exon_seq_list, outname, "fasta")
        mafft_cline = MafftCommandline(input=outname)
        stdout, stderr = mafft_cline()
        with open("%s%s/aligned/%s_aligned.fasta" % (path, outpathname, exon),
                  "w") as handle:
            handle.write(stdout)
        align = AlignIO.read("%s%s/aligned/%s_aligned.fasta" %
                             (path, outpathname, exon), "fasta")
