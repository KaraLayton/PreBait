#!/usr/bin/python
from Bio import SeqIO
# import os
import re
import sys

usage="""
BetterBest.py --exdir Exonerate/Output/Directory
This script parses the exonerate results from blasting the query(aplysia exons)
to target (Nudi transcriptomes). This is to map the (assumed) introns in the CDS
from the transcriptomes. The script then picks which transcriptome sequence had
the better hit to each aplysia target (gene) and only saves this result. These
'better hits' will then be used as a query to reciprically blast the other
transcriptomes to find all of the hits for a given target. This works better than
comparing the results when aplysia was the query because it is too divergent and
all hits are not captured in the filtered exonerate results.

"""
####Set variables
organisms = ['Chr_wes', 'Chr_mag']
outpathname = "BetterChrTopHit"
prefix="exons_"
#Set length threshold and percent ID threshold. Anything below these thresholds will be thrown out.
L_threshold = int("120")


args = sys.argv[1:]
#Print usage if no input is put in by user
if not args:
    print usage
    sys.exit(1)
if args[0] == '--exdir':
    path = args[1]


input_dictionary = {}
naughty_list = []
dup_tracker_list = []
for organism in organisms:
    # input_seq_iterator reads fasta sequences
    for record in SeqIO.parse("%s%s%s.fasta" % (path,prefix,organism), "fasta"):
        #Dont add exon to the record_dict if it is too short
        if len(record) < 120:
            continue
        # Pulls out the score, length, sequence and qlength from names in output
        # All of this info came from blasting aplysia (query) exons to target transcriptome
        # qlength is query exon length, length is target exon length
        record_dict = {}
        smatch = re.search(r"SS(\d*)", record.name)
        lmatch = re.search(r"LL(\d*)/(\d*)", record.name)
        tmatch = re.search(r"TT(.*):HH", record.name)
    	exon_name = record.name.split(":")[0]
        txt_seq = str(tmatch.group(1))
        record_dict["exon_name"] = exon_name
        record_dict["species"] = organism
        record_dict["score"] = int(smatch.group(1))
        record_dict["length"] = int(lmatch.group(1))
        record_dict["qlength"] = int(lmatch.group(2))
        record_dict["seq"] = record
        # If the exon length is longer than the query(ie the real exon length)
        # then we know something is wrong, so it is thrown out later
        if record_dict["qlength"] < record_dict["length"]:
            naughty_list.append(exon_name)
        # Checks if exon for this organism is in dictionary already.
        if txt_seq in dup_tracker_list:
            continue
        else:
            # If exon not in dictionary add. If there is already one, pick which one has the higher score
            try:
                if input_dictionary[exon_name]["score"] < record_dict["score"]:
                    input_dictionary[exon_name] = record_dict
                    dup_tracker_list.append(txt_seq)
                else:
                    pass
            except:
                input_dictionary[exon_name] = record_dict
                dup_tracker_list.append(txt_seq)

# Seperate the best hits into the organisms then export them as fasta files
# Filters for exons longer than 120bp and makes sure that they were not flagged
# in the 'naughty_list' for being longer than expected.
Chr_wes_best = []
Chr_mag_best = []
for value in input_dictionary.values():
    if value['length'] > L_threshold and value['exon_name'] not in naughty_list:
        if value['species'] == 'Chr_wes':
            Chr_wes_best.append(value["seq"])
        else:
            Chr_mag_best.append(value["seq"])
outname = "%s/Chr_wes_best.fasta" % (path)
SeqIO.write(Chr_wes_best, outname, "fasta")
outname = "%s/Chr_mag_best.fasta" % (path)
SeqIO.write(Chr_mag_best, outname, "fasta")
