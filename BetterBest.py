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
outpathname = "BetterChrTopHit"

with open(exon_list_path, "rU") as listhandle:
    Aply_exons = [x.strip(">\n") for x in listhandle.readlines()]

input_dictionary = {}
naughty_list = []
dup_tracker_list = []
for input_string in inputs:
    inx = "%s%s" % (target_set, input_string)
    handle = open("%s%s.fasta" % (path, inx))
    input_seq_iterator = SeqIO.parse(handle, "fasta")

    for record in input_seq_iterator:
        record_dict = {}
        smatch = re.search(r"SS(\d*)", record.name)
        lmatch = re.search(r"LL(\d*)/(\d*)", record.name)
        tmatch = re.search(r"(TTChr_..._\d*)", record.name)
        exon_name = record.name.split(":")[0]
        txt_seq = str(tmatch.group(1))
        record_dict["exon_name"] = exon_name
        record_dict["species"] = input_string
        # record_dict["txt_seq"]=txt_seq
        record_dict["score"] = int(smatch.group(1))
        record_dict["length"] = int(lmatch.group(1))
        record_dict["qlength"] = int(lmatch.group(2))
        record_dict["seq"] = record
        if record_dict["qlength"] < record_dict["length"]:
            naughty_list.append(exon_name)
        if txt_seq in dup_tracker_list:
            pass
        else:
            #if exon not in dictionary add. If there is already one, pick which one has the higher score
            try:
                if input_dictionary[exon_name]["score"] < record_dict["score"]:
                    input_dictionary[exon_name] = record_dict
                    dup_tracker_list.append(txt_seq)
                else:
                    pass
            except:
                input_dictionary[exon_name] = record_dict
                dup_tracker_list.append(txt_seq)
    handle.close()
Chr_wes_best = []
Chr_mag_best = []
for value in input_dictionary.values():
    if value['length'] > 120 and value['exon_name'] not in naughty_list:
        if value['species'] == 'Chr_wes':
            Chr_wes_best.append(value["seq"])
        else:
            Chr_mag_best.append(value["seq"])
outname = "%s%s/Chr_wes_best.fasta" % (path, outpathname)
SeqIO.write(Chr_wes_best, outname, "fasta")
outname = "%s%s/Chr_mag_best.fasta" % (path, outpathname)
SeqIO.write(Chr_mag_best, outname, "fasta")
print len(naughty_list)
