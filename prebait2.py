#!/usr/bin/env python3
import os
import subprocess
import sys
import re


from Bio import SeqIO
# from Bio.Alphabet import SingleLetterAlphabet
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from pathlib import Path


#For Venn Diagrams
from matplotlib_venn import venn2, venn3
import numpy
import scipy
from matplotlib import pyplot as plt
# from os.path import basename




def param2dict(param_path):
    """Read the parameter file and save the values in a dictionary for later use."""
    with open(param_path, mode = 'r') as param_handle:
        plines = param_handle.readlines()
        #Removes trailing white space and description snipet from param file
        param_list = [p.split('#', 1)[0].strip(' ') for p in plines]
        paramdict_keys = ['Threads','TxtPath','GeneTarget_path', 'Genome', 'MinBaitLength', 'MinBaitDiv']
        param_dict = dict(zip(paramdict_keys, param_list))
    return param_dict

def run_command(command, command_out_path=None):
    """Runs command in subprocess. Prints error if commmand did not run. 
    if outpath given it will print output there, otherwise this function returns it. """
    command_output = subprocess.getstatusoutput(command)
    if command_output[0]:
        print(f"{command} Failed")
        sys.exit(1)
    if command_out_path:
        with open(command_out_path, "w") as out_handle:
            out_handle.write(command_output[1])
    return command_output[1]

def run_command_parallel(command, num_threads):
    p = Pool(num_threads)
    p.map(run_command, command)
    return 

def seqname_lister(seq_file):
    """Reads all of the names in a fasta file and returns a list"""
    records = SeqIO.parse(seq_file,"fasta")
    seqname_list = [r.id for r in records]
    return seqname_list

def DNAorAA(seq_file):
    """Reads a file and guesses if DNA or RNA.
    Exits with error message if invalid residues are present"""
    records = list(SeqIO.parse(seq_file,"fasta"))
    seqset = set(''.join([str(record.seq) for record in records]))
    dna = set("NATGC")
    nondna = seqset-dna
    if len(nondna) > 2:
        alphabet = 'aa'
    elif 'X' in nondna:
        print(f"ERROR:{seq_file} has '-' or 'X' in it. Please fix and run again")
        sys.exit(1)
    else:
        alphabet = 'dna'
    return alphabet

def blaster(param_dict, query, blastdb_prefix, blast_out_path = None):
    """Runs a blast search and returns the names of the hit sequences from the blast database.
    If out_path specified saves results to file, otherwise returns as string"""
    alphabet = DNAorAA(query)
    if alphabet == 'dna':
        blast_prog = 'tblastx'
    elif alphabet == 'aa':
        blast_prog = 'tblastn'
    threads = param_dict["Threads"]
    #create blast database if there is none
    if not Path(blastdb_prefix+".nin").exists():
        dbcommand = f"makeblastdb -in {blastdb_prefix} -parse_seqids -dbtype nucl"
        run_command(dbcommand)
    #Run Blast
    bcommand = f"{blast_prog} -query {query} -db {blastdb_prefix} -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads {threads}"
    blast_output = run_command(bcommand,command_out_path=blast_out_path)
    return blast_output

def geneset2genome_blast(param_dict, geneset):
    """ Blasts a geneset to a genome cds like this one:
    #http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/
    reformats the blast results to the refseq protein IDs and removes duplicate hits. The blast cutoff is e-20
    Output is written to file and returned by function
    """
    rawblast_output = blaster(param_dict, query = geneset, blastdb_prefix = param_dict["Genome"])
    blast_output = re.sub(r".*cds_(.*\..*)_.*",r"\1",rawblast_output)
    unique_blast_lines = list(set(blast_output.split("\n")))
    blast_outpath = Path(geneset).with_suffix('.txt')
    with open(blast_outpath,'w') as out_handle:
        out_handle.write("\n".join(unique_blast_lines))
    return unique_blast_lines

def isafasta(fasta_to_check_path):
    """Returns True if file is parsable as a fasta by BioPython"""
    records = SeqIO.parse(fasta_to_check_path,"fasta") 
    try:
        return any(records)
    except:
        return False

def venndiagramer(venndict, outfile_path):
    """Makes a venn diagram of dictionary values as sets with the keys as the names.
    Detects if there are 2 or three sets. Returns error if there are not. """
    vkeys = list(venndict.keys())
    vin = [set(v) for v in venndict.values()]
    if len(vkeys) == 2:
        venn2(vin,vkeys)
    elif len(vkeys) ==3:
        venn3(vin,vkeys)
    else:
        print("Cannot make venn diagram. Need 2 or 3 groups to compare")
        return 
    plt.savefig(str(outfile_path))
    return

def combine_genesets(param_dict):
    """Steps 1-3 in the prebait workflow.
    Takes gene target fasta sets and blasts them against the reference genome.
    The blast hits are saved as the refseq protein id. These protein Ids are 
    compared across the sets and then returned as a set. A venn
    diagram is produced comparing the overlap between the genesets. """
    geneset_path = Path(param_dict["GeneTarget_path"])
    #list of geneset file paths excluding hiden files
    paths = [str(p) for p in geneset_path.iterdir() if not p.name.startswith(".")]
    geneset_dict = {}
    for fa_path in paths:
        if isafasta(fa_path):
            set_name = Path(fa_path).stem
            cds_set = geneset2genome_blast(param_dict, geneset=fa_path)
            geneset_dict[set_name] = cds_set
    venndiagramer(geneset_dict, outfile_path=geneset_path.joinpath('Venn.pdf'))
    the_set = []
    for listi in geneset_dict.values():
        the_set += listi
    return set(the_set)






param_path = '/Users/josec/Desktop/git_repos/PreBait/Example_paramfile.txt'
param_dict = param2dict(param_path)

target_genes = combine_genesets(param_dict)
print(len(target_genes))
print(target_genes)

# set1 = set(['XP_005090798.1', 'XP_005111329.1', 'XP_005107541.1', 'XP_005109146.1', 'NP_001191497.1', 'XP_012936479.1'])
# print(set1)
# set2 = set(['XP_005106305.2', 'XP_005101496.1', 'XP_005100335.1', 'XP_005110818.1', 'XP_005099957.1', 'XP_005113254.1', 'XP_012940806.1', 'XP_005111045.2', 'XP_005093366.1', 'NP_001191541.1', 'XP_005100202.2', 'XP_005112795.1', 'XP_005094159.1', 'XP_005098277.1', 'XP_012935999.1', 'XP_012941843.1'])
# set3 = set(['C', 'D',' E', 'F', 'G'])
# names = ['Fast5_Chr_wes_dna_test', 'Nudis2_Chr_wes_aa_test', 'Set33']
# venn3([set1, set2, set3], names)
# print([set1, set2, set3])
# plt.savefig("/Users/josec/Desktop/PyPreBait/venn.pdf")



