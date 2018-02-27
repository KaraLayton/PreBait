#!/usr/bin/env python3
import os
import subprocess
# import sys

from Bio import SeqIO
# from Bio.Alphabet import SingleLetterAlphabet
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from multiprocessing import Pool
# from os.path import basename




def param2dict(param_path):
    """Read the parameter file and save the values in a dictionary for later use."""
    with open(param_path, mode = 'r') as param_handle:
        plines = param_handle.readlines()
        #Removes trailing white space and description snipet from param file
        param_list = [p.split('#', 1)[0].strip(' ') for p in plines]
        paramdict_keys = ['Threads','TxtPath','TargetPath', 'Genome', 'MinBaitLength', 'MinBaitDiv']
        param_dict = dict(zip(paramdict_keys, param_list))
    return param_dict

def run_command(command,command_out_path='none' change this!):
    command_output = subprocess.getstatusoutput(command)
    with open(command_out_path, "w") as out_handle:
        out_handle.write(command_output[1])
        if command_output[0]:
            print(f"{command} Failed")
            sys.exit(1)

def seqname_lister(seq_file):
    """Reads all of the names in a fasta file and returns a list"""
    records = SeqIO.parse(seq_file,"fasta")
    return [r.id for r in records]

def DNAorAA(seq_file):
    """Reads the first sequence in a file and guesses if DNA or RNA"""
    dna = set("NATGC")
    records = list(SeqIO.parse(seq_file,"fasta"))
    aminoacids = set(records[0].seq) - dna
    if aminoacids:
        alphabet = 'aa'
    else:
        alphabet = 'dna'
    return alphabet

def blaster(param_dict, query, blastdb_prefix):
    """Runs a blast search and returns the names of the hit sequences from the blast database."""
    alphabet = DNAorAA(query)
    if alphabet == 'dna':
        blast_prog = 'tblastx'
    elif alphabet == 'aa':
        blast_prog = 'tblastn'
    threads = param_dict["Threads"]
    #create blast database if there is none
    if not os.path.isfile(blastdb_prefix+".nin"):
        dbcommand = f"makeblastdb -in {blastdb_prefix} -parse_seqids -dbtype nucl"
        run_command(dbcommand)
    #Run Blast
    bcommand = f"{blast_prog} -query {query} -db {blastdb_prefix} -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads {threads}"
    run_command(bcommand,command_out_path='/Users/josec/Desktop/PyPreBait/Targets/Fast5_Chr_wes_dna_test.fa.txt')
    return





param_path = '/Users/josec/Desktop/git_repos/PreBait/Example_paramfile.txt'
param_dict = param2dict(param_path)
print(param_dict['Threads'])

name_duds=seqname_lister('/Users/josec/Desktop/NudiPreBait2/1_Generate_Target_Sets/AgalmaRuns/AgalmaCDSfast5/AgalmaOutput/Fast52genes.fa')[:10]
blaster(param_dict, '/Users/josec/Desktop/PyPreBait/Targets/Fast5_Chr_wes_dna_test.fa', param_dict["Genome"])

# class TargetSet:
#     """Takes a fasta file of target seqs (like agalma output) and a reference genome.
#     Blasts the targets to the genome and retrieves the names of the genome hits and the sequences """
#     def __init__(self,target_file,param_dict):
#         self.threads = param_dict["Threads"]
#         self.target_file = target_file
#         self.genome_file = param_dict["Genome"]

#     def targetNamelist(self):
#         self.target_file # fasta read file here
#         name_list="FillTHISIN"
#         return name_list
    
#     def genomeOrthologs(self):
#         pass

