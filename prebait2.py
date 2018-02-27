#!/usr/bin/env python3
# import os
# import subprocess
# import sys

# from Bio import SeqIO
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





param_path = '/Users/josec/Desktop/git_repos/PreBait/Example_paramfile.txt'
param_dict = param2dict(param_path)
print(param_dict['Threads'])


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

