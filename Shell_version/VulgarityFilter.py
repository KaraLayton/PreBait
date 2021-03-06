#!/usr/bin/env python
import re
import sys

usage = """
Usage: VulgarityFilter.py --in exonerate_outfile.txt

Parses the output of exonerate with: --ryo '%qi\t%pi\t%qas\t%V\tEND\n'
output: inputname.fa of a fasta with each exon listed as a seperate sequence IF the exon is >200bp.
example exonerate command:
exonerate --model est2genome Test44.fasta /Users/josec/Desktop/NudiSilicoTest/Exonerate/acl_ref_AplCal3.0_chrUn.fa -Q DNA -T DNA --showvulgar F --showalignment F --percent 90 --verbose 0 --ryo '%qi\t%pi\t%qas\t%V\tEND\n' --fsmmemory 20G --bestn 1 > exonerate_outfile.txt


This chops up the sequences in the exonerate_outfile into exons and saves them (if > 200bp) as a fasta with headers:
>Gene_exonNumber
EXONSEQUENCE

"""


def openfile(infile):
    # Read exonerate output and make list of output for each target (gene)
    # The results are stored in target_list
    with open(infile, 'Ur') as inhandle:
        target_list = inhandle.read().split('END')
    return target_list


def parser(target):
    # This takes the results and parses them into the target_dict
    target_dict = {}
    tlist = [x.replace('\n', '') for x in target.split('\t')]
    target_dict['Target'] = tlist[0]
    target_dict['Percent'] = tlist[1]
    target_dict['CDS'] = tlist[2]
    # Extract exon lengths from vulgar output, then add to target_dict
    # Each exon length is preceded by 'M '
    # Exon lengths are all in sequential list.
    vulgar_raw = tlist[3]
    vlist = re.findall(r'M\s\d*', vulgar_raw)
    vulgar = [int(x.replace('M ', '')) for x in vlist]
    target_dict['Vulgar'] = vulgar
    return target_dict


def splitter(CDS, vlist):
    # Cuts up the CDS into exon sequences and stores them in a list
    # This is done from the vulgar list of exon lengths
    # CDS is the sequence of concatenated exons
    counter = 0
    exonseq_list = []
    for v in vlist:
        exonseq_list.append(CDS[counter:counter + v])
        counter += v
    return exonseq_list


def writer(target_dict):
    # Writes the exon sequences into the outfile.fa
    # Each exon sequence is named by the order in the target (gene)
    outstring = ""
    vlist = target_dict['Vulgar']
    CDS = target_dict['CDS']
    exonseq_list = splitter(CDS, vlist)
    for exon_number, exonseq in enumerate(exonseq_list):
        #Only keep the exon if it is over 200bp
        if len(exonseq)>199:
            outstring += (">%s_%d\n%s\n") % (target_dict['Target'], exon_number,
                                         exonseq)
        else:
            pass
    return outstring


def main():
    # Engine for the program.
    args = sys.argv[1:]
    #args=['--in','/Users/josec/Desktop/NudiPreBait/NudiSilicoTest/test444out.txt'] #For testing
    #Print usage if no input is put in by user
    if not args:
        print usage
        sys.exit(1)
    if args[0] == '--in':
        infile = args[1]
    outfile_name = infile.replace('.txt', '_200.fa')
    outfile = open(outfile_name, 'w')
    target_list = openfile(infile)
    for target in target_list:
        # check if target file is blank
        if len(target) > 2:
            target_dict = parser(target)
            # Filter to remove low ID hits
            if float(target_dict['Percent']) < 98:
                continue
            outstring = writer(target_dict)
            outfile.write(outstring)
    outfile.close()


if __name__ == '__main__':
    main()
