# PreBait

Walkthrough of bait design from Agalma output.

###########################################################

Clone the scripts into the PreBait folder with this:

git clone https://github.com/ignacio3437/PreBait.git

Make sure to set the number of threads and folder location in Start2Finish.sh

###########################################################

Directories starting with a number are the steps taken to get to the final baitset file that was sent to MycroArray for bait design and synthesis. Each of these directories has an outfolder which is the input for the next step.

Detailed descriptions:
1- Single copy orthologous protein coding phylogenetically informative genes were identified using the Agalma pipeline. The pipeline was run twice. Nudis2AA was run with the defaults (AA search) on all of the Nudibranch transcriptomes. CDSFast5 was run on a smaller subset using DNA data. The Teasdale 500 are the genes used in the Teasdale etal 2016 study. These genes were also targeted in this study.

2- To standardize across the different target sets and to find intron boundries in our transcriptome data, we used Aplysia californica as a reference genome. Each gene from each target was blasted to the Aply genome and if there was a hit '<'e-20, the Aply protein accession for that gene was saved in a list text file.

3- A venn diagram of the 3 target sets was created. From this point forward the target sets were combined to a single file. The fasta file of the Aply sequences was created.

4- Exonerate was used to identify the introns in the Aplysia target genes. Exonerate was run on multiple cores using GNU Parallel. The exonerate query sequence was the aply_targets and the exonerate target file was the Aplysia genome. The custom --ryo format was then parsed in the next step.

5- The exonerate output was parsed by VulgarityFilter.py. Takes the aply target genes and cuts them at the intron boundries into exons. Exons shorter than 200bp were filtered and thrown out. The sequences in the fasta file were named
'>'AplyProteinID_exonNumber
actagcagEXONSEQUENCEactgcatgc

6- Each exon of the Aplysia targets was used as a target sequence in exonerate. The queries were the Chromodoris transcriptomes. This pulled out the exons of the genes that we want to target in our exon capture experiment. The 'ryo' format was modified slightly and the output is a fasta file, not a text file.

7- At this point we would like to make exon alignments of the Chromodoris spp. to see if the exons are informative at the 'genus level'. The previous step does not pull out every exon for BOTH transcriptomes due to the fact that the query(Aplysia) is too distant from the transcriptome query. Additionally, exonerate sometimes pulls out several hits for each target per transcriptome. The BetterBest.py script picks the single best sequence for each exon target from all of the exonerate hits from both transcriptomes. This 'best' hit will then be blasted to the reciprocal transcriptome to then make exon alignments. A hit is assigned as 'best' based on the exonerate score and length. It removes hits that are longer than expected.

8- 
