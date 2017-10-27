######Set the working directory and number of threads######
wd='/Users/josec/Desktop/NudiPreBaitTest'
threads='7'


##STEP 1 - Generate the target sets

#Running Agalma to find target genes. 
#scripts for running agalma are here:
$wd/1_Generate_Target_Sets/AgalmaRuns/Nudis2AA/nudis2Commands.sh
$wd/1_Generate_Target_Sets/AgalmaRuns/AgalmaCDSfast5/fast5Commands.sh
#I put the output of these scripts in Output_1 so no need to rerun these two shell scripts

#Teasdale data from supp materials from Teasdale et al 2016
cd $wd/1_Generate_Target_Sets/Teasdale500/
#Make a list of geneids for Lotia gigantea. From table S1
#Lgig500_geneID.txt
#Convert to genbank accession numbers using this website: 
#https://biodbnet-abcc.ncifcrf.gov/db/db2dbRes.php
#List of accessions to download from genbank
#Lgig500_Accession.txt
#download fasta sequences from Accession list by running 
python /Users/josec/Desktop/NudiPreBaitTemplate/PreBait_git/fetcher.py > $wd/1_Generate_Target_Sets/Output_1/Teasdale500_Lgig_dna.fa


##STEP 2 - Find the Aplysia sequence for each target gene identified in Step 1

#Set up working and out directory for blast search. The blast query is the output from step1 so we will copy those files over
mkdir $wd/2_Blast2AplyCDS
mkdir $wd/2_Blast2AplyCDS/Output_2
mkdir $wd/2_Blast2AplyCDS/Blasting
cp $wd/1_Generate_Target_Sets/Output_1/* $wd/2_Blast2AplyCDS/Blasting
cd $wd/2_Blast2AplyCDS/Blasting

#Download CDS of Aplysia genome 3.0 from here:
#http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/    
unzip /Users/josec/Desktop/NudiPreBaitTest/2_Blast2AplyCDS/Aply3.0CDS/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna.zip
#Create blast database of this fasta file
makeblastdb -in $wd/2_Blast2AplyCDS/Aply3.0CDS/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna -parse_seqids -dbtype nucl
#Path to the Aplysia CDS database
aply_cds=$wd'/2_Blast2AplyCDS/Aply3.0CDS/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna'

#tblastx for DNA queries and tblastn for protein queries. max evalue is e-20
#List of Aplysia gene names (of the hits) printed to Nudis2_aply_blast.txt. Extraneous parts of seq name trimmed with sed command
tblastn -query Nudis2_Chr_wes_aa.fa -db $aply_cds -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads $threads | sed 's/.*cds_\(.*\..*\)_.*/\1/g' > Nudis2_aply_blast.txt
#Remove duplicate Aplysia Hits
sort Nudis2_aply_blast.txt | uniq > Nudis2_aply_blastu.txt


#Run again with tblastx for the other two DNA target files
tblastx -query Fast5_Chr_wes_dna.fa -db $aply_cds -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads $threads | sed 's/.*cds_\(.*\..*\)_.*/\1/g' > Fast5_aply_blast.txt
sort Fast5_aply_blast.txt | uniq > Fast5_aply_blastu.txt
tblastx -query Teasdale500_Lgig_dna.fa -db $aply_cds -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads $threads | sed 's/.*cds_\(.*\..*\)_.*/\1/g' > Teas_aply_blast.txt
sort Teas_aply_blast.txt | uniq > Teas_aply_blastu.txt
#Copy to output folder and rename file
for f in *u.txt ; do cp $f $wd/2_Blast2AplyCDS/Output_2/${f%aply_blastu.txt}genelist.txt;done


##STEP 3 - Create a shared database and count genes from each target set.

# Set up
mkdir $wd/3_Aplysia_Targets/Output_3
mkdir $wd/3_Aplysia_Targets/VennDiagram
cp $wd/2_Blast2AplyCDS/Output_2/* $wd/3_Aplysia_Targets/VennDiagram
cd $wd/3_Aplysia_Targets/VennDiagram
#Combine target sets and remove duplicates.
cat *.txt | sort | uniq > targets_aply.txt
# Add header line to each list
for f in *genelist.txt;do echo ${f%_genelist.txt} | cat - $f > temp && mv temp $f;done
#Execute R script to make venndiagram.pdf and move it to outfolder
cp /Users/josec/Desktop/NudiPreBaitTemplate/PreBait_git/VennDiagram.R $wd/3_Aplysia_Targets/VennDiagram/
Rscript VennDiagram.R 
cp VennDiagramBaits.pdf $wd/3_Aplysia_Targets/Output_3
#This pulls out the sequences from the Aplasia_CDS file and sends them to the output of Step_3. 
#while loop is a bit slow takes 20min.... maybe write python script for this step?
while read p ; do sed -n -e "/$p/,/>/ p" $aply_cds | sed ';$d' >>$wd/3_Aplysia_Targets/Output_3/targets.aply.fa ;done < targets_aply.txt
#clean up the outputs
sed -i '' "s/>.*cds_\(.\{2,20\}\..\).*/>\1/g" $wd/3_Aplysia_Targets/Output_3/targets.aply.fa


##STEP 4 - Map the introns onto the Aplysia Genome.

#Set up
cd /Users/josec/Desktop/NudiPreBaitTest/4_Map_Introns/
mkdir Output_4
cp /Users/josec/Desktop/NudiPreBaitTest/3_Aplysia_Targets/Output_3/targets.aply.fa /Users/josec/Desktop/NudiPreBaitTest/4_Map_Introns/
touch /Users/josec/Desktop/NudiPreBaitTest/4_Map_Introns/Output_4/targets_ex_out.txt
#Download Genome of Aplysia genome 3.0 from here:
#http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_genomic.fna.gz
unzip /Users/josec/Desktop/NudiPreBaitTest/4_Map_Introns/acl_ref_AplCal3.0_chrUn.fa.zip
#Run Exonerate. Run on multiple threads with GNU Parallel
#Takes a few hours. Put Output in the template directory so you dont have to run it. 
parallel --jobs $threads exonerate --model est2genome -q targets.aply.fa -t acl_ref_AplCal3.0_chrUn.fa -Q DNA -T DNA --showvulgar F --showalignment F --verbose 0 --ryo \"%qi\\t%pi\\t%qas\\t%V\\tEND\\n\" --fsmmemory 20G --bestn 1 --querychunktotal $threads --querychunkid  >> /Users/josec/Desktop/NudiPreBaitTest/4_Map_Introns/Output_4/targets_ex_out.txt ::: $(eval echo "{1..$threads}")
#The formating of the exonerate output with the --ryo flag is interpreted by VulgarityFilter.py


##STEP 5 - Cut the Aplysia target fasta into exons using custom python script
cp /Users/josec/Desktop/NudiPreBaitTest/4_Map_Introns/Output_4/targets_ex_out.txt /Users/josec/Desktop/NudiPreBaitTest/5_Aplysia_Exons
python $wd/PreBait_git/VulgarityFilter.py --in $wd/5_Aplysia_Exons/targets_ex_out.txt
