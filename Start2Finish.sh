######Set the working directory and number of threads######
WD='/Users/josec/Desktop/NudiPreBait2'
threads='7'


##STEP 1 - Generate the target sets

#Running Agalma to find target genes.
#scripts for running agalma are here:
$WD/1_Generate_Target_Sets/AgalmaRuns/Nudis2AA/nudis2Commands.sh
$WD/1_Generate_Target_Sets/AgalmaRuns/AgalmaCDSfast5/fast5Commands.sh
#I put the output of these scripts in Output_1 so no need to rerun these two shell scripts

#Teasdale data from supp materials from Teasdale et al 2016
cd $WD/1_Generate_Target_Sets/Teasdale500/
#Make a list of geneids for Lotia gigantea. From table S1
#Lgig500_geneID.txt
#Convert to genbank accession numbers using this website:
#https://biodbnet-abcc.ncifcrf.gov/db/db2dbRes.php
#List of accessions to download from genbank
#Lgig500_Accession.txt
#download fasta sequences from Accession list by running
python $WD/PreBait/fetcher.py > $WD/1_Generate_Target_Sets/Output_1/Teasdale500_Lgig_dna.fa


##STEP 2 - Find the Aplysia sequence for each target gene identified in Step 1

#Set up working and out directory for blast search. The blast query is the output from step1 so we will copy those files over
mkdir $WD/2_Blast2AplyCDS/Output_2
mkdir $WD/2_Blast2AplyCDS/Blasting
cp $WD/1_Generate_Target_Sets/Output_1/* $WD/2_Blast2AplyCDS/Blasting
cd $WD/2_Blast2AplyCDS/Aply3.0CDS

#Download CDS of Aplysia genome 3.0 from here:
#curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna.gz
#Create blast database of this fasta file
makeblastdb -in $WD/2_Blast2AplyCDS/Aply3.0CDS/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna -parse_seqids -dbtype nucl
#Path to the Aplysia CDS database
aply_cds=$WD'/2_Blast2AplyCDS/Aply3.0CDS/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna'

#tblastx for DNA queries and tblastn for protein queries. max evalue is e-20
#List of Aplysia gene names (of the hits) printed to Nudis2_aply_blast.txt. Extraneous parts of seq name trimmed with sed command
cd $WD/2_Blast2AplyCDS/Blasting
tblastn -query Nudis2_Chr_wes_aa.fa -db $aply_cds -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads $threads | sed 's/.*cds_\(.*\..*\)_.*/\1/g' > Nudis2_aply_blast.txt
#Remove duplicate Aplysia Hits
sort Nudis2_aply_blast.txt | uniq > Nudis2_aply_blastu.txt


#Run again with tblastx for the other two DNA target files
tblastx -query Fast5_Chr_wes_dna.fa -db $aply_cds -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads $threads | sed 's/.*cds_\(.*\..*\)_.*/\1/g' > Fast5_aply_blast.txt
sort Fast5_aply_blast.txt | uniq > Fast5_aply_blastu.txt
tblastx -query Teasdale500_Lgig_dna.fa -db $aply_cds -evalue 1e-20 -outfmt '6 sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads $threads | sed 's/.*cds_\(.*\..*\)_.*/\1/g' > Teas_aply_blast.txt
sort Teas_aply_blast.txt | uniq > Teas_aply_blastu.txt
#Copy to output folder and rename file
for f in *u.txt ; do cp $f $WD/2_Blast2AplyCDS/Output_2/${f%aply_blastu.txt}genelist.txt;done


##STEP 3 - Create a shared database and count genes from each target set.

# Set up
mkdir $WD/3_Aplysia_Targets/Output_3
mkdir $WD/3_Aplysia_Targets/VennDiagram
cp $WD/2_Blast2AplyCDS/Output_2/* $WD/3_Aplysia_Targets/VennDiagram
cd $WD/3_Aplysia_Targets/VennDiagram
#Combine target sets and remove duplicates.
cat *.txt | sort | uniq > targets_aply.txt
# Add header line to each list
for f in *genelist.txt;do echo ${f%_genelist.txt} | cat - $f > temp && mv temp $f;done
#Execute R script to make venndiagram.pdf and move it to outfolder
cp $WD/PreBait/VennDiagram.R $WD/3_Aplysia_Targets/VennDiagram/
Rscript VennDiagram.R
cp VennDiagramBaits.pdf $WD/3_Aplysia_Targets/Output_3
#This pulls out the sequences from the Aplasia_CDS file and sends them to the output of Step_3.
#while loop is a bit slow takes 20min.... maybe write python script for this step?
#I put this file in the google folder so no need to execute....
while read p ; do sed -n -e "/$p/,/>/ p" $aply_cds | sed ';$d' >>$WD/3_Aplysia_Targets/Output_3/targets.aply.fa ;done < targets_aply.txt
#clean up the outputs
sed -i '' "s/>.*cds_\(.\{2,20\}\..\).*/>\1/g" $WD/3_Aplysia_Targets/Output_3/targets.aply.fa


##STEP 4 - Map the introns onto the Aplysia Genome.

#Set up
cd $WD/4_Map_Introns/
mkdir Output_4
cp $WD/3_Aplysia_Targets/Output_3/targets.aply.fa $WD/4_Map_Introns/
touch $WD/4_Map_Introns/Output_4/targets_ex_out.txt
#Download Genome of Aplysia genome 3.0 from here:
#curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_genomic.fna.gz
unzip $WD/4_Map_Introns/acl_ref_AplCal3.0_chrUn.fa.zip

#http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_genomic.fna.gz
unzip $WD/4_Map_Introns/acl_ref_AplCal3.0_chrUn.fa.zip
#Run Exonerate. Run on multiple threads with GNU Parallel
#Takes quite a few hours. More than 3 on my computer.
#I put a copy of this file in the Step 5 directory so you dont have to run it.
parallel --jobs $threads exonerate --model est2genome -q targets.aply.fa -t acl_ref_AplCal3.0_chrUn.fa -Q DNA -T DNA --showvulgar F --showalignment F --softmasktarget T --verbose 0 --ryo \"%qi\\t%pi\\t%qas\\t%V\\tEND\\n\" --fsmmemory 20G --bestn 1 --querychunktotal $threads --querychunkid  >> $WD/4_Map_Introns/Output_4/targets_ex_out.txt ::: $(eval echo "{1..$threads}")
#The formating of the exonerate output with the --ryo flag is interpreted by VulgarityFilter.py


##STEP 5 - Cut the Aplysia target fasta into exons using custom python script

cp $WD/4_Map_Introns/Output_4/targets_ex_out.txt $WD/5_Aplysia_Exons
cd $WD/5_Aplysia_Exons
python $WD/PreBait/VulgarityFilter.py --in $WD/5_Aplysia_Exons/targets_ex_out.txt
#Generate file for HybPiper Reference
python $WD/Prebait/VulgarityFilter2.py --in /Users/josec/Desktop/NudiPreBait2/5_Aplysia_Exons/targets_ex_out.txt


##STEP 6 - Cut the Nudi transcriptomes to match the aplysia exons
#Old step 6 Clean transcriptomes we are skipping because it does not make a difference. In other words, the transcriptome cleaning did not remove the "best hit" so the betterbest.py will still pick the same sequence, regardless of if there are other spurious sequences in the transcriptomes.

#Set up
cp $WD/1_Generate_Target_Sets/Transcriptomes/*.fasta $WD/6_Nudi_Exons
#Renaming file for clarity
cp $WD/5_Aplysia_Exons/targets_ex_out_200.fa $WD/6_Nudi_Exons/Aply_exons_200.fa
cd $WD/6_Nudi_Exons
#Run exonerate, output is a fasta file!
for txtm in Chr*.fasta; do parallel --jobs $threads exonerate --model est2genome -q Aply_exons_200.fa -t $txtm -Q DNA -T DNA --querychunktotal $threads --showvulgar F --showalignment F --verbose 0 --fsmmemory 20G --bestn 1 --ryo '\>%qi:TT%ti:HH%qab/qae:SS%s:PP%ps:LL%tal/%ql\\n%tas\\n' --querychunkid >> exons_$txtm ::: $(eval echo "{1..$threads}"); done


##STEP 7- Pick the 'best' hit for each exon from the (2) exonerate runs.

#Set up
for exon_file in exons*; do cp $exon_file $WD/7_Best_Exons/${exon_file%_trinity.fasta}.fasta;done
cd $WD/7_Best_Exons
#Run custom python script to save best hit.
python $WD/PreBait/BetterBest.py --exdir $PWD/


##STEP 8- Reciprical Blast

#Set up
cp $WD/7_Best_Exons/*best.fasta $WD/8_Recip_Blast/
cd $WD/8_Recip_Blast/
mkdir $WD/8_Recip_Blast/blastdb
cp $WD/1_Generate_Target_Sets/Transcriptomes/*.fasta $WD/8_Recip_Blast/blastdb

#Make a combined file of all the Chr exons names to pull out the ones that do not have a reciprical blast hit
cat Chr_mag_best.fasta Chr_wes_best.fasta | grep ">" | sed 's/:.*//' | sed 's/\>//' > Chr_exon_names.txt
#Make one line fasta for both best fasta to pull out 'singletons' later
cat Chr_mag_best.fasta Chr_wes_best.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | sed 's/:HH.*//' > Chr_best_combined.fa
#Make Blast database of each transcriptome
for txtm in $WD/8_Recip_Blast/blastdb/*.fasta; do makeblastdb -in $txtm -parse_seqids -dbtype nucl;done
#Reciprical Megablast and return top blast hit csv format
blastn -task dc-megablast -query Chr_mag_best.fasta -db blastdb/Chr_wes_trinity.fasta -evalue 1e-20 -outfmt '10 qseqid qseq sseqid sseq pident' -max_target_seqs 1 -max_hsps 1 -num_threads $threads >> Chr_Chr_blast.txt
blastn -task dc-megablast -query Chr_wes_best.fasta -db blastdb/Chr_mag_trinity.fasta -evalue 1e-20 -outfmt '10 qseqid qseq sseqid sseq pident' -max_target_seqs 1 -max_hsps 1 -num_threads $threads >> Chr_Chr_blast.txt


##STEP 9- Exon Filtering

#Set Up
cp $WD/8_Recip_Blast/Chr_Chr_blast.txt $WD/9_Exon_Filtering
cp $WD/8_Recip_Blast/Chr_best_combined.fa $WD/9_Exon_Filtering
cp $WD/8_Recip_Blast/Chr_exon_names.txt $WD/9_Exon_Filtering
cd $WD/9_Exon_Filtering

# Parse the output and generate alignments with RBlast_Parser.py
python $WD/PreBait/RBlast_Parser.py --btxt Chr_Chr_blast.txt > stats.txt
#Pull out singletons (ie seqs that did not have a hit in the reciprical blast)
while read p; do grep -A 1 $p: Chr_best_combined.fa; done <singletons.txt | sed 's/:.*//' > singletons.fa
#Combine singletons with Filtered blast results.
cat singletons.fa Fil_Bhits.fasta > Ktar.fa


##STEP 10- Ktar Final outfiles
cd $WD/10_Ktar
cp $WD/9_Exon_Filtering/stats.txt $PWD
cp $WD/9_Exon_Filtering/Ktar.fa $PWD

#Make file of just the names of the exons targeted. Exons are named after the aplysia protein they translate into.
cat Ktar.fa | grep ">" | sed 's/\>//' > Exon_names.txt
#Make a file listing the genes targeted
cat Exon_names.txt | sed 's/\..*//' | sort | uniq > Gene_names.txt
#Add to stats file number of genes targeted
wc -l Gene_names.txt | sed 's/ Gene_names.txt//' | sed 's/    /Number of genes kept: /'>> stats.txt


#All Done!
