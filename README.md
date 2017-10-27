# PreBait
Selection of code snippets used to find exons from transcriptomes that match a set of divergent targets using exonerate. The intron boundaries are mapped and the transcriptome contigs are cut up as putative exons. These exons are then filtered and reciprocally searched for in each of the transcriptomes. Finally an exon alignment is generated for each target gene. This will then be used to design baits for exon capture.

Make sure to set the number of threads and folder location in Start2Finish.sh 