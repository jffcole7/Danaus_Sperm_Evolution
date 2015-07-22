#!/bin/bash

#before running this script, please read readme.txt to find out what software and hardware are  required to run this script
#these are the genomic data and annotations necessary to run gKaKs for Danaus butterflies
#the following files may need to be unzipped for the program to work
#this contains the link to the cds for Danaus plexippus
wget "http://monarchbase.umassmed.edu/download/Dp_geneset_OGS2_cds.fasta.gz"
wget "http://monarchbase.umassmed.edu/download/Dp_geneset_OGS2.gff3.gz"
wget "http://monarchbase.umassmed.edu/download/Chry.ctg.v0.fa.gz"
wget "http://monarchbase.umassmed.edu/download/Gili.ctg.v0.fa.gz"

#this aligns chrysippus assembly to plexippus CDS, callculates dNdS ratios for all orthologous genes, using the YN method
#it may be tweaked however necessary for certain needs
perl gKaKs1.3_tweaked.pl \
-query_seq="Dp_geneset_OGS2_cds.fasta" \
-gff="Dp_geneset_OGS2.gff3" \
-hit_seq="Chry.ctg.v0.fa" \
-spe=6 \
-kaks_file="chrys-dnds.txt" \
-program=/usr/bin/KaKs_Calculator \
-KaKsCms=YN \

perl gKaKs1.3_tweaked.pl \
-query_seq="Dp_geneset_OGS2_cds.fasta" \
-gff="Dp_geneset_OGS2.gff3" \
-hit_seq="Gili.ctg.v0.fa.gz" \
-spe=6 \
-kaks_file="gili-dnds.txt" \
-program=/usr/bin/KaKs_Calculator \
-KaKsCms=YN \

#before taking kaks files into R using the R script, the header needs to be adjusted to be "R-friendly"
#this can be done either manually or using this simple set of scripts
perl -ni -e 'print unless $. ==1' gili-dnds.txt
perl -pi -e 'print "Sequence	Method	Ka	Ks	KaKs	P-value(fisher)	length	s-sites	n-sites	fold-sites	substitutions	s-substitions	n-substitutions	fold-s-substitutions	fold-n-substitutions	divergence-time	sub-rate-ratio	GC	ML-score	AIC	akaike-weight	model"' gili-dnds.txt

