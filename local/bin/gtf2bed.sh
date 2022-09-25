#!/bin/bash

#***************************************************************************************************************#
#                                                                                                               #
#			Naif converter from an input GTF to an output in a slightly modified BED format structured as follow:     #
#				$1: chromosome (e.g. "chr1", "chr5", "chrX", ...)                                                       #
#				$2: starting position of the feature in the chromosome (the chromosme start at position 0)              #
#				$3: ending position of the feature in the chromosome                                                    #
#				$4: name: the gene_id in the GTF file is here considered as the name of the feature                     #
#				$5: score is a field present both in GTF and BED files                                                  #
#				$6: strand (can be "+", "-" or "." if there's no strand)                                                #
#				$7: the entire "attribute" field of the GTF file (column 9 of the GTF file) excluding the gene_id			  #
#                                                                                                               #     
#			The converter also provides te uhmbern of rows in both input and output files                             #
#                                                                                                               #
#***************************************************************************************************************#

echo " "
echo Give me the gtf \file \(include the .gtf at the end\)\:
read VARGTF
echo " "
echo Number of rows \in GTF \file\:
wc -l $VARGTF | awk '{print $1}'
echo " "
echo Give me the name of the .bed \file you\'re creating \(include .bed at the end\)\:
read VARBED
echo " " 

awk -F"\t" '{split($9,a," "); print $1 "\t" $4 "\t" $5 "\t" substr(a[2], 2, 15) "\t" $6 "\t" $7 "\t" substr($9, length(a[2])+10, length($9))}' $VARGTF > $VARBED

echo Here\'s a preview of your BED \file:
head $VARBED
echo " "
echo Number of rows \in BED \file\:
wc -l $VARBED | awk '{print $1}'
echo " "

