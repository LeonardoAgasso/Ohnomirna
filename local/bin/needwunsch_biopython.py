#!/usr/bin/env python3

#All-vs-all Needleman-Wunsch alignment (+score) given a list of sequence in fasta format

#idea coming from:
#https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-218

#The Bio.Align module contains the PairwiseAligner class for global and local alignments
#using the Needleman-Wunsch, Smith-Waterman, Gotoh (three-state), and Waterman-Smith-Beyer
#global and local pairwise alignment algorithms, with numerous options to change the alignment parameters.

from Bio import SeqIO
from Bio import Align
from Bio.Align.substitution_matrices import Array

seed_length = 6
final_notseed_bps = 1
#ORDERED nitrogenous bases to be substituted in seed
nonseed_nb = "ACGU"
seed_nb = "ZVRB"
all_nb = nonseed_nb+seed_nb

def seed_modify(sequence):
	l1 = ConvertToList(sequence)
	for i in myRange(len(sequence)-seed_length-final_notseed_bps, len(sequence)-final_notseed_bps-1, 1):
		for j in range(len(nonseed_nb)):
			if(l1[i] == nonseed_nb[j]):
				l1[i] = seed_nb[j]
	mod_sequence = ConvertToString(l1)
	return mod_sequence

def seed_restore(sequence):
	l1 = ConvertToList(sequence)
	for i in myRange(len(sequence)-seed_length-final_notseed_bps, len(sequence)-final_notseed_bps-1, 1):
		match l1[i]:
			case "Z":
				l1[i] = "A"
			case "B":
				l1[i] = "U"
			case "R":
				l1[i] = "G"
			case "V":
				l1[i] = "C"
	rest_sequence = ConvertToString(l1)
	return rest_sequence

def ConvertToList(string):
	list1 = []
	list1[:0] = string
	return list1

def ConvertToString(list):
	str1 = ""
	for elem in list:
		str1 += elem
	return str1

def myRange(start,end,step):
	i = start
	while i < end:
		yield i
		i += step
		yield end

def Seed_Sensistive_SimMatr():
	counts = Array(nonseed_nb+seed_nb, dims=2) #define the empty 8x8 substitution-matrix
	#print(counts)

	#parameters to set
	outseed_match = 1
	outseed_mismatch = -1
	#outseed_gap = 0

	seed_match = 3
	seed_mismatch = -3
	#seed_gap = -2

	# set the entrances of the subs matrix
	for i in range(len(all_nb)):
		for j in range(len(all_nb)):
			if(all_nb[i] in seed_nb and all_nb[j] in seed_nb): #only when confronting two nbs both in the seed
				if(all_nb[i]==all_nb[j]):
					counts[i][j]=seed_match
				else:
					counts[i][j]=seed_mismatch
			else:
				if(all_nb[i]==all_nb[j]):
					counts[i][j]=outseed_match
				else:
					counts[i][j]=outseed_mismatch
	print(counts)
	return counts

#_____________MAIN_____________
names = []
sequences = []

seq = "ACGTACGTACGTACGTACGT"
print(seq)

seq_mod = seed_modify(seq)
print(seq_mod)

aligner = Align.PairwiseAligner()
aligner.mode = "global" # can be changed to "local", default is global
aligner.gap_score = -5

aligner.substitution_matrix = Seed_Sensistive_SimMatr()   #p.97 Biopython tutorial


# Store ordered info about mirna ids
with open("mockmirs_prova.fa") as handle:
	for record in SeqIO.parse(handle, "fasta"):	#count how many sequences
		names.append(record.id)
		sequences.append(str(record.seq))

print(sequences[1], sequences[2])

for i in range(len(sequences)):
	for j in range(len(sequences)):
		if(i!=j):	#don't align a mirna with itself
			seq_mod_i=seed_modify(sequences[i])
			seq_mod_j=seed_modify(sequences[j])

			alignments = aligner.align(seq_mod_i,seq_mod_j)	#list of alignments, notice we are using global alignments
			print(">{}: {}\t vs. \t{}: {}".format(names[i],seq_mod_i,names[j],seq_mod_j))
			for alignment in alignments[0]:
				print(alignment)

			score = aligner.score(seq_mod_i, seq_mod_j)
			print("Score: ", score, "\n")
