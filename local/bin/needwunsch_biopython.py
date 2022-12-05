#All-vs-all Needleman-Wunsch alignment (+score) given a list of sequence in fasta format

from Bio import SeqIO
from Bio import pairwise2
from Bio import Align
from Bio.Align import PairwiseAligner	#consider using this one instead of deprecated pairwise2

from Bio.pairwise2 import format_alignment


seed_length = 6
final_notseed_bps = 1


def seed_modify(sequence):
	l1 = ConvertToList(sequence)
	for i in myRange(len(sequence)-seed_length-final_notseed_bps, len(sequence)-final_notseed_bps-1, 1):
		match l1[i]:
			case "A":
				l1[i] = "Z"
			case "U":
				l1[i] = "B"
			case "G":
			 	l1[i] = "R"
			case "C":
			 	l1[i] = "V"
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


#_____________MAIN_____________
names = []
sequences = []

# Store ordered info about mirna ids
with open("mockmirs_prova.fa") as handle:
	for record in SeqIO.parse(handle, "fasta"):	#count how many sequences
		names.append(record.id)
		sequences.append(str(record.seq))

for i in range(len(sequences)):
	for j in range(len(sequences)):
		if(i!=j):	#don't align a mirna with itself
			alignments = pairwise2.align.globalxx(seed_modify(sequences[i]), seed_modify(sequences[j]),  one_alignment_only=True)	#list of alignments, notice we are using global alignments
			print(">{}: {}\t vs. \t{}: {}".format(names[i],sequences[i],names[j],sequences[j]))
			for a in alignments:
				print(format_alignment(*a))
