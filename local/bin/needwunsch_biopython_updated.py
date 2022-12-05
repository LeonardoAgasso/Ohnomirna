#All-vs-all Needleman-Wunsch alignment (+score) given a list of sequence in fasta format

from Bio import SeqIO
from Bio import Align
from Bio.Align.substitution_matrices import Array

#seed portion in a mature miRNA (to be better defined)
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

def ScoreMatrix():
	matr = Array(nonseed_nb+seed_nb, dims=2) #define the empty 8x8 score-matrix
	print(matr)

	#parameters to set
	outseed_match = 1
	outseed_mismatch = -1
	#outseed_gap = 0     gap are not a parametere in score matrices

	seed_match = 4
	seed_mismatch = -4
	#seed_gap = -2		gap are not a parametere in score matrices

	# set the entrances of the score matrix
	for i in range(len(all_nb)):
		for j in range(len(all_nb)):
			if(all_nb[i] in seed_nb and all_nb[j] in seed_nb): #only when confronting two nitrogenous bases both in the seed
				if(all_nb[i]==all_nb[j]):
					counts[i][j]=seed_match
				else:
					counts[i][j]=seed_mismatch
			else:
				if(all_nb[i]==all_nb[j]):
					counts[i][j]=outseed_match
				else:
					counts[i][j]=outseed_mismatch

	return matr

def test_seedModifier(testseq):
	if(testseq==""):
		seq = "ACGTACGTACGTACGTACGT" #sequenza di prova per testare la sostituzione delle basi nel seed
		print("using predefined testseq:")
	print(seq)
	seq_mod = seed_modify(seq)
	print(seq_mod)



#_____________MAIN_____________
names = []
sequences = []

aligner = Align.PairwiseAligner()
aligner.mode = "global" # can be changed to "local"


# Store ordered info about mirna ids
with open("mockmirs_prova.fa") as handle:
	for record in SeqIO.parse(handle, "fasta"):	#count how many sequences
		names.append(record.id)
		sequences.append(str(record.seq))

for i in range(len(sequences)):
	for j in range(len(sequences)):
		if(i!=j):	#don't align a mirna with itself
			alignments = aligner.align(sequences[i], sequences[j])	#list of alignments, notice we are using global alignments
			print(">{}: {}\t vs. \t{}: {}".format(names[i],sequences[i],names[j],sequences[j]))
			for alignment in alignments:
				print(alignment)

			score = aligner.score(sequences[i], sequences[j])
			print("Score: ", score, "\n")
