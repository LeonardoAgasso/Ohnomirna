#!/usr/bin/env python3
#
#All-vs-all Needleman-Wunsch alignment (+score) given a list of sequence in fasta format
#Output in a quasi-FASTA format
#
#Author: Leonardo Agasso, 2022
#
#idea from:
#https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-218


import sys
import errno
import warnings
import numpy as np


from optparse import OptionParser
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


#Seed parameters
seed_length = 6
final_notseed_bps = 1	#how many final nucleotides of the mature miRNA are not included in the seed

#Ordered nitrogenous bases to be substituted in seed
nonseed_nt = "ACGU"
seed_nt = "ZVRB"
all_nt = nonseed_nt+seed_nt


def ignore_broken_pipe(func):
	#avoid broken pipe error
	try:
		func()
	except IOError as e:
		if e.errno == errno.EPIPE:
			sys.exit(0)
		else:
			raise

def parse_args():
	
	parser = OptionParser(usage=format_usage('''
		%prog [OPTIONS] <FASTA >qFASTA

		Perform an all-vs-all alignment (using the Needleman-Wunsch algorithm) between the sequences in FASTA,
		returns a qFASTA file where the header contains the names of the two aligned sequences
		while the "body" contains the alignment score and the alignment itself.\033[0m
		Available seed types are (blue identifies the seed):

		  •		\033[7m'8mer'\033[0m 		            \033[7m'7mer-m8'\033[0m:
		  		    ..\033[34mNNNNNNN\033[0mN-5'		..\033[34mNNNNNNN\033[0mN-5'
				      |||||||			  |||||||
				ORF...\033[33mNNNNNNN\033[31mA\033[0m...	    ORF...\033[33mNNNNNNN\033[0mN...
				      87654321                    87654321
		 
		  •		\033[7m'7mer-A1'\033[0m                   \033[7m'6mer'\033[0m (default) :
		  		    ..N\033[34mNNNNNN\033[0mN-5'		..N\033[34mNNNNNN\033[0mN-5'
				       ||||||			   ||||||
				ORF...N\033[33mNNNNNN\033[31mA\033[0m...	    ORF...N\033[33mNNNNNN\033[0mN...
				      87654321                    87654321
	'''
	#parse command line arguments (colored using ANSI escape sequences)
	#Seed types that can be added in future
	#		• 'GUM' : to be defined
	#		• 'GUT' : to be defined
	#		• 'LP' : to be defined
	#		• 'BM' : to be defined
	#		• 'BT' : to be defined
	'''
		In-seed nucleotides are converted in the the following way:		
		  A ---> Z
		  C ---> V
		  G ---> R
		  U ---> B	
	'''	

		
	))

	parser.add_option('-t', '--seed-type', default='6mer', metavar='STRING',
						help='Specify the format of the seed you want to consider (default: %default)')

	parser.add_option('-g', '--gap-open', default=-5, metavar='INT',
						help='Specify the gap opening penalty (default: %default)')

	parser.add_option('-e', '--gap-extend', default=-1, metavar='INT',
						help='Specify the gap extension penalty (default: %default)')

	parser.add_option('-x', '--mismatch-score', default=-1, metavar='INT',
						help='Specify the mismatch score (default: %default)')
	
	parser.add_option('-m', '--match-score', default=1, metavar='INT',
						help='Specify the match score (default: %default)')

	parser.add_option('-y', '--seed-mismatch-score', default=-3, metavar='INT',
						help='Specify the seed mismatch score (default: %default)')

	parser.add_option('-n', '--seed-match-score', default=3, metavar='INT',
						help='Specify the seed match score (default: %default)')

	parser.add_option('-p', '--penalize-end-gaps', default=False, action='store_true',
						help='Penalize end gaps (default: %default)')

	parser.add_option('-o', '--one-alignment', default=False, action='store_true',
						help='Return only one alignment per pair (default: %default)')

	parser.add_option('-s', '--score-only', default=False, action='store_true',
						help='Return only the alignment score (default: %default)')

	parser.add_option('-r', '--remove-seed', default=False, action='store_true',
						help='Do not consider the seed weighted differently in the alignment (default: %default)')


	return parser.parse_args()

def format_usage(usage):
	
	def prefix_length(line):
		length = 0
		while length < len(line) and line[length] in (' ', '\t'):
			length += 1
		return length

	lines = usage.split('\n')
	while len(lines) and len(lines[0].strip()) == 0:
		del lines[0]
	while len(lines) and len(lines[-1].strip()) == 0:
		del lines[-1]

	plen = min(prefix_length(l) for l in lines if len(l.strip()) > 0)
	return '\n'.join(l[plen:] for l in lines)

def seed_modify(sequence):
	#convert the seed nucleotides of a miRNA to a new alphabet (AUCG--->ZVRB)
	l1 = ConvertToList(sequence)
	for i in myRange(len(sequence)-seed_length-final_notseed_bps, len(sequence)-final_notseed_bps-1, 1):
		for j in range(len(nonseed_nt)):
			if(l1[i] == nonseed_nt[j]):
				l1[i] = seed_nt[j]
	mod_sequence = ConvertToString(l1)
	return mod_sequence

def seed_restore(sequence):
	#restore the seed nucleotide to the standard alphabet (ZVRB--->AUCG) to display them
	l1 = ConvertToList(sequence)
	for i in myRange(len(sequence)-seed_length-final_notseed_bps, len(sequence)-final_notseed_bps-1, 1):
		for j in range(len(nonseed_nt)):
			if(l1[i] == seed_nt[j]):
				l1[i] = nonseed_nt[j]
	rest_sequence = ConvertToString(l1)
	return rest_sequence

def ConvertToList(string):
	#Convert a string to an ordered list of characters
	list1 = []
	list1[:0] = string
	return list1

def ConvertToString(list):
	#Convert a list to a string
	str1 = ""
	for elem in list:
		str1 += elem
	return str1

def myRange(start,end,step):
	#Define a range for the for-loop
	i = start
	while i < end:
		yield i
		i += step
		yield end

def subs_matrix(seed_flag):
	#Visualize to better handle the substitution matrix in an intuitive way (since the matrix is symmetric just the superior triangular matrix has to be filled)
	if seed_flag:
		subs_matrix = np.array([
		#	A	|		C	|		G	|		U	|
		[	3,			-1,			-1,			-1,		],	# A

		[	0,			3,			-1,			-1,		],	# C

		[	0,			0,			3,			-1,		],	# G

		[	0,			0,			0,			3,		],	# U
		])
	
	else:
		subs_matrix = np.array([
		#	A	|		C	|		G	|		U	|		Z	|		V	|		R	|		B
		[	3,			-1,			-1,			-1,			3,			-1,			-1,			-1	],	# A

		[	0,			3,			-1,			-1,			-1,			3,			-1,			-1	],	# C

    	[	0,			0,			3,			-1,			-1,			-1,			3,			-1	],	# G

    	[	0,			0,			0,			3,			-1,			-1,			-1,			3	],	# U

    	[	0,			0,			0,			0,			3,			-1,			-1,			-1	],	# Z

    	[	0,			0,			0,			0,			0,			3,			-1,			-1	],	# V

		[	0,			0,			0,			0,			0,			0,			3,			-1	],	# R

    	[	0,			0,			0,			0,			0,			0,			0,			3	]	# B
		])

	return subs_matrix

def Sim_matr_to_dic(string1,string2,seed_flag):
	#convert the substitution matrix into a dictionary
	sm = subs_matrix(seed_flag)
	list1 = ConvertToList(string1)

	if seed_flag:
		list3 = list1
		string3 = string1
	else:
		list2 = ConvertToList(string2)
		list3 = list1+list2
		string3 = string1+string2

	dic_subs_matrix = {}
	for i in range(0,len(string3)):
		for j in range(i,len(string3)):
			dic_subs_matrix[list3[i],list3[j]] = sm[i,j]

	return dic_subs_matrix

def print_fasta_header(seq, names, i, j, options):
	
	global seq_i, seq_j 
	
	if options.remove_seed:
		seq_i = seq[i]
		seq_j = seq[j]
		#print every alignment in a quasi FASTA format
		print(">{}: {}\t | \t{}: {}".format(names[i], seq_i, names[j], seq_j))
					
	else:
		seq_i = seed_modify(seq[i])
		seq_j = seed_modify(seq[j])
		#print every alignment in a quasi FASTA format
		print(">{}: {}\t | \t{}: {}".format(names[i], seed_restore(seq_i), names[j], seed_restore(seq_j)))

def normalized_score(score, options, seq_i, seq_j):
	#normalize the score between 0 and 1
	tot_len = len(seq_i)+len(seq_j)
	max_score = (seed_length*options.seed_match_score)+((0.5)*(tot_len)*options.match_score)
	score = score/max_score
	return score

def seed_definition(options):
	global seed_length
	global final_notseed_bps

	final_notseed_bps = 1

	#switch structure to define the seed length
	if options.seed_type=='8mer' or options.seed_type=='7mer-m8':
		seed_length = 7
	elif options.seed_type=='6mer' or options.seed_type=='7mer-A1':
		seed_length = 6

#function to print a matrix as a picture
def print_matrix(matrix):
	for i in range(0,len(matrix)):
		for j in range(0,len(matrix[i])):
			print(matrix[i][j], end="\t\t")
		print()

#function to save a matrix as a png. Darker colors represent higher values and axes values are given by a dictionary
def print_matrix_png(matrix, dic):
	import matplotlib.pyplot as plt

	fig = plt.figure()
	fig.set_size_inches(10.0, 10.0)
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.2, bottom=0.05)
	cax = ax.matshow(matrix, interpolation='nearest', cmap=plt.cm.Blues)
	fig.colorbar(cax)

	#Set the ticks to be at the edges of the matrix
	ax.set_xticks(np.arange(len(dic)), minor=False)
	ax.set_yticks(np.arange(len(dic)), minor=False)

	#Set the ticks labels to be the dictionary keys
	ax.set_xticklabels(dic.values(), minor=False)
	ax.set_yticklabels(dic.values(), minor=False)

	#Rotate the x axis labels
	plt.xticks(rotation=90)

	#Save the figure as a png
	fig.savefig('alignments_score_matrix.png')




def main():
	names = []
	sequences = []

	#Ignore the deprecationwarning from Bio.pairwise2
	warnings.filterwarnings("ignore", category=PendingDeprecationWarning)

	# Parse command line arguments
	options, args = parse_args()

	# Define seed properties
	seed_definition(options)

	# Define similarity matrix
	dic_sssmatr = Sim_matr_to_dic(nonseed_nt, seed_nt, options.remove_seed)

	#Define dictionary to store the names of the mirnas associated with the matrix indices
	dic_names = {}
	
	# Store ordered info about mirna IDs
	with open("/dev/stdin","r") as handle:
		for record in SeqIO.parse(handle, "fasta"):	#count how many sequences
			names.append(record.id)
			sequences.append(str(record.seq))
		
		#Define matrix to store the scores of every single alignment
		conf_matrix = np.zeros((len(names),len(names)))

	#for-loop that performs the all-vs-all confrontation
		for i in range(len(sequences)):
			for j in range(len(sequences)):
				if(i!=j):	#don't align a mirna with itself
					
					print_fasta_header(sequences, names, i, j, options)
					
					a = pairwise2.align.globalds(seq_i, 
												seq_j, 
												dic_sssmatr, 
												options.gap_open, 
												options.gap_extend, 
												one_alignment_only = options.one_alignment, 
												score_only = options.score_only, 
												penalize_end_gaps = options.penalize_end_gaps, 
												penalize_extend_when_opening = True)
					
					norm_score = round(normalized_score(a[0][2], options, seq_i, seq_j),3)
					print(format_alignment(*a[0]), ' Normalized score=' , norm_score)

					#fill the matrix with the scores
					conf_matrix[i,j] = a[0][2]
					#fill the dictionary with the names
					dic_names[i] = names[i]
					dic_names[j] = names[j]

	print_matrix_png(conf_matrix, dic_names)



if __name__=='__main__':

	#Ignore the deprecationwarning from Bio.pairwise2
	warnings.filterwarnings("ignore", category=PendingDeprecationWarning) 
	
	ignore_broken_pipe(main)
	