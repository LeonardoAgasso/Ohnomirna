#Substitution matrix for modified Needleman-Wunsch algorithm

from Bio import Align
from Bio.Align.substitution_matrices import Array

nonseed_nb = "ACGT"
seed_nb = "ZBRV"
all_nb = nonseed_nb+seed_nb # string: "ACGTZBRV"

counts = Array(nonseed_nb+seed_nb, dims=2) #define the empty matrix
print(counts)

#parameters to set
outseed_match = 1
outseed_mismatch = -1
#outseed_gap = 0

seed_match = 4
seed_mismatch = -4
#seed_gap = -2

# set the entrances of the subs matrix
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



print(counts)
