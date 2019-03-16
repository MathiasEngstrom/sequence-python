from sequence import Sequence
import matplotlib.pyplot as plt

# ---------------Initialising instances of Sequence to test the class methods--------------------- #

dna1 = Sequence('ATTGG')
dna1_copy = Sequence('ATTGG')
dna1_complement = Sequence('taacc')  # Using lower cases to confirm that the class is insensitive to case
dna2 = Sequence('AGTAT')
invalid_dna = Sequence('ATTGU')
rna1 = Sequence('AUUGG', 'RNA')
rna1_complement = Sequence('uaacc', 'rna')

# -------------------Asserting that the methods return the expected values------------------------ #

# Testing of n_bases()
assert dna1.n_bases() == 5

# Testing of __is_valid() and get_is_valid()
assert dna1.get_is_valid()
assert rna1.get_is_valid()
assert not invalid_dna.get_is_valid()

# Testing of __eq__()
assert dna1 == dna1_copy
assert dna1 != dna2

# Testing of complement()
assert dna1.complement() == dna1_complement and dna1.complement().get_type() == 'DNA'
assert rna1.complement('rna') == rna1_complement and rna1.complement('rna').get_type() == 'RNA'

# Testing of first_unmatched_basepair()
assert dna1.first_unmatched_basepair(dna1_copy) == -1
assert dna1.first_unmatched_basepair(dna2) == 1

# Testing of swap_mutation()
assert dna1.swap_mutation(dna2) == 3

# -----------------------------Working with genomes from input files------------------------------- #

# Reading genome_01.dat, printing number of bases and asserting that the Sequence object has correct format
genome1 = Sequence.sequence_from_file('genome_01.dat')
print('Number of bases in genome_01: ', genome1.n_bases())
assert genome1.get_is_valid

# Splitting genome1 into list of Sequence objects and calculating the length of the first gene
gene_list1 = genome1.gene_separation()
n_bases_gene1 = gene_list1[0].n_bases()
print('Number of bases in the first gene of genome_01: ', n_bases_gene1)

# Making a list of gene lengths from the list of genes gene_list1 and plotting it in a histogram
gene_lengths = []
for i in range(len(gene_list1)):
    gene_lengths.append(gene_list1[i].n_bases())
plt.hist(gene_lengths)
plt.xlabel('Gene length')
plt.ylabel('Number of genes')
plt.title('Histogram of gene lengths')
plt.savefig('script_project_hist.png')
plt.clf()

# Reading the genome_02.dat file and splitting into genes
genome2 = Sequence.sequence_from_file('genome_02.dat')
gene_list2 = genome2.gene_separation()

# Calculating number of swap mutations between each of the genes in the two genomes, and plotting them in a scatter plot
n_swaps = []
for i in range(len(gene_list1)):
    n_swaps.append(gene_list1[i].swap_mutation(gene_list2[i]))
plt.scatter(gene_lengths, n_swaps)
plt.xlabel('Gene length in number of base pairs')
plt.ylabel('Number of swap mutations')
plt.title('Scatter plot of gene lengths and number of swap mutations ')
plt.savefig('script_project_scatter.png')
