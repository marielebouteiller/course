import re


def PeptideEncoding(dna, protein):
    peptides = []
    l_dna = len(dna)
    l_protein = 3*len(protein)
    rev = Reverse(dna)
    for i in range(l_dna-l_protein):
        if ProteinTranslation(Transcription(dna[i:i+l_protein])) == protein:
            peptides.append(dna[i:i+l_protein])
        if ProteinTranslation(Transcription(rev[i:i+l_protein])) == protein:
            pep = rev[i:i+l_protein]
            peptides.append(Reverse(pep))
    return peptides


def Transcription(dna):
    rna = ''
    for i in dna:
        if i == 'T':
            rna += 'U'
        else:
            rna += i
    return rna


def Reverse(dna):
    reverse_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev = ''
    for i in dna:
        rev = reverse_bases[i]+rev
    return rev


def ProteinTranslation(RNA):
    l = len(RNA)
    protein = ''
    for i in range(int(l/3)):
        codon = RNA[3*i:3*i+3]
        if RNA_table[codon] == '':
            break
        protein += RNA_table[codon]
    return protein


def PrintListOfPeptides(pep):
    s = ''
    for p in pep:
        s += p + '\n'
    return s[:-1]


RNA = open('/Users/dtj848/Documents/DOC/python/genome sequencing(BioinfII)/RNA_codon_table_1.txt').read().split('\n')
RNA_table = {}
for l in RNA:
    pattern = re.search(r'(\w{3,3})\s(.*)', l)
    RNA_table[pattern.group(1)] = pattern.group(2)


# dna1 = 'GAAACT'
# print(Transcription(dna1))
# print(ProteinTranslation(Transcription(dna1)))
# print(ProteinTranslation(Transcription(Reverse(dna1))))
# print(Reverse(dna1))

# dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
# protein = 'MA'
# print(PeptideEncoding(dna, protein))
# print(PrintListOfPeptides(PeptideEncoding(dna, protein)))

# input = open('/Users/dtj848/Downloads/dataset_96_7.txt').read().split('\n')
# dna_input = input[0]
# protein_input = input[1]
# f = open('data.txt', 'w')
# print(PrintListOfPeptides(PeptideEncoding(dna_input, protein_input)))
# f.write(PrintListOfPeptides(PeptideEncoding(dna_input, protein_input)))
# f.close()

Tyrocidine_B1 = 'VKLFPWFNQY'
#prot2 = 'QIQVLEG'
input = open(
    '/Users/dtj848/Downloads/Bacillus_brevis_genome.txt').read().split('\n')
genome = ''
for i in input:
    genome += i
# print(genome)
print(Tyrocidine_B1)
print(PrintListOfPeptides(PeptideEncoding(genome, Tyrocidine_B1)))
