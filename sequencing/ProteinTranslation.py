import re


def ProteinTranslation(RNA):
    l = len(RNA)
    protein = ''
    for i in range(int(l/3)):
        codon = RNA[3*i:3*i+3]
        if RNA_table[codon] == '':
            break
        protein += RNA_table[codon]
    return protein


RNA = open('/Users/dtj848/Documents/DOC/python/genome sequencing(BioinfII)/RNA_codon_table_1.txt').read().split('\n')
RNA_table = {}
for l in RNA:
    # i either have one letter in the second group or nothing, so * is for one element or nothing
    pattern = re.search(r'(\w{3,3})\s(.*)', l)
    RNA_table[pattern.group(1)] = pattern.group(2)

# input = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
# print(ProteinTranslation(input))

RNA_input = open('/Users/dtj848/Downloads/dataset_96_4.txt').read()
f = open('data.txt', 'w')
f.write(ProteinTranslation(RNA_input))
f.close()
