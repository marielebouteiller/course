import re


def LinearSpectrum(peptide):
    PrefixMass = []
    PrefixMass.append(0)
    for i in range(len(peptide)):
        PrefixMass.append(PrefixMass[-1] + int(AminoAcidMass[peptide[i]]))
    LinearSpectrum = []
    LinearSpectrum.append(0)
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum


f = open('/Users/dtj848/Documents/DOC/python/genome sequencing(BioinfII)/integer_mass_table.txt').read().split('\n')
AminoAcidMass = {}
for l in f:
    pattern = re.search(r'(\w)\s(\d*)', l)
    AminoAcidMass[pattern.group(1)] = pattern.group(2)
# print(AminoAcidMass)

# peptide = 'NQEL'
# print(LinearSpectrum(peptide))


peptide3 = 'IAKPAKLLLWHVYFVDQFDEAGGFACDIQNPQYKI'
peptide2 = 'DHTYQFTHSNWIIAWLPLFEECQSVIQFQHIEQNLCDCCYKLYA'

output = LinearSpectrum(peptide3)
s = ''
for i in output:
    s += str(i) + ' '

f_output = open('data.txt', 'w')
f_output.write(s[:-1])
f_output.close()
