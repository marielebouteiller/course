import re


def CyclicSpectrum(peptide):
    PrefixMass = []
    PrefixMass.append(0)
    for i in range(len(peptide)):
        PrefixMass.append(PrefixMass[-1] + int(AminoAcidMass[peptide[i]]))
    CycloSpectrum = []
    CycloSpectrum.append(0)
    peptideMass = PrefixMass[-1]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            CycloSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(peptide):
                overlap = peptideMass-(PrefixMass[j]-PrefixMass[i])
                CycloSpectrum.append(overlap)
    CycloSpectrum.sort()
    return CycloSpectrum


f = open('/Users/dtj848/Documents/DOC/python/genome sequencing(BioinfII)/integer_mass_table.txt').read().split('\n')
AminoAcidMass = {}
for l in f:
    pattern = re.search(r'(\w)\s(\d*)', l)
    AminoAcidMass[pattern.group(1)] = pattern.group(2)
print(AminoAcidMass)
peptide = 'LEQN'
print(CyclicSpectrum(peptide))


peptide3 = 'YALGDNLIFGSACY'
#peptide2 = 'DHTYQFTHSNWIIAWLPLFEECQSVIQFQHIEQNLCDCCYKLYA'

# output = CyclicSpectrum(peptide3)
# s = ''
# for i in output:
#     s += str(i) + ' '
# # print(s[:-1])

# f_output = open('data.txt', 'w')
# f_output.write(s[:-1])
# f_output.close()
