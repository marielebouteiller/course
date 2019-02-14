import re


def LinearScoring(peptide, spectrum):
    theo_spectrum = LinearSpectrum(peptide)
    count = 0
    for i in spectrum:
        if i in theo_spectrum:
            count += 1
            theo_spectrum.remove(i)
    return count


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

# peptide = 'NQEL'
# spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
# print(LinearSpectrum(peptide))
# print(LinearScoring(peptide, spectrum))


def parser(filename):
    input = open(filename).read()[:-1]
    input.split('\n')
    peptide = input.split('\n')[0]
    spectrum = [int(i) for i in input.split('\n')[1].split(' ')]
    return [peptide, spectrum]


peptide = parser('/Users/dtj848/dataset_4913_1.txt')[0]
spectrum = parser('/Users/dtj848/dataset_4913_1.txt')[1]

print(LinearScoring(peptide, spectrum))
