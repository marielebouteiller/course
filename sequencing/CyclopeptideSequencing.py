def CyclopeptideSequencing(spectrum):
    CandidatePeptides = [[]]
    FinalPeptides = []
    ParentMass = Mass_spectrum(spectrum)
    while bool(CandidatePeptides) != False:
        CandidatePeptides_totest = Expand(CandidatePeptides)
        CandidatePeptides = CandidatePeptides_totest.copy()
        for pep in CandidatePeptides_totest:
            if Mass_peptide(pep) == ParentMass:
                # a revoir le not in
                if CyclicSpectrum(pep) == spectrum and bool(pep in FinalPeptides) == False:
                    FinalPeptides.append(pep)
                CandidatePeptides.remove(pep)
            elif IsPeptideInSpectrum(pep, spectrum) == False:  # a revoir
                CandidatePeptides.remove(pep)
    return FinalPeptides


def CyclopeptideSequencingNoCopy(spectrum):
    CandidatePeptides = [[]]
    FinalPeptides = []
    ParentMass = Mass_spectrum(spectrum)

    while CandidatePeptides:
        CandidatePeptides = Expand(CandidatePeptides)
        valid_peptides = []

        for i, pep in enumerate(CandidatePeptides):
            if Mass_peptide(pep) == ParentMass \
                    and CyclicSpectrum(pep) == spectrum \
                    and pep not in FinalPeptides:
                FinalPeptides.append(pep)

            elif IsPeptideInSpectrum(pep, spectrum):
                valid_peptides.append(i)

        CandidatePeptides = [CandidatePeptides[k] for k in valid_peptides]

    return FinalPeptides


def IsPeptideInSpectrum(peptide, spectrum):
    # need to check the spectrum of the linear peptide and not the cyclopeptide as it is meant to be expanded so it is not linear
    l = LinearSpectrum(peptide)
    output = True
    for i in l:
        if bool(i in spectrum) == False:
            output = False
    return output


def CyclicSpectrum(peptide):
    PrefixMass = []
    PrefixMass.append(0)
    for i in range(len(peptide)):
        PrefixMass.append(PrefixMass[-1] + int(peptide[i]))
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


def LinearSpectrum(peptide):
    PrefixMass = []
    PrefixMass.append(0)
    for i in range(len(peptide)):
        PrefixMass.append(PrefixMass[-1] + int(peptide[i]))
    LinearSpectrum = []
    LinearSpectrum.append(0)
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum


def Mass_spectrum(spectrum):
    mass = max(spectrum)
    return(mass)


def Mass_peptide(peptide):
    mass = 0
    for i in peptide:
        mass += i
    return mass


def Expand(list_of_peptides):
    new_list_of_peptides = []
    for l in list_of_peptides:
        for a in AminoAcidMass:
            l.append(a)
            new_list_of_peptides.append(l)
            l = l[:-1]
    return new_list_of_peptides


def PrintFinalPeptides(peptides):
    output = ''
    for pep in peptides:
        for p in pep:
            output += str(p)+'-'
        output = output[:-1]
        output += ' '
    return output[:-1]


AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113,
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

# input = []
# input2 = [0, 113, 128, 186, 241, 299, 314, 427]
# # spectrum = [0, 71, 87, 101, 158, 188, 259, 444]
# # print(Mass(input2))
# # print(IsPeptideInSpectrum(input2, spectrum))
# # print(CyclicSpectrum(input2))
# print(CyclopeptideSequencingNoCopy(input2))
# print(CyclopeptideSequencing(input2))
# print(Mass_peptide(input2))
# sol = [[113, 128, 186], [113, 186, 128], [128, 113, 186],
#      [128, 186, 113], [186, 113, 128], [186, 128, 113]]
# print(PrintFinalPeptides(CyclopeptideSequencing(input2)))

# input = open('/Users/dtj848/dataset_100_6.txt').read()
# input = input[:-1]
# input2 = input.split()
# input3 = []
# for i in input2:
#     input3.append(int(i))

# print(input3)


def parser(filename):
    input = open(filename).read()[:-1]
    return [int(i) for i in input.split(' ')]


input = parser('/Users/dtj848/dataset_100_6.txt')


# je dois avoir mon input as a list of integer, not a giant string
# f = open('data.txt', 'w')
# f.write(PrintFinalPeptides(CyclopeptideSequencing(input3)))
# f.close()
# print(CyclopeptideSequencing(input3))
