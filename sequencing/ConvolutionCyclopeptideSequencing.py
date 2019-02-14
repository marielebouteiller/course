import os


def ConvolutionCyclopeptideSequencing(M, N, Spectrum):
    Leaderboard = [[]]
    LeaderPeptide = []
    ParentMass = Mass_spectrum(spectrum)
    while Leaderboard:
        Leaderboard_totest = Expand(Leaderboard, Spectrum, M)
        Leaderboard = Leaderboard_totest.copy()
        if not LeaderPeptide:
            score_leader_peptide = 0
        else:
            score_leader_peptide = CyclopeptideScoring(
                LeaderPeptide[0], spectrum)
        for peptide in Leaderboard_totest:
            if Mass_peptide(peptide) == ParentMass and CyclopeptideScoring(peptide, spectrum) > score_leader_peptide:
                LeaderPeptide = []
                LeaderPeptide.append(peptide.copy())
                score_leader_peptide = CyclopeptideScoring(
                    peptide.copy(), spectrum)
            elif Mass_peptide(peptide) == ParentMass and CyclopeptideScoring(peptide, spectrum) == score_leader_peptide:
                LeaderPeptide.append(peptide.copy())
            elif Mass_peptide(peptide) > ParentMass:
                Leaderboard.remove(peptide)
        if not Leaderboard:  # if Laederboard is empty at that point, I just skip the next Trim obviously
            break
        Leaderboard = Trim(Leaderboard, spectrum, N)
        # if len(Leaderboard[0]) == 9:
        #     print(Leaderboard)
    return LeaderPeptide


def Expand(list_of_peptides, spectrum, M):
    new_list_of_peptides = []
    ext_alphabet = find_ext_alphabet(Convolution(spectrum), M)
    for l in list_of_peptides:
        for a in ext_alphabet:
            l.append(a)
            new_list_of_peptides.append(l)
            l = l[:-1]
    return new_list_of_peptides


def Trim(Leaderboard, spectrum, N):
    Scores = []
    if len(Leaderboard) < N:
        Leaderboard_to_return = Leaderboard
    else:
        for pep in Leaderboard:
            Scores.append(LinearScoring(pep, spectrum))
        Scores.sort(reverse=True)
        Score_min = Scores[N-1]
        valid_pep = []
        for i, pep in enumerate(Leaderboard):
            if LinearScoring(pep, spectrum) >= Score_min:
                valid_pep.append(i)
        Leaderboard_to_return = []
        for k in valid_pep:
            Leaderboard_to_return.append(Leaderboard[k])
    # Leaderboard = [Leaderboard[i] for i in valid_pep]
    return Leaderboard_to_return


def Mass_spectrum(spectrum):
    mass = max(spectrum)
    return(mass)


def Mass_peptide(peptide):
    mass = 0
    for i in peptide:
        mass += i
    return mass


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
        PrefixMass.append(PrefixMass[-1] + int(peptide[i]))
    LinearSpectrum = []
    LinearSpectrum.append(0)
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum


def Convolution(spectrum):
    output = []
    for mass in spectrum:
        for mass_prime in spectrum:
            if mass_prime > mass:
                output.append(mass_prime-mass)
    # output.sort()
    return output


def find_ext_alphabet(conv, M):
    conv_restricted = []
    for i in conv:
        if i < 201 and i > 56:
            conv_restricted.append(i)
    conv_count = {}
    for i in conv_restricted:
        if i not in conv_count.keys():
            conv_count[i] = 0
        conv_count[i] += 1
    count = list(conv_count.values())
    count.sort(reverse=True)
    threshold = count[M]
    ext_alphabet = []
    for i in conv_restricted:
        if conv_count[i] >= threshold:
            if i not in ext_alphabet:
                ext_alphabet.append(i)
    ext_alphabet.sort()
    return ext_alphabet


def CyclopeptideScoring(peptide, spectrum):
    theo_spectrum = CyclicSpectrum(peptide)
    count = 0
    for i in spectrum:
        if i in theo_spectrum:
            count += 1
            theo_spectrum.remove(i)
    return count


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


# def parser(filename):
#     filestream = open(filename).read()
#     spectrum = []
#     try:
#         for i in filestream.split('\n')[2].split(' '):
#             spectrum.append(int(i))
#     except ValueError:
#         pass

#     N = int(filestream.split('\n')[1])
#     M = int(filestream.split('\n')[0])
#     return [M, N, spectrum]


# M = 20
# N = 60
# spectrum = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208,
#             228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]


FILENAME = os.path.expanduser('~/Downloads/spectrum25.txt')
# FILENAME = os.path.expanduser('~/Downloads/data_test.txt')
# spectrum = parser(FILENAME)[2]
# spectrum.sort()
# N = parser(FILENAME)[1]
# M = parser(FILENAME)[0]

N = 1000
M = 20
spectrum25 = open(FILENAME).read().split()
spectrum = []
for i in spectrum25:
    spectrum.append(int(i))
print(spectrum)

# print(Convolution(Spectrum))
print(find_ext_alphabet(Convolution(spectrum), M))
# print(ConvolutionCyclopeptideSequencing(M, N, spectrum))


result = ConvolutionCyclopeptideSequencing(M, N, spectrum)

output = ''
for r in result:
    for i in r:
        output += str(i)+'-'
    output = output[:-1]
    output += ' '
output = output[:-1]

f_output = open('data.txt', 'w')
f_output.write(output)
f_output.close()


print(len(result))
# for i in result:
#     print(i)
#     print(CyclopeptideScoring(i, spectrum))


# exp_result = [113, 115, 114, 128, 97, 163, 131, 129, 129, 147, 57, 57, 129]
# print(exp_result)
# print(CyclopeptideScoring(exp_result, spectrum))


# output = ''
# for i in result:
#     output += str(i)+'-'
# print(output[:-1])
# print(result)
# for i in result:
#     print(CyclopeptideScoring(i, Spectrum))
# print(len(result))
