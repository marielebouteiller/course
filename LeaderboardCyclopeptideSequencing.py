import os

# version used to return all peptides with the highest score


def LeaderboardCyclopeptideSequencing(spectrum, N):
    Leaderboard = [[]]
    LeaderPeptide = []
    ParentMass = Mass_spectrum(spectrum)
    while Leaderboard:
        Leaderboard_totest = Expand(Leaderboard)
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
    return LeaderPeptide

# used to return only one leader peptide
# def LeaderboardCyclopeptideSequencing(spectrum, N):
#     Leaderboard = [[]]
#     LeaderPeptide = []
#     ParentMass = Mass_spectrum(spectrum)
#     while Leaderboard:
#         Leaderboard_totest = Expand(Leaderboard)
#         Leaderboard = Leaderboard_totest.copy()
#         for peptide in Leaderboard_totest:
#             if Mass_peptide(peptide) == ParentMass and LinearScoring(peptide, spectrum) > LinearScoring(LeaderPeptide, spectrum):
#                 # if not copy then the Leader peptide is also expanded at the next step
#                 LeaderPeptide = peptide.copy()

#             elif Mass_peptide(peptide) > ParentMass:
#                 Leaderboard.remove(peptide)
#         if not Leaderboard:  # if Laederboard is empty at that point, I just skip the next Trim obviously
#             break
#         Leaderboard = Trim(Leaderboard, spectrum, N)
#     return LeaderPeptide


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


def Expand(list_of_peptides):
    new_list_of_peptides = []
    for l in list_of_peptides:
        for a in AminoAcidMass:
            l.append(a)
            new_list_of_peptides.append(l)
            l = l[:-1]
    return new_list_of_peptides


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


AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113,
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


# N = 10
# Spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
# print(LeaderboardCyclopeptideSequencing(Spectrum, N))

def parser(filename):
    filestream = open(filename).read()[:-1]
    spectrum = []
    try:
        for i in filestream.split('\n')[1].split(' '):
            spectrum.append(int(i))
    except ValueError:
        pass

    N = int(filestream.split('\n')[0])
    return [N, spectrum]


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


# print(CyclopeptideScoring(result[0], Spectrum))
# print(CyclopeptideScoring(result[1], Spectrum))
# print(LinearScoring(result[0], Spectrum))
# print(LinearScoring(result[1], Spectrum))
FILENAME = os.path.expanduser('~/Downloads/dataset_102_10.txt')
Spectrum = parser(FILENAME)[1]
N = parser(FILENAME)[0]

result = LeaderboardCyclopeptideSequencing(Spectrum, N)
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

output = ''
for i in result:
    output += str(i)+'-'
print(output[:-1])
print(result)
for i in result:
    print(CyclopeptideScoring(i, Spectrum))
print(len(result))
