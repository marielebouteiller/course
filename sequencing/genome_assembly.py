def Composition(Text, k):
    kmer = []
    for i in range(len(Text)-k+1):
        kmer.append(Text[i:i+k])
    return kmer

def PathToGenome(path):
    gen = str(path[0])
    k = len(path[0])
    for i in range(1, len(path)):
        gen += str(path[i][k-1])
    return gen