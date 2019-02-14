# from a collection of kmer, in the right order with k-1 overlap, PathToGenome reconstructs the DNA sequence


def PathToGenome(path):
    gen = str(path[0])
    k = len(path[0])
    for i in range(1, len(path)):
        gen += str(path[i][k-1])
    return gen


#kmer = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
path = open('dataset_198_3.txt').read().split()
print(PathToGenome(path))

txt = PathToGenome(path)
f = open('data.txt', 'w')
f.write(txt)
f.close()
