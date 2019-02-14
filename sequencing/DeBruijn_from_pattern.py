# from a collection of kmers, DeBruijn_from_patterns returns the deBruijn graph, with all the directed edges between kmers


def DeBruijn_from_patterns(Patterns):
    graph = {}
    for p in Patterns:
        if not p[:-1] in graph.keys():
            graph[p[:-1]] = []
        graph[p[:-1]].append(p[1:])
    return graph


#Patterns = ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
# print(DeBruijn_from_patterns(Patterns))


def print_the_collection(dico):
    output = ''
    for k in dico:
        output += k + ' -> ' + dico[k][0]
        if len(dico[k]) > 1:
            for i in range(1, len(dico[k])):
                output += ','+dico[k][i]
        output += '\n'
    return output[:-1]


# print(print_the_collection(DeBruijn_from_patterns(Patterns)))

collectionOfPatterns = open('dataset_200_8.txt').read().split()
txt = print_the_collection(DeBruijn_from_patterns(collectionOfPatterns))
f = open('data.txt', 'w')
f.write(txt)
f.close()
