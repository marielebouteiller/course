# from a list of kmer, StringReconstruction finds the dna sequence
# 1st: construction of the DeBruijn graph from the patterns: DeBruijn_from_pattern
# 2nd: find the eulerian path from this DeBruijn graph : EulerianPath
# 3rd: find the sequence from the path: PathToGenome


def StringReconstruction(Patterns):
    return PathToGenome(EulerianPath(DeBruijn_from_patterns(Patterns)))


def DeBruijn_from_patterns(Patterns):
    graph = {}
    for p in Patterns:
        if not p[:-1] in graph.keys():
            graph[p[:-1]] = []
        graph[p[:-1]].append(p[1:])
    return graph


def EulerianPath(graph):  # graph is a dictionary
    cycle = []
    unused_edges = {}
    unused_edges = graph
    startingPoint = IdentifyStartAndStop(graph)[0]
    endingPoint = IdentifyStartAndStop(graph)[1]
    cycle.append(startingPoint)  # i manually do the first loop
    cycle.append(unused_edges[startingPoint][0])
    unused_edges[startingPoint].remove(unused_edges[startingPoint][0])
    if bool(unused_edges[startingPoint]) == False:
        unused_edges.pop(startingPoint)
    while cycle[-1] != endingPoint:
        cycle.append(unused_edges[cycle[-1]][0])
        unused_edges[cycle[-2]].remove(unused_edges[cycle[-2]][0])
        if bool(unused_edges[cycle[-2]]) == False:
            unused_edges.pop(cycle[-2])

    while bool(unused_edges) == True:
        i = 0  # find the index where the node has some unused edges
        while not cycle[i] in unused_edges.keys():
            i += 1
        new_start_node = cycle[i]
        loop = CreateNewCycle(unused_edges, new_start_node)[1]
        del cycle[i]
        cycle[i:i] = loop
    return cycle


def CreateNewCycle(graph1, node):
    cycle = [node]
    cycle.append(graph1[node][0])
    graph1[node].remove(graph1[node][0])
    if bool(graph1[node]) == False:
        graph1.pop(node)
    while cycle[-1] != node:
        cycle.append(graph1[cycle[-1]][0])
        graph1[cycle[-2]].remove(graph1[cycle[-2]][0])
        if bool(graph1[cycle[-2]]) == False:
            graph1.pop(cycle[-2])
    return graph1, cycle


def IdentifyStartAndStop(graph2):
    countInOut = {}
    for k in graph2:
        countInOut[k] = []
        countInOut[k].append(len(graph2[k]))
        countInOut[k].append(0)
    for k in graph2:
        for i in graph2[k]:
            if not i in graph2:
                countInOut[i] = []
                countInOut[i].append(0)
                countInOut[i].append(0)
            countInOut[i][1] += 1
    for k in countInOut:
        if countInOut[k][1]-countInOut[k][0] == -1:
            starting = k
        if countInOut[k][1]-countInOut[k][0] == 1:
            ending = k
    return [starting, ending]


def PathToGenome(path):
    gen = str(path[0])
    k = len(path[0])
    for i in range(1, len(path)):
        gen += str(path[i][k-1])
    return gen


#collection_4mers = ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']
# print(StringReconstruction(collection_4mers))

# collection = open(
#     '/Users/dtj848/Downloads/dataset_203_7.txt').read().split('\n')
# del collection[0]  # the first line is k, length of the nucleotides
# del collection[-1]  # the last line is empty
# f = open('data.txt', 'w')
# f.write(StringReconstruction(collection))
# f.close()

input = ['AAAT','AATG','ACCC','ACGC','ATAC','ATCA','ATGC','CAAA','CACC','CATA','CATC','CCAG','CCCA','CGCT','CTCA','GCAT','GCTC','TACG','TCAC','TCAT','TGCA']
print(StringReconstruction(input))