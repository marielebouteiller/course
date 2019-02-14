# from a collection of (k,d)mers, reconstruct the deBruijn graph, and find a path compatible with the pairs
# I check in the deBruijnPathFromGraph function, that the new node is compatible for both sequences, but I did not add that on the CreateNewcycle function, so not properly done
# suggest in the course to find an eulerian path, and then check if the path is comptatible with both pair of each (k,d)mer, ie first_seq and second_seq are overlapping


def GenomeFromPath(path, d):  # path is a list of tuples
    first_seq = str(path[0][0])
    for i in range(1, len(path)):
        first_seq += str(path[i][0][-1])
    second_seq = str(path[0][1])
    for i in range(1, len(path)):
        second_seq += str(path[i][1][-1])
    k = len(path[0][0])
    return first_seq[:(k+d+1)]+second_seq


# graph is a dictionary, find the pah and create the sequence together
def deBruijnPathFromGraph(graph, d):
    cycle = []  # list of tuples
    unused_edges = {}
    unused_edges = graph
    startingPoint = IdentifyStartAndStop(graph)[0]
    endingPoint = IdentifyStartAndStop(graph)[1]
    k = len(startingPoint[0]) + 1
    cycle.append(startingPoint)  # i manually do the first loop
    first_seq = startingPoint[0]
    second_seq = startingPoint[1]
    # d!=0, so i do not have any constraint on the first edge
    cycle.append(unused_edges[startingPoint][0])
    first_seq += unused_edges[startingPoint][0][0][-1]
    second_seq += unused_edges[startingPoint][0][1][-1]
    unused_edges[startingPoint].remove(unused_edges[startingPoint][0])
    if bool(unused_edges[startingPoint]) == False:
        unused_edges.pop(startingPoint)

    for i in range(d):  # before i need to check the consistency between the 2 sequences
        # while cycle[-1] != endingPoint:
        first_seq += unused_edges[cycle[-1]][0][0][-1]
        second_seq += unused_edges[cycle[-1]][0][1][-1]
        cycle.append(unused_edges[cycle[-1]][0])
        unused_edges[cycle[-2]].remove(unused_edges[cycle[-2]][0])
        if bool(unused_edges[cycle[-2]]) == False:
            unused_edges.pop(cycle[-2])

    while cycle[-1] != endingPoint:
        j = 0
        while True:
            # check if consistent with the 2nd kmer of the pair
            first_seq += unused_edges[cycle[-1]][j][0][-1]
            second_seq += unused_edges[cycle[-1]][j][1][-1]
            cycle.append(unused_edges[cycle[-1]][j])
            if first_seq[-1] == second_seq[-(k+d+1)]:
                break
            else:
                cycle.pop(-1)
                first_seq = first_seq[:-1]
                second_seq = second_seq[:-1]
                j += 1

        unused_edges[cycle[-2]].remove(unused_edges[cycle[-2]][j])
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
# graph with keys as tuple and values as lists of tuple


# same as for the StringeReconstruction from collection of normal kmers
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

# to get the graph from the list of pairs


def deBruijnGraphFromPairs(kmer_pairs):
    deBruijn_pairs = {}
    for l in kmer_pairs:
        if not (Prefix(l[0]), Prefix(l[1])) in deBruijn_pairs.keys():
            deBruijn_pairs[(Prefix(l[0]), Prefix(l[1]))] = []
        deBruijn_pairs[(Prefix(l[0]), Prefix(l[1]))].append(
            (Suffix(l[0]), Suffix(l[1])))
    return deBruijn_pairs


def Prefix(Pattern):
    return Pattern[0:len(Pattern)-1]


def Suffix(Pattern):
    return Pattern[1:len(Pattern)]


def ReadPairs(text):  # to get a list from the list of k,d-mers, each element of the list is a list with the 2 kmer
    list_pairs = []
    for l in text:
        a, b = l.split('|')
        pair = []
        pair.append(a)
        pair.append(b)
        list_pairs.append(pair)
    return list_pairs


text = open('/Users/dtj848/Downloads/dataset_204_16.txt').read().split('\n')
del text[0]
del text[-1]
d = 200
readingFile = ReadPairs(text)
graph = deBruijnGraphFromPairs(ReadPairs(text))
path = deBruijnPathFromGraph(graph, d)
genome = GenomeFromPath(path, d)

f = open('data.txt', 'w')
f.write(genome)
f.close()

# TryText = ['GAGA|TTGA', 'TCGT|GATG', 'CGTG|ATGT', 'TGGT|TGAG','GTGA|TGTT', 'GTGG|GTGA', 'TGAG|GTTG', 'GGTC|GAGA', 'GTCG|AGAT']
# print(ReadPairs(TryText))
# print(TryText)
# graph = deBruijnGraphFromPairs(ReadPairs(TryText))
# print(graph)
# print(IdentifyStartAndStop(graph))

# path = deBruijnPathFromGraph(graph, 2)
# print(path)
# genome = GenomeFromPath(path, 2)
# print(genome)

input=['ACC|ATA','ACT|ATT','ATA|TGA','ATT|TGA','CAC|GAT','CCG|TAC','CGA|ACT','CTG|AGC','CTG|TTC','GAA|CTT','GAT|CTG','GAT|CTG','TAC|GAT','TCT|AAG','TGA|GCT','TGA|TCT','TTC|GAA']
graph = deBruijnGraphFromPairs(ReadPairs(input))
path = deBruijnPathFromGraph(graph, 1)
print(path)
genome = GenomeFromPath(path, 1)
print(genome)