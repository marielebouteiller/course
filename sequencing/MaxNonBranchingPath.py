# from a collection of kmers, draw a deBruijn graph, find the contigs = maximal non branching path


# collect all the max non branching paths: a list of lists, from list of kmers
def MaximalNonBranchingPath(collection):
    graph = DeBruijn_from_patterns(collection)
    unused_edges = graph.copy()
    InOut = in_out(graph)  # this is a dictionary
    PathsCollection = []
    while bool(unused_edges) == True:
        l = []  # list of nodes != [1,1]
        # pour sortir de cette boucle while quand il ne reste que des loop [1,1]
        for k in unused_edges:
            if InOut[k] != [1, 1]:
                l.append(k)
        if bool(l) == False:
            break
        i = 0
        while InOut[list(unused_edges.keys())[i]] == [1, 1] or InOut[list(unused_edges.keys())[i]][1] == 0:
            i += 1
        node = list(unused_edges.keys())[i]
        path = [node]
        path.append(unused_edges[node][0])
        unused_edges[node].remove(unused_edges[node][0])
        if bool(unused_edges[node]) == False:
            unused_edges.pop(node)
        while path[-1] != node and InOut[path[-1]] == [1, 1]:
            path.append(unused_edges[path[-1]][0])
            unused_edges[path[-2]].remove(unused_edges[path[-2]][0])
            if bool(unused_edges[path[-2]]) == False:
                unused_edges.pop(path[-2])
        PathsCollection.append(path)

    while bool(unused_edges) == True:  # find the loop paths
        node1 = list(unused_edges.keys())[0]
        path = [node1]
        path.append(unused_edges[node1][0])
        unused_edges[node1].remove(unused_edges[node1][0])
        if bool(unused_edges[node1]) == False:
            unused_edges.pop(node1)
        while path[-1] != node1 and InOut[path[-1]] == [1, 1]:
            path.append(unused_edges[path[-1]][0])
            unused_edges[path[-2]].remove(unused_edges[path[-2]][0])
            if bool(unused_edges[path[-2]]) == False:
                unused_edges.pop(path[-2])
        PathsCollection.append(path)

    return PathsCollection


# build the graph (dictionary) from the collection of kmers
def DeBruijn_from_patterns(Patterns):
    graph = {}
    for p in Patterns:
        if not p[:-1] in graph.keys():
            graph[p[:-1]] = []
        graph[p[:-1]].append(p[1:])
    return graph


def in_out(graph2):  # for each node, return the numbers of in-edges and out-edges
    countInOut = {}
    for k in graph2:
        countInOut[k] = []
        countInOut[k].append(0)
        countInOut[k].append(len(graph2[k]))
    for k in graph2:
        for i in graph2[k]:
            if not i in graph2:
                countInOut[i] = []
                countInOut[i].append(0)
                countInOut[i].append(0)
            countInOut[i][0] += 1
    return countInOut


def buildOneNonBranchingPath(graph1, node):
    # until reaching a node != 1,1
    # or until the ending node=first node
    path = [node]
    InOut1 = in_out(graph1)
    path.append(graph1[node][0])
    graph1[node].remove(graph1[node][0])
    if bool(graph1[node]) == False:
        graph1.pop(node)
    while path[-1] != node and InOut1[path[-1]] == [1, 1] and InOut1[path[-1]][1] != 0:
        path.append(graph1[path[-1]][0])
        graph1[path[-2]].remove(graph1[path[-2]][0])
        if bool(graph1[path[-2]]) == False:
            graph1.pop(path[-2])
    return graph1, path


def PathToGenome(paths):  # returns a lists of sequence
    gen = []
    for p in paths:
        contig = str(p[0])
        k = len(p[0])
        for i in range(1, len(p)):
            contig += str(p[i][k-1])
        gen.append(contig)
    gen.sort()
    return gen


def PrintTheContigs(l):
    output = ''
    for i in l:
        output += i
        output += '\n'
    return output[:-1]


input = ['ATG', 'ATG', 'TGT', 'TGG', 'CAT', 'GGA', 'GAT', 'AGA']
input1 = ['ATG', 'ATG', 'TGT', 'TGG', 'CAT',
          'GGA', 'GAT', 'AGA', 'BEF', 'EFH', 'FHB', 'HBE']
input2 = ['12', '23', '34', '35', '67', '76']
# coll = MaximalNonBranchingPath(input1)
# print(coll)
# # p1 = PathToGenome(coll)
# # print(p1)
# print(PrintTheContigs(coll))
print(PrintTheContigs(PathToGenome(MaximalNonBranchingPath(input1))))

# coll = open('/Users/dtj848/Downloads/dataset_205_5.txt').read().split('\n')
# del coll[-1]
# f = open('data.txt', 'w')
# f.write(PrintTheContigs(PathToGenome(MaximalNonBranchingPath(coll))))
# f.close()

# # print(len('data.txt'))
# coll_mine = open('/Users/dtj848/data.txt').read().split('\n')
# print(len(coll_mine))
# coll_ref = open(
#     '/Users/dtj848/dataset_output_reference.txt').read().split('\n')
# print(len(coll_ref))
