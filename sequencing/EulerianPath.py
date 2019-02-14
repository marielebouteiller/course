import re
# Eulerian Path: from a graph -organised with ''->'' lines, find the eulerian path which travels through every node


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


def PrintEulerianPath(eulCycle):
    output_cycle = ''
    for i in eulCycle:
        output_cycle += str(i)+'->'
    return output_cycle[:-2]


def createDictFromGraph(graph):  # import dataset and transform into a dictionnary
    graph_in_dict = {}
    for l in graph:
        if not l:
            continue
        pattern = re.search(r'(\d+)\s->\s(.+)', l)
        graph_in_dict[pattern.group(1)] = pattern.group(2).split(',')
    return graph_in_dict


#input = {0: [2], 1: [3], 2: [1], 3: [0, 4], 6: [3, 7], 7: [8], 8: [9], 9: [6]}
# print(IdentifyStartAndStop(input))
# print(PrintEulerianPath(EulerianPath(input)))

graph = open('/Users/dtj848/Downloads/dataset_203_6.txt').read().split('\n')
graph_dict = createDictFromGraph(graph)
f = open('data.txt', 'w')
f.write(PrintEulerianPath(EulerianPath(graph_dict)))
f.close()
