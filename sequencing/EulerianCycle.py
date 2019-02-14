import re
# Eulerian cycle: from a graph -organised with ''->'' lines, find one eulerian cycle which travels through every node


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


def EulerianCycle(graph):  # graph is a dictionary
    cycle = []
    unused_edges = {}
    unused_edges = graph
    new_node = list(unused_edges.keys())[0]
    cycle.append(new_node)  # i manually do the first loop
    cycle.append(unused_edges[new_node][0])
    unused_edges[new_node].remove(unused_edges[new_node][0])
    if bool(unused_edges[new_node]) == False:
        unused_edges.pop(new_node)
    while cycle[-1] != new_node:
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


def PrintEulerianCycle(eulCycle):
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


graph = open('/Users/dtj848/Downloads/dataset_203_2.txt').read().split('\n')
graph_dict = createDictFromGraph(graph)
f = open('data.txt', 'w')
f.write(PrintEulerianCycle(EulerianCycle(graph_dict)))
f.close()
