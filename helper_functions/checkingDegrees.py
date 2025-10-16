import networkx as nx

totalDegreeSet = set()

# you can put whatever graph file to read off all of the possible degrees
with open("k4ek5_13.g6", "rb") as file:
    for line in file:
        graph = nx.from_graph6_bytes(line[:-1])
        degreeSet = set()
        for v in graph.nodes():
            degreeSet.add(len(list(nx.neighbors(graph,v))))
        totalDegreeSet = totalDegreeSet.union(degreeSet)

print(totalDegreeSet)
