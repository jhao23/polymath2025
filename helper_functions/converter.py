import networkx as nx
import pynauty as pn
import os

# This program helps compile the graph files produced by gluing programs and removes isomorphic copies of graphs.
# For some reason, this program starts eating up a lot of memory as it goes on, but I can't seem to find a memory leak right now.

def convertNetworkXGraph(nx_graph):
    pn_graph = pn.Graph(len(nx_graph))
    for n in nx_graph:
        list_neighbors = []
        for nbs in nx_graph[n]:
            list_neighbors.append(nbs)
        pn_graph.connect_vertex(n, list_neighbors)
    return pn_graph

totalhits = 0
totalscans = 0

for j in range(1500, 1807): # <- These bounds can be changed to however many files you want to process
    hits = total = 0
    hash_set = set()
    file_to_write = open("../14-graphsComp/9-graphs/14-9-" + str(j) + "VEgraphs.g6", "wb") # <- File name can also be modified
    
    for i in range(0,40000): # <- These bounds can also be changed depending on how many gluing problems were put into one file
        filename = "14-9-"+ str(j)+ "file" + str(i) + ".g6"
        if os.path.exists(filename):
            with open(filename, "rb") as f:
                for line in f:
                    graph = nx.from_graph6_bytes(line[10:-1])
                    pn_graph = convertNetworkXGraph(graph)
                    cert = pn.certificate(pn_graph)
                    total += 1

                    if cert not in hash_set:
                        hits += 1
                        hash_set.add(cert)
                        graph6bytes = nx.to_graph6_bytes(graph)
                        file_to_write.write(graph6bytes)
                    if total % 10000 == 0:
                        print("File " + str(j) + ": " + str(hits) + " out of " + str(total))
            #os.remove(filename) <- If you want to get rid of the intermediate graph files, uncomment here

    print("File " + str(j) + ": " + str(len(hash_set)))
    file_to_write.close()

    if not len(hash_set): # <- Makes sure empty files aren't clogging up the directory
        os.remove("../14-graphsComp/9-graphs/14-9-" + str(j) + "VEgraphs.g6")
    totalscans += total
    totalhits += hits

print("Before: " + str(totalscans) + "   After: " + str(totalhits))
