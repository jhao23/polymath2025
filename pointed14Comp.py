import networkx as nx
import pickle

deg6graphsComp = nx.read_graph6("ref_graphs/k3k5e_6.g6")
deg7graphsComp = nx.read_graph6("ref_graphs/k3k5e_7.g6")
deg8graphsComp = nx.read_graph6("ref_graphs/k3k5e_8.g6")
deg9graphsComp = nx.read_graph6("ref_graphs/k3k5e_9.g6")
deg10graphsComp = nx.read_graph6("ref_graphs/k3k5e_10.g6")

deg6graphs = []
deg7graphs = []
deg8graphs = []
deg9graphs = []
deg10graphs = []

for graph in deg6graphsComp:
    deg6graphs.append(nx.complement(graph))
for graph in deg7graphsComp:
    deg7graphs.append(nx.complement(graph))
for graph in deg8graphsComp:
    deg8graphs.append(nx.complement(graph))
for graph in deg9graphsComp:
    deg9graphs.append(nx.complement(graph))
for graph in deg10graphsComp:
    deg10graphs.append(nx.complement(graph))

shared_neighborhood6 = set()
shared_neighborhood7 = set()
shared_neighborhood8 = set()
shared_neighborhood9 = set()
shared_neighborhood10 = set()

num_6_pointed_graphs = num_7_pointed_graphs = num_8_pointed_graphs = num_9_pointed_graphs = num_10_pointed_graphs = 0

graphs = nx.read_graph6("ref_graphs/neighborhood2.g6")

for i in range(len(graphs)):
    G = graphs[i]
    for v in G:
        neighbors = list(nx.neighbors(G, v))
        possible_K = nx.induced_subgraph(G, neighbors)
        if len(neighbors) == 6:
            for shared_neighborhood in deg6graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood6.add(shared_neighborhood)
                    num_6_pointed_graphs += 1
        elif len(neighbors) == 7:
            for shared_neighborhood in deg7graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood7.add(shared_neighborhood)
                    num_7_pointed_graphs += 1
        elif len(neighbors) == 8:
            for shared_neighborhood in deg8graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood8.add(shared_neighborhood)
                    num_8_pointed_graphs += 1
        elif len(neighbors) == 9:
            for shared_neighborhood in deg9graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood9.add(shared_neighborhood)
                    num_9_pointed_graphs += 1
        elif len(neighbors) == 10:
            for shared_neighborhood in deg10graphs:
                if nx.is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood10.add(shared_neighborhood)
                    num_10_pointed_graphs += 1

print(len(shared_neighborhood6), len(shared_neighborhood7), len(shared_neighborhood8), len(shared_neighborhood9), len(shared_neighborhood10))
print(num_6_pointed_graphs, num_7_pointed_graphs, num_8_pointed_graphs, num_9_pointed_graphs, num_10_pointed_graphs)

list_6_neighborhoods = list(shared_neighborhood6)
list_7_neighborhoods = list(shared_neighborhood7)
list_8_neighborhoods = list(shared_neighborhood8)
list_9_neighborhoods = list(shared_neighborhood9)
list_10_neighborhoods = list(shared_neighborhood10)

pointed_graphs_6 = []
pointed_graphs_7 = []
pointed_graphs_8 = []
pointed_graphs_9 = []
pointed_graphs_10 = []

for nbh in list_6_neighborhoods:
    pointed_graphs_6.append([])
for nbh in list_7_neighborhoods:
    pointed_graphs_7.append([])
for nbh in list_8_neighborhoods:
    pointed_graphs_8.append([])
for nbh in list_9_neighborhoods:
    pointed_graphs_9.append([])
for nbh in list_10_neighborhoods:
    pointed_graphs_10.append([])

for i in range(len(graphs)):
    G = graphs[i]
    for v in G:
        neighbors = list(nx.neighbors(G, v))
        possible_K = nx.induced_subgraph(G, neighbors)
        size_neighbors = len(neighbors)
        match size_neighbors :
            case 6:
                for j in range(len(list_6_neighborhoods)):
                    shared_neighborhood = list_6_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_6[j].append( (i, v, G_copy) )
            case 7:
                for j in range(len(list_7_neighborhoods)):
                    shared_neighborhood = list_7_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_7[j].append( (i, v, G_copy) )
            case 8:
                for j in range(len(list_8_neighborhoods)):
                    shared_neighborhood = list_8_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_8[j].append( (i, v, G_copy) )
            case 9:
                for j in range(len(list_9_neighborhoods)):
                    shared_neighborhood = list_9_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_9[j].append( (i, v, G_copy) )
            case 10:
                for j in range(len(list_10_neighborhoods)):
                    shared_neighborhood = list_10_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_10[j].append( (i, v, G_copy) )

test_6_count = test_7_count = test_8_count = test_9_count = test_10_count = 0
for same_nbh in pointed_graphs_6:
    for graph in same_nbh:
        test_6_count += 1
for same_nbh in pointed_graphs_7:
    for graph in same_nbh:
        test_7_count += 1
for same_nbh in pointed_graphs_8:
    for graph in same_nbh:
        test_8_count += 1
for same_nbh in pointed_graphs_9:
    for graph in same_nbh:
        test_9_count += 1
for same_nbh in pointed_graphs_10:
    for graph in same_nbh:
        test_10_count += 1
print(len(pointed_graphs_6), len(pointed_graphs_7), len(pointed_graphs_8), len(pointed_graphs_9), len(pointed_graphs_10))
print(test_6_count, test_7_count, test_8_count, test_9_count, test_10_count)

file_to_write = open("pointed_graphs_6_comp", "wb")
pickle.dump(pointed_graphs_6, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_7_comp", "wb")
pickle.dump(pointed_graphs_7, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_8_comp", "wb")
pickle.dump(pointed_graphs_8, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_9_comp", "wb")
pickle.dump(pointed_graphs_9, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_10_comp", "wb")
pickle.dump(pointed_graphs_10, file_to_write)
file_to_write.close()
