import networkx as nx
import pickle

deg3graphs = nx.read_graph6("ref_graphs/k5k3e_3.g6")
deg4graphs = nx.read_graph6("ref_graphs/k5k3e_4.g6")
deg5graphs = nx.read_graph6("ref_graphs/k5k3e_5.g6")
deg6graphs = nx.read_graph6("ref_graphs/k5k3e_6.g6")
deg7graphs = nx.read_graph6("ref_graphs/k5k3e_7.g6")

shared_neighborhood3 = set()
shared_neighborhood4 = set()
shared_neighborhood5 = set()
shared_neighborhood6 = set()
shared_neighborhood7 = set()

num_3_pointed_graphs = num_4_pointed_graphs = num_5_pointed_graphs = num_6_pointed_graphs = num_7_pointed_graphs = 0

graphs = nx.read_graph6("ref_graphs/k4ek5_14.g6")

for i in range(len(graphs)):
    G = graphs[i]
    for v in G:
        neighbors = list(nx.neighbors(G, v))
        possible_K = nx.induced_subgraph(G, neighbors)
        if len(neighbors) == 3:
            for shared_neighborhood in deg3graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood3.add(shared_neighborhood)
                    num_3_pointed_graphs += 1
        elif len(neighbors) == 4:
            for shared_neighborhood in deg4graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood4.add(shared_neighborhood)
                    num_4_pointed_graphs += 1
        elif len(neighbors) == 5:
            for shared_neighborhood in deg5graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood5.add(shared_neighborhood)
                    num_5_pointed_graphs += 1
        elif len(neighbors) == 6:
            for shared_neighborhood in deg6graphs:
                if nx.vf2pp_is_isomorphic(shared_neighborhood, possible_K):
                    shared_neighborhood6.add(shared_neighborhood)
                    num_6_pointed_graphs += 1
        elif len(neighbors) == 7:
            if nx.is_isomorphic(deg7graphs, possible_K):
                shared_neighborhood7.add(deg7graphs)
                num_7_pointed_graphs += 1

print(len(shared_neighborhood3), len(shared_neighborhood4), len(shared_neighborhood5), len(shared_neighborhood6), len(shared_neighborhood7))
print(num_3_pointed_graphs, num_4_pointed_graphs, num_5_pointed_graphs, num_6_pointed_graphs, num_7_pointed_graphs)

list_3_neighborhoods = list(shared_neighborhood3)
list_4_neighborhoods = list(shared_neighborhood4)
list_5_neighborhoods = list(shared_neighborhood5)
list_6_neighborhoods = list(shared_neighborhood6)
list_7_neighborhoods = list(shared_neighborhood7)

pointed_graphs_3 = []
pointed_graphs_4 = []
pointed_graphs_5 = []
pointed_graphs_6 = []
pointed_graphs_7 = []

for nbh in list_3_neighborhoods:
    pointed_graphs_3.append([])
for nbh in list_4_neighborhoods:
    pointed_graphs_4.append([])
for nbh in list_5_neighborhoods:
    pointed_graphs_5.append([])
for nbh in list_6_neighborhoods:
    pointed_graphs_6.append([])
for nbh in list_7_neighborhoods:
    pointed_graphs_7.append([])

for i in range(len(graphs)):
    G = graphs[i]
    for v in G:
        neighbors = list(nx.neighbors(G, v))
        possible_K = nx.induced_subgraph(G, neighbors)
        size_neighbors = len(neighbors)
        match size_neighbors :
            case 3:
                for j in range(len(list_3_neighborhoods)):
                    shared_neighborhood = list_3_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_3[j].append( (i, v, G_copy) )
            case 4:
                for j in range(len(list_4_neighborhoods)):
                    shared_neighborhood = list_4_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_4[j].append( (i, v, G_copy) )
            case 5:
                for j in range(len(list_5_neighborhoods)):
                    shared_neighborhood = list_5_neighborhoods[j]
                    if nx.is_isomorphic(possible_K, shared_neighborhood):
                        isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                        for original_v in isomorphism:
                            isomorphism[original_v] = "k" + str(isomorphism[original_v])
                        G_copy = nx.relabel_nodes(G, isomorphism)
                        pointed_graphs_5[j].append( (i, v, G_copy) )
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

test_3_count = test_4_count = test_5_count = test_6_count = test_7_count = 0
for same_nbh in pointed_graphs_3:
    for graph in same_nbh:
        test_3_count += 1
for same_nbh in pointed_graphs_4:
    for graph in same_nbh:
        test_4_count += 1
for same_nbh in pointed_graphs_5:
    for graph in same_nbh:
        test_5_count += 1
for same_nbh in pointed_graphs_6:
    for graph in same_nbh:
        test_6_count += 1
for same_nbh in pointed_graphs_7:
    for graph in same_nbh:
        test_7_count += 1
print(len(pointed_graphs_3), len(pointed_graphs_4), len(pointed_graphs_5), len(pointed_graphs_6), len(pointed_graphs_7))
print(test_3_count, test_4_count, test_5_count, test_6_count, test_7_count)

file_to_write = open("pointed_graphs_3", "wb")
pickle.dump(pointed_graphs_3, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_4", "wb")
pickle.dump(pointed_graphs_4, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_5", "wb")
pickle.dump(pointed_graphs_5, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_6", "wb")
pickle.dump(pointed_graphs_6, file_to_write)
file_to_write.close()

file_to_write = open("pointed_graphs_7", "wb")
pickle.dump(pointed_graphs_7, file_to_write)
file_to_write.close()
