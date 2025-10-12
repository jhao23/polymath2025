import networkx as nx
import itertools
import pickle
import GluingClassesJ as gcj

graphs = nx.read_graph6("ref_graphs/k4ek5_15.g6")

deg4graphs = nx.read_graph6("ref_graphs/k5k3e_4.g6")
deg5graphs = nx.read_graph6("ref_graphs/k5k3e_5.g6")
deg6graphs = nx.read_graph6("ref_graphs/k5k3e_6.g6")

shared_neighborhood_4 = set()
shared_neighborhood_5 = set()
shared_neighborhood_6 = set()

num_4_pointed_graphs = num_5_pointed_graphs = num_6_pointed_graphs = 0

for i in range(len(graphs)):
    G = graphs[i]
    for v in G:
        neighbors = list(nx.neighbors(G, v))
        possible_K = nx.induced_subgraph(G, neighbors)
        if len(neighbors) == 4:
            for shared_neighborhood in deg4graphs:
                if nx.is_isomorphic(possible_K, shared_neighborhood):
                    shared_neighborhood_4.add(shared_neighborhood)
                    num_4_pointed_graphs += 1
        elif len(neighbors) == 5:
            for shared_neighborhood in deg5graphs:
                if nx.is_isomorphic(possible_K, shared_neighborhood):
                    shared_neighborhood_5.add(shared_neighborhood)
                    num_5_pointed_graphs += 1
        else:
            for shared_neighborhood in deg6graphs:
                if nx.is_isomorphic(possible_K, shared_neighborhood):
                    shared_neighborhood_6.add(shared_neighborhood)
                    num_6_pointed_graphs += 1

print(len(shared_neighborhood_4), len(shared_neighborhood_5), len(shared_neighborhood_6))
print(num_4_pointed_graphs, num_5_pointed_graphs, num_6_pointed_graphs)

list_4_neighborhoods = list(shared_neighborhood_4)
list_5_neighborhoods = list(shared_neighborhood_5)
list_6_neighborhoods = list(shared_neighborhood_6)

pointed_graphs_4 = []
pointed_graphs_5 = []
pointed_graphs_6 = []

for nbh in list_4_neighborhoods:
    pointed_graphs_4.append([])
for nbh in list_5_neighborhoods:
    pointed_graphs_5.append([])
for nbh in list_6_neighborhoods:
    pointed_graphs_6.append([])

for i in range(len(graphs)):
    G = graphs[i]
    for v in G:
        neighbors = list(nx.neighbors(G, v))
        possible_K = nx.induced_subgraph(G, neighbors)
        if len(neighbors) == 4:
            for j in range(len(list_4_neighborhoods)):
                shared_neighborhood = list_4_neighborhoods[j]
                if nx.is_isomorphic(possible_K, shared_neighborhood):
                    isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                    for original_v in isomorphism:
                        isomorphism[original_v] = "k" + str(isomorphism[original_v])
                    G_copy = nx.relabel_nodes(G, isomorphism)
                    pointed_graphs_4[j].append((i, v, G_copy))
        elif len(neighbors) == 5:
            for j in range(len(list_5_neighborhoods)):
                shared_neighborhood = list_5_neighborhoods[j]
                if nx.is_isomorphic(possible_K, shared_neighborhood):
                    isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                    for original_v in isomorphism:
                        isomorphism[original_v] = "k" + str(isomorphism[original_v])
                    G_copy = nx.relabel_nodes(G, isomorphism)
                    pointed_graphs_5[j].append((i, v, G_copy))
        else:
            for j in range(len(list_6_neighborhoods)):
                shared_neighborhood = list_6_neighborhoods[j]
                if nx.is_isomorphic(possible_K, shared_neighborhood):
                    isomorphism = nx.vf2pp_isomorphism(possible_K, shared_neighborhood)
                    for original_v in isomorphism:
                        isomorphism[original_v] = "k" + str(isomorphism[original_v])
                    G_copy = nx.relabel_nodes(G, isomorphism)
                    pointed_graphs_6[j].append((i, v, G_copy))
test_4_count = 0
test_5_count = 0
test_6_count = 0
for same_nbh in pointed_graphs_4:
    for graph in same_nbh:
        test_4_count += 1
for same_nbh in pointed_graphs_5:
    for graph in same_nbh:
        test_5_count += 1
for same_nbh in pointed_graphs_6:
    for graph in same_nbh:
        test_6_count += 1

print(len(pointed_graphs_4), len(pointed_graphs_5), len(pointed_graphs_6))
print(test_4_count, test_5_count, test_6_count)

def glue(G, a, H, b, M, g_map, h_map):
    # Copy G and H into a new graph
    glued_graph = G.copy()
    glued_graph.add_nodes_from(H.nodes())
    glued_graph.add_edges_from(H.edges())

    # Connect b to all vertices in G and a to all vertices in H
    for g_vertex in G.nodes():
        glued_graph.add_edge(b, g_vertex)

    for h_vertex in H.nodes():
        glued_graph.add_edge(a, h_vertex)

    # Add edges between vertices of G and H based on a succesful gluing represented by M
    for x in g_map:
        for y in h_map:
            if M.matrix[g_map[x]][h_map[y]].exists == gcj.EdgeExists.TRUE:
                glued_graph.add_edge(x, y)

    return glued_graph


def new_no_edge_glue(G, a_old, H, b_old, perm):
    G_k = list(nx.neighbors(G,a_old))
    map1 = {a_old:"a"}
    g_map = {}
    # for i in range(len(G_k)):
    #     map1[G_k[i]] = "k"+str(i)
    index = 0
    for i in G.nodes():
        if i not in G_k and i != a_old:
            new_name = "g"+str(i)
            map1[i] = new_name
            g_map[new_name] = index
            index += 1
    nx.relabel_nodes(G, map1, False)


    H_k = list(nx.neighbors(H,b_old))
    map2 = {b_old:"b"}
    h_map = {}
    index = 0
    for i in H_k:
        map2[i] = perm[i]
    for i in H.nodes():
        if i not in H_k and i != b_old:
            new_name = "h"+str(i)
            map2[i] = new_name
            h_map[new_name] = index
            index += 1
    H_new = nx.relabel_nodes(H, map2)


    M = gcj.PotentialEdgeMatrix(len(g_map), len(h_map))
    for i in range(len(g_map)):
        for j in range(len(h_map)):
            M.matrix[i][j].set_exists(gcj.EdgeExists.FALSE)

    glued_graph = glue(G, "a", H_new, "b", M, g_map, h_map)
    return glued_graph, G, H_new

real_gluing_problems = []
all_pointed_graphs = [pointed_graphs_4, pointed_graphs_5, pointed_graphs_6]

for same_deg in all_pointed_graphs:
    for same_nbh in same_deg:
        for pair in itertools.combinations_with_replacement(same_nbh, 2):
            first_pointed_graph = pair[0]
            second_pointed_graph = pair[1]
            first_v = first_pointed_graph[1]
            second_v = second_pointed_graph[1]
            G = first_pointed_graph[2]
            H = second_pointed_graph[2]

            a_neighbors = list(nx.neighbors(G, first_v))
            G_k = nx.induced_subgraph(G, a_neighbors)
            b_neighbors = list(nx.neighbors(H, second_v))
            H_k = nx.induced_subgraph(H, b_neighbors)

            all_automorphisms = nx.vf2pp_all_isomorphisms(H_k, G_k)
            for automorphism in all_automorphisms:
                G_copy = G.copy()
                H_copy = H.copy()
                glued_graph, G_lbl, H_lbl = new_no_edge_glue(G_copy, first_v, H_copy, second_v, automorphism)
                #print(G_lbl.nodes(), H_lbl.nodes())
                is_new = True
                for old_problems in real_gluing_problems:
                    if nx.vf2pp_is_isomorphic(glued_graph, old_problems[3]):
                        is_new = False
                        break
                if is_new:
                    real_gluing_problems.append((G_lbl, H_lbl, automorphism, glued_graph))
    print(len(real_gluing_problems))

file_to_write = open("15gluingproblems", "wb")
pickle.dump(real_gluing_problems, file_to_write)
file_to_write.close()

