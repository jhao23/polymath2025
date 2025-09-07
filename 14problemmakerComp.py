#will contain a elements of the form (i,j, permutation, Graph) where i,j are indices used, 
#permutation is the automorphism, Graph is the glued graph
import itertools
import pickle
import networkx as nx
import GluingClassesJComp as gcjc
import bz2
# Description: We will use a matrix M to glue two pointed graphs together. 
# Input: Two pointed graphs (G, a) and (H, b) as well as the intersection K. 
#        Also a d' x d' matrix M and g_map, h_map mappings from the matrix to the vertices of G and H. 
# Output: A graph, based on glueing along M. 
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

    # Add edges between vertices of G and H based on a succesfull gluing represented by M
    for x in g_map:
        for y in h_map:
            if M.matrix[g_map[x]][h_map[y]].exists == gcjc.EdgeExists.TRUE:
                glued_graph.add_edge(x, y)
    
    return glued_graph

#function to do the basic gluing with no extra edges
#G,H are nx graphs, a_old and b_old are integers of the vertices that define pointed graphs.
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
    #I think you should not do this in place - double check if this causes issues
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


    M = gcjc.PotentialEdgeMatrix(len(g_map), len(h_map))
    for i in range(len(g_map)):
        for j in range(len(h_map)):
            M.matrix[i][j].set_exists(gcjc.EdgeExists.FALSE)
    
    glued_graph = glue(G, "a", H_new, "b", M, g_map, h_map)
    return glued_graph, G, H_new

def save_problems(problem_dict, deg, count, num_col):
    file_name = "14-problemsComp/14gluingproblems"+str(deg)+"-"+str(count) + ".bz2"
    problem_list = []
    for graph_hash in problem_dict:
        for problem in problem_dict[graph_hash]:
            problem_list.append(problem)
    print("Deg:", deg, " Count:", count, "Num problems:", len(problem_list), "Num Collisions:", num_col)
    
    with bz2.BZ2File(file_name, "wb") as f:
        pickle.dump(problem_list, f)

    return 

'''
file_to_read = open("pointed_graphs_6_comp", "rb")
pointed_graphs_6 = pickle.load(file_to_read)
file_to_read.close()

file_to_read = open("pointed_graphs_7_comp", "rb")
pointed_graphs_7 = pickle.load(file_to_read)
file_to_read.close()
'''
file_to_read = open("pointed_graphs_8_comp", "rb")
pointed_graphs_8 = pickle.load(file_to_read)
file_to_read.close()
'''
file_to_read = open("pointed_graphs_9_comp", "rb")
pointed_graphs_9 = pickle.load(file_to_read)
file_to_read.close()

file_to_read = open("pointed_graphs_10_comp", "rb")
pointed_graphs_10 = pickle.load(file_to_read)
file_to_read.close()
'''
print("Loaded Pointed Graphs!")

real_gluing_problems = []
#all_pointed_graphs = [pointed_graphs_6, pointed_graphs_10, pointed_graphs_9, pointed_graphs_7, pointed_graphs_8]
#deg_list = [6,10,9,7,8]
all_pointed_graphs = [pointed_graphs_8]
deg_list = [8]
for i in range(0,1):
    neighborhood_index = 0
    same_deg = all_pointed_graphs[i]
    deg_num = deg_list[i]
    part_num = 4039
    new_gluing_problems = {}
    num_added = 0
    num_collisions = 0
    for same_nbh in same_deg:
        if (neighborhood_index < 31):
            neighborhood_index += 1
            continue
        pair_index = 0
        for pair in itertools.combinations_with_replacement(same_nbh, 2):
            if (neighborhood_index == 31) and (pair_index < 250000):
                pair_index += 1
                continue #so i stopped at 5200000 pairs, and i need to not process them again
            first_pointed_graph = pair[0]
            second_pointed_graph = pair[1]
            first_vtx = first_pointed_graph[1]
            second_vtx = second_pointed_graph[1]
            G = first_pointed_graph[2]
            H = second_pointed_graph[2]

            a_neighbors = list(nx.neighbors(G, first_vtx))
            G_k = nx.induced_subgraph(G, a_neighbors)
            b_neighbors = list(nx.neighbors(H, second_vtx))
            H_k = nx.induced_subgraph(H, b_neighbors)
            
            all_automorphisms = nx.vf2pp_all_isomorphisms(H_k, G_k)
            for automorphism in all_automorphisms:
                G_copy = G.copy()
                H_copy = H.copy()
                glued_graph, G_lbl, H_lbl = new_no_edge_glue(G_copy, first_vtx, H_copy, second_vtx, automorphism)
                #print(G_lbl.nodes(), H_lbl.nodes())
                # is_new = True
                # for old_problems in new_gluing_problems:
                #     if nx.vf2pp_is_isomorphic(glued_graph, old_problems[3]):
                #         is_new = False
                #         break
                # if is_new:
                #     new_gluing_problems.append((G_lbl, H_lbl, automorphism, glued_graph))
                new_hash = nx.weisfeiler_lehman_graph_hash(glued_graph)
                if new_hash not in new_gluing_problems:
                    new_gluing_problems[new_hash] = [(G_lbl, H_lbl, automorphism, glued_graph)]
                    num_added += 1
                else:
                    is_new = True
                    for old_problems in new_gluing_problems[new_hash]:
                        if nx.vf2pp_is_isomorphic(glued_graph, old_problems[3]):
                            is_new = False
                            break 
                    if is_new:
                        new_gluing_problems[new_hash].append((G_lbl, H_lbl, automorphism, glued_graph))
                        num_added += 1
                        num_collisions += 1
                
                if num_added >= 40000:
                    save_problems(new_gluing_problems, deg_num, part_num, num_collisions)
                    new_gluing_problems = {}
                    num_added = 0
                    num_collisions = 0
                    part_num += 1
            if pair_index % 10000 == 0:
                print("Pair " + str(pair_index) + " done.")
            pair_index += 1
        
        print("Neighborhood " + str(neighborhood_index) + " done processing.")
        neighborhood_index += 1

    if num_added > 0:
        save_problems(new_gluing_problems, deg_num, part_num, num_collisions)
