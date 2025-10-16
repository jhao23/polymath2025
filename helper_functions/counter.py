import networkx as nx
import os

total_graphs = 0
total_checked = 0

#The filenames and ranges can be altered for the glued graph files
for i in range(0, 1807):
    problem_count = 0
    checked_count = 0
    for j in range(0, 25000):
        filename1 = "../14-graphsComp/14-9-" + str(i) + "file" + str(j) + ".g6"

        if os.path.exists(filename1):
            
            with open(filename1, "rb") as file1:
                for line in file1:
                    problem_count += 1

    print("File " + str(i) + ":", problem_count, end="")
    
    filename2 = "../14-graphsComp/9-graphs/14-9-" + str(i) + "VEgraphs.g6"
    
    if os.path.exists(filename2):
        with open(filename2, "rb") as file2:
            for line in file2:
                checked_count += 1

    print(" -> " + str(checked_count))

    total_graphs += problem_count
    total_checked += checked_count

print("Total graphs: " + str(total_graphs) + ", After removing isomorphic copies: " + str(total_checked))
