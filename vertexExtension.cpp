/*
    Vertex Extension Code (for (K5, J5) graphs)

    This code basically just tries to add on a vertex to a graph and calculates all possible ways to add edges onto the vertex
    without violating the (K5, J5) restriction.

    Made by : Jason Hao (based off of code made by Elisha Kahan)
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <bit>
#include <bitset>
#include <algorithm>
#include <chrono>

using namespace std;

// Note!! : This does require compiling this in C++20 or later because bit_width() and popcount() requires it

struct Subgraph{
    int type;
    // Regular K5s -> 0
    // Complement K5s -> 1
    // Complement J5s -> 2 

    unsigned int vertices;
};

bool compare_subgraph(Subgraph const &graphA, Subgraph const &graphB)
{
    return graphA.vertices < graphB.vertices;
}

// This turns a string of bits representing vertices to a vector of individual vertices represented as bits
vector<unsigned int> bit_to_vect(unsigned int vertices)
{
    vector<unsigned int> vertexVector;
    int numvertices = popcount(vertices);
    for (int i = 0; i < numvertices; i++)
    {
        unsigned int vertex = 1 << (bit_width(vertices) - 1);
        vertexVector.push_back(vertex);
        vertices ^= vertex;
    }
    return vertexVector;
}

vector<int> bit_to_index(unsigned int vertices)
{
    vector<int> indexVector;
    int numvertices = popcount(vertices);
    for (int i = 0; i < numvertices; i++)
    {
        int index = bit_width(vertices) - 1;
        indexVector.push_back(index + 1);
        vertices ^= 1 << index;
    }
    return indexVector;
}

vector<unsigned int> complement(vector<unsigned int> &graph)
{
    vector<unsigned int> compGraph;
    unsigned int pointer = (1 << (graph.size() - 1));
    for (unsigned int row : graph)
    {
        compGraph.push_back((~row ^ pointer) & ((1 << (graph.size())) - 1));
        pointer >>= 1;
    }
    return compGraph;
}

// seems to break with size 32 and above, but we are guaranteed to not hit that case (otherwise, R(K_5, J_5) is 33, which is not correct)
vector<unsigned int> expand_intervals(vector<vector<unsigned int>> &intervals, int size)
{
    vector<unsigned int> expanded_graphs;
    for (vector<unsigned int> &interval : intervals)
    {
        vector<int> common_vertices = bit_to_index( (interval[0] | ~(interval[1])) & ((1 << size) - 1) );
        for (unsigned int i = 0; i <= (1 << (size - common_vertices.size())) - 1; i++)
        {
            unsigned int fill_in = i;
            for (int j = common_vertices.size() - 1; j >= 0; j--)
            {
                fill_in = ( (fill_in >> (common_vertices[j] - 1)) << (common_vertices[j])) | (fill_in & ((1 << (common_vertices[j] - 1)) - 1) );
            }
            expanded_graphs.push_back(fill_in | interval[0]);
        }

        interval[0] <<= 1;
        interval[1] = (interval[1] << 1) + 1;
    }
    return expanded_graphs;
}

vector<unsigned int> decodeG6(string graph6)
{
    int num_vertices = graph6[0] - 63;
    int size = num_vertices * (num_vertices - 1) / 2;
    vector<unsigned short> bitVector(size + 6, 0); //i want to instantiate the entire array as zeroes

    int index = 0;
    for (int j = 1; j < graph6.length(); j++)
    {
        int bits = graph6[j] - 63;
        for (int k = 5; k >= 0; k--)
        {
            unsigned short bit = (bits & (1 << k)) >> k;
            bitVector[index] = bit;
            index++;
        }
    }

    // This part actually makes the adjacency matrix
    // The binary representation of each unsigned int represents a row of the adjacency matrix
    vector<unsigned int> graph(num_vertices, 0);

    index = 0;
    for (int column = 1; column < num_vertices; column++)
    {
        for (int row = 0; row < column; row++)
        {
            graph[row] |= bitVector[index] << (num_vertices - column - 1);
            graph[column] |= bitVector[index] << (num_vertices - row - 1);
            index += 1;
        }
    }

    return graph;
}

// This is just the Bron-Kerbosch algorithm with pivoting using bit operations, that's it
// Wikipedia link: https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
void max_clique(unsigned int r, unsigned int p, unsigned int x, vector<Subgraph> &maxSubgraphs, vector<unsigned int> &graph, int graph_size, int target_size)
{
    if (p == 0 && x == 0)
    {
        if (popcount(r) == target_size)
        {
            maxSubgraphs.emplace_back(0, r);
        }
        return;
    }

    unsigned int pivot = bit_width(p | x) - 1;

    for (int index : bit_to_index(p & ~(graph[graph_size - pivot - 1])) )
    {
        unsigned int pointer = 1 << (index - 1);
        unsigned int neighbors = graph[graph_size - index];
        max_clique(r | pointer, p & neighbors, x & neighbors, maxSubgraphs, graph, graph_size, target_size);
        p &= ~pointer;
        x |= pointer;
    }
}

void max_set(unsigned int r, unsigned int p, unsigned int x, vector<Subgraph> &maxSubgraphs, vector<unsigned int> &max_sets_n, vector<unsigned int> &graph, int graph_size, int target_size)
{
    if (p == 0 && x == 0)
    {
        if (popcount(r) == target_size)
        {
            maxSubgraphs.emplace_back(1, r);
        }
        else if (popcount(r) == (target_size - 1))
        {
            for (unsigned int max_set : max_sets_n)
            {
                if (popcount(r & max_set) == (target_size - 2))
                {
                    maxSubgraphs.emplace_back(2, r | max_set);
                }
            }
            max_sets_n.push_back(r);
        }
        return;
    }

    unsigned int pivot = bit_width(p | x) - 1;
    for (int index : bit_to_index( p & (graph[graph_size - pivot - 1] | (1 << pivot)) ))
    {
        unsigned int pointer = 1 << (index - 1);
        unsigned int neighbors = ~graph[graph_size - index] ^ pointer;
        max_set(r | pointer, p & neighbors, x & neighbors, maxSubgraphs, max_sets_n, graph, graph_size, target_size);
        p &= ~pointer;
        x |= pointer;
    }
}

vector<Subgraph> max_subgraphs(vector<unsigned int> &graph, int size)
{
    vector<Subgraph> subgraphs;
    vector<unsigned int> max_sets;
    max_clique(0, (1 << size) - 1, 0, subgraphs, graph, size, 4);
    max_set(0, (1 << size) - 1, 0, subgraphs, max_sets, graph, size, 4);
    sort(subgraphs.begin(), subgraphs.end(), compare_subgraph);
    return subgraphs;
}

// we're only checking the subgraphs with the new vertex because all of the intervals already don't contain the old subgraphs
vector<Subgraph> max_subgraphs_last(vector<unsigned int> &old_graph, vector<unsigned int> &graph, int size)
{
    vector<Subgraph> newsubgraphs;
    vector<Subgraph> new_sets;
    vector<unsigned int> max_sets;
    vector<unsigned int> compGraph = complement(old_graph);

    max_clique(0, (1 << (size - 1)) - 1, 0, new_sets, compGraph, size - 1, 3);
    for (Subgraph &subgraph : new_sets)
    {
        max_sets.push_back(subgraph.vertices << 1);
    }

    max_clique(1, graph[size - 1], 0, newsubgraphs, graph, size, 4);
    max_set(1, (~graph[size - 1] & ~1) & ((1 << size) - 1), 0, newsubgraphs, max_sets, graph, size, 4);
    sort(newsubgraphs.begin(), newsubgraphs.end(), compare_subgraph);
    return newsubgraphs;
}

void vertex_extend(vector<vector<unsigned int>> &intervals, vector<Subgraph> &max_subgraphs, vector<unsigned int> &graph, int size)
{
    for (Subgraph &subgraph : max_subgraphs)
    {
        if (intervals.empty())
        {
            break;
        }
        vector<vector<unsigned int>> new_intervals;
        if (subgraph.type == 0)
        {
            for (vector<unsigned int> &interval : intervals)
            {
                if ((subgraph.vertices & interval[1]) == subgraph.vertices)
                {
                    vector<unsigned int> diff_vertices = bit_to_vect(subgraph.vertices & ~interval[0]);
                    int num_diff = diff_vertices.size();
                    if (num_diff >= 1)
                    {
                        new_intervals.push_back({interval[0], interval[1] & ~diff_vertices[0]});
                    }
                    if (num_diff >= 2)
                    {
                        new_intervals.push_back({interval[0] | diff_vertices[0], interval[1] & ~diff_vertices[1]});
                    }
                    if (num_diff >= 3)
                    {
                        new_intervals.push_back({interval[0] | diff_vertices[0] | diff_vertices[1], interval[1] & ~diff_vertices[2]});
                    }
                    if (num_diff == 4)
                    {
                        new_intervals.push_back({interval[0] | diff_vertices[0] | diff_vertices[1] | diff_vertices[2], interval[1] & ~diff_vertices[3]});
                    }
                }
                else
                {
                    new_intervals.push_back(interval);
                }
            }
            intervals = new_intervals;
        }
        else if (subgraph.type == 1)
        {
            for (vector<unsigned int> &interval : intervals)
            {
                int count = popcount(subgraph.vertices & interval[0]);
                if (count == 0)
                {
                    vector<unsigned int> inter_vertices = bit_to_vect(subgraph.vertices & interval[1]);
                    int num_inter = inter_vertices.size();

                    if (num_inter >= 2)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[0] | inter_vertices[1], interval[1]});
                    }
                    if (num_inter >= 3)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[0] | inter_vertices[2], interval[1] & ~inter_vertices[1]});
                        new_intervals.push_back({interval[0] | inter_vertices[1] | inter_vertices[2], interval[1] & ~inter_vertices[0]});
                    }
                    if (num_inter == 4)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[0] | inter_vertices[3], interval[1] & ~(inter_vertices[1] | inter_vertices[2])});
                        new_intervals.push_back({interval[0] | inter_vertices[1] | inter_vertices[3], interval[1] & ~(inter_vertices[0] | inter_vertices[2])});
                        new_intervals.push_back({interval[0] | inter_vertices[2] | inter_vertices[3], interval[1] & ~(inter_vertices[0] | inter_vertices[1])});
                    }
                }
                else if (count == 1)
                {
                    vector<unsigned int> inter_vertices = bit_to_vect(subgraph.vertices & interval[1] & ~interval[0]);
                    int num_inter = inter_vertices.size();
                    
                    if (num_inter >= 1)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[0], interval[1]});
                    }
                    if (num_inter >= 2)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[1], interval[1] & ~inter_vertices[0]});
                    }
                    if (num_inter == 3)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[2], interval[1] & ~(inter_vertices[0] | inter_vertices[1])});
                    }
                }
                else
                {
                    new_intervals.push_back(interval);
                }
            }
            intervals = new_intervals;
        }
        else
        {
            for (vector<unsigned int> &interval : intervals)
            {
                if ((subgraph.vertices & interval[0]) == 0)
                {
                    vector<unsigned int> inter_vertices = bit_to_vect(subgraph.vertices & interval[1]);
                    int num_inter = inter_vertices.size();
                    if (num_inter >= 1)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[0], interval[1]});
                    }
                    if (num_inter >= 2)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[1], interval[1] & ~inter_vertices[0]});
                    }
                    if (num_inter >= 3)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[2], interval[1] & ~(inter_vertices[0] | inter_vertices[1])});
                    }
                    if (num_inter == 4)
                    {
                        new_intervals.push_back({interval[0] | inter_vertices[3], interval[1] & ~(inter_vertices[0] | inter_vertices[1] | inter_vertices[2])});
                    }
                }
                else
                {
                    new_intervals.push_back(interval);
                }
            }
            intervals = new_intervals;
        }
    }
    
    if (!intervals.empty())
    {
        vector<unsigned int> expanded_graphs = expand_intervals(intervals, size);
        for (unsigned int gluing : expanded_graphs)
        {
            vector<unsigned int> new_graph;
            for (int i = 0; i < graph.size(); i++)
            {
                new_graph.push_back( (graph[i] << 1) | ((gluing & (1 << (size - i - 1))) >> (size - i - 1)) );
            }
            new_graph.push_back(gluing << 1);
            
            if (new_graph.size() >= 29)
            {
                cout << "Size: " << new_graph.size() << " Graph: "; //print the graph or something
                for (unsigned int row : new_graph)
                {
                    cout << row << " ";
                }
                cout << "\n";
            }
            vector<vector<unsigned int>> interval_copy = intervals;
            vector<Subgraph> new_subgraphs = max_subgraphs_last(graph, new_graph, size + 1);
            vertex_extend(interval_copy, new_subgraphs, new_graph, size + 1);
        }
    }
}

int main()
{
    ifstream graphFile;
    graphFile.open("14-graphsComp/7-graphs/14-7-0VEgraphs.g6"); // <- put whatever name of the file/path which has your graphs

    string graph6string;
    int index = 0;

    const auto start{chrono::steady_clock::now()};

    while (getline(graphFile, graph6string))
    {
        vector<unsigned int> graph = decodeG6(graph6string.substr(10));
        // Cutting off >>graph6<< header on each line. Otherwise, use graph6string instead of graph6string.substr(10) if no header
        int graph_size = graph.size();
        vector<vector<unsigned int>> intervals = {{0, (1 << graph_size) - 1}};
        vector<Subgraph> subgraphs = max_subgraphs(graph, graph_size);
        vertex_extend(intervals, subgraphs, graph, graph_size);
        if (index % 10000 == 0)
        {
            cout << "Index " << index << " done.\n";
        }
        index++;
    }

    const auto finish{chrono::steady_clock::now()};
    const chrono::duration<double> elapsed_seconds{finish - start};
    cout << "Elapsed Time: " << elapsed_seconds << "\n";

    return 0;
}
/*
FlO[O
Fhogg
FgqPg
FlgGg
FhMIG
FhcYG
FhELO
FtTgO */
