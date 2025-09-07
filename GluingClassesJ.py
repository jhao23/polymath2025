import enum
class EdgeExists(enum.Enum):
    UNKNOWN = 0 
    TRUE = 1 
    FALSE = 2 

class ClauseType(enum.Enum):
    #CLIQUE_ALL_EDGES refers to all vertices in G, H, and K all being connected together
    CLIQUE_ALL_EDGES = 0 
    INDEP = 1 
    #CLIQUE_SAME_SET_MISSING_EDGE means there is 1 edge missing - so now can only form J6, and only forms
    # a J6 if ALL edges between G and H are chosen
    CLIQUE_SAME_SET_MISSING_EDGE = 2 


# Define a class to represent variables as described on page 7
class PotentialEdge:
        
    # Construct a new PotentialEdge between vertex G_vert in G and H_vert in H
    def __init__(self, G_vert, H_vert):
        # Set variables
        self.G_vertex = G_vert
        self.H_vertex = H_vert
        # Set current value of the variable
        self.exists = EdgeExists.UNKNOWN
        # Create sets of clique clauses and independent-set clauses
        self.clique_all_edges_clauses = []
        self.ind_set_clauses = []
        self.clique_same_set_clauses = []

    def __str__(self):
        if self.exists == EdgeExists.UNKNOWN:
            return "("+str(self.G_vertex) + "," + str(self.H_vertex)+ "," + "U)"
        elif self.exists == EdgeExists.TRUE:
            return "("+str(self.G_vertex) + "," + str(self.H_vertex)+ "," + "T)"
        else:
            return "("+str(self.G_vertex) + "," + str(self.H_vertex)+ "," + "F)"
        
    # Set the value of the variable
    # NOTE: This should only be used when we change the value from UNKNOWN to TRUE or FALSE
    # INPUT: new_value for exists to be set to
    def set_exists(self, new_value):
        if self.exists != EdgeExists.UNKNOWN:
            raise ValueError("Value was not EdgeExists.Unknown")
        self.exists = new_value
        # Decrease the number of unknown for each clause this variable is in
        for clause in self.clique_all_edges_clauses:
            clause.decr_num_unknown(new_value)
        for clause in self.ind_set_clauses:
            clause.decr_num_unknown(new_value)
        for clause in self.clique_same_set_clauses:
            clause.decr_num_unknown(new_value)
        
        
    # Add a clause that the variable is in
    # INPUT: a clause that this variable is in and that should be added to it's list,
    #        an enum that is the type of clause
    def add_clause(self, clause, clause_type):
        if clause_type == ClauseType.CLIQUE_ALL_EDGES:
           self.clique_all_edges_clauses.append(clause)
        elif clause_type == ClauseType.INDEP:
           self.ind_set_clauses.append(clause)
        else:
           self.clique_same_set_clauses.append(clause)



# Define a class to represent clauses as described on page 7
class Clause:
        
    # Construct a new Clause
    # INPUT: a list of PotentialEdges/variables that are contained in this Clause 
    #        an enum that is the type of clause
    def __init__(self, variables, clause_type):
        # Set variables
        self.potential_edges = variables
        self.num_unknown = len(variables)
        self.clause_type = clause_type
        # Number of variables whose value is undesired 
        #        (i.e. number of TRUEs if it's a clique clause, FALSEs if an independent set clause)
        self.num_undesired = 0 
            
        # Add clauses to each potential edge
        for pot_edge in variables:
            pot_edge.add_clause(self, clause_type)

    def __str__(self):
        output = ""
        for var in self.potential_edges:
            output += str(var) + ", "
        output += str(self.clause_type)
        return output
        
        
    # Decrease the number of unknowns, called when a PotentialEdge's value goes from UNKNOWN to TRUE or FALSE
    # INPUT: the new value that the PotentialEdge was changed to
    def decr_num_unknown(self, new_value):
        # Decrease the number of unknowns by 1
        self.num_unknown -= 1 
        # Update number of variables with undesired value
        if self.clause_type == ClauseType.INDEP and new_value == EdgeExists.FALSE:
            self.num_undesired += 1 
        elif self.clause_type != ClauseType.INDEP and new_value == EdgeExists.TRUE:
            self.num_undesired += 1 
        
        
    # Determine whether all PotentialEdges in the clause are set to a value (i.e. are not UNKNOWN)
    # OUTPUT: a boolean, True if there are no UNKNOWNs, False otherwise
    def is_full(self):
        return self.num_unknown == 0 
        
    # Determine whether the clique causes a FAIL state
    #       i.e. if all variables are TRUE when this is a clique clause or all variables FALSE when this is an independent set clause
    # OUTPUT: a boolean, True if in a FAIL state, False otherwise
    def in_fail_state(self):
        if self.clause_type == ClauseType.INDEP:
            return self.num_unknown == 0  and self.num_undesired == len(self.potential_edges)
        elif self.clause_type == ClauseType.CLIQUE_ALL_EDGES:
            #if there are k potential edges, and k-1 set to True, then 1 is set to true - this is a J
            return self.num_unknown <= 1 and self.num_undesired >= len(self.potential_edges) - 1
        else:
            return self.num_unknown == 0  and self.num_undesired == len(self.potential_edges)

    def is_satisfied(self):
        if self.clause_type == ClauseType.CLIQUE_ALL_EDGES:
            return self.num_undesired + self.num_unknown <= len(self.potential_edges) - 2
        else:
            return self.num_undesired + self.num_unknown <= len(self.potential_edges) - 1

# Construct a class to represent the matrix of variables
class PotentialEdgeMatrix:
    
    # Construct a new Matrix of Potential Edges as described in the paper
    # INPUT: num_rows = |VG|-|VK|-1, num_cols = |VH|-|VK|-1
    def __init__(self, num_rows, num_cols):
        # Create the matrix
        # NOTE: Should probably change to a more efficient data structure, since python lists are LinkedLists
        self.matrix = []
        for row in range(0 , num_rows):
            current_row = []
            for col in range(0 , num_cols):
                current_row.append(PotentialEdge(row, col))
            self.matrix.append(current_row)
    def __str__(self):
        output = "["
        for list in self.matrix:
            output += "["
            for elem in list:
                output += str(elem)
                output += ", "
            output += "],"
        output += "]"

        return output


