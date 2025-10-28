#ifndef STARTER_H
#define STARTER_H

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <utility> // For std::move
#include <random>  // For random initialization

using namespace std;

// --- Helper trim (Defined in Header as inline) ---
inline string trim(const string &str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (string::npos == first) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

// --- Graph Node Class ---
class Graph_Node
{
private:
    string Node_Name;
    vector<int> Children;
    vector<string> Parents;
    int nvalues;
    vector<string> values;
    vector<float> CPT;

public:
    // Constructors (Defined in Header)
    Graph_Node() : Node_Name(""), nvalues(0) {}
    Graph_Node(string name, int n, vector<string> vals) : Node_Name(std::move(name)), nvalues(n), values(std::move(vals)) {}

    // Getters (Defined in Header)
    string get_name() const { return Node_Name; }
    vector<int> get_children() const { return Children; }
    vector<string> get_Parents() const { return Parents; }
    vector<float> get_CPT() const { return CPT; }
    int get_nvalues() const { return nvalues; }
    vector<string> get_values() const { return values; }

    // Setters/Modifiers (Defined in Header)
    void set_CPT(const vector<float>& new_CPT) { CPT = new_CPT; }
    void set_Parents(const vector<string>& Parent_Nodes) { Parents = Parent_Nodes; }
    int add_child(int new_child_index) 
    {
        for (int child : Children) if (child == new_child_index) return 0;
        Children.push_back(new_child_index);
        return 1;
    }
    void set_name(const string &nm) { Node_Name = nm; }
    void set_values(const vector<string> &vals) { values = vals; nvalues = (int)vals.size(); }
};

// --- Network Class ---
class network
{
    list<Graph_Node> Pres_Graph;

public:
    // Methods (Defined in Header)
    int addNode(const Graph_Node& node) { Pres_Graph.push_back(node); return 0; }
    int netSize() const { return (int)Pres_Graph.size(); }
    
    // Declarations for complex methods
    int get_index(const string& val_name) const;
    list<Graph_Node>::iterator get_nth_node(int n);
    list<Graph_Node>::iterator search_node(const string& val_name);
    list<Graph_Node> &get_nodes_list() { return Pres_Graph; }
};

// --- Utility: map value to index (Defined in Header as inline) ---
inline int get_value_index(const vector<string> &values, const string &val)
{
    for (int i = 0; i < (int)values.size(); i++)
    {
        if (values[i] == val) return i;
    }
    return -1; // not found, or "?"
}

int network::get_index(const string& val_name) const
{
    int count = 0;
    for (const auto& node : Pres_Graph)
    {
        if (node.get_name() == val_name) return count;
        count++;
    }
    return -1;
}

list<Graph_Node>::iterator  network::get_nth_node(int n)
{
    int count = 0;
    for (list<Graph_Node>::iterator listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
    {
        if (count == n) return listIt;
        count++;
    }
    return Pres_Graph.end();
}

list<Graph_Node>::iterator network::search_node(const string& val_name)
{
    for (list<Graph_Node>::iterator listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
    {
        if (listIt->get_name() == val_name) return listIt;
    }
    return Pres_Graph.end();
}

// --- Variable Elimination Factor Class ---
class Factor {
public:
    vector<int> variables; // Node indices that this factor depends on
    vector<int> card;      // Cardinality (number of values) for each variable in the factor
    vector<int> stride;    // Stride for calculating index from assignment
    vector<float> values;  // The actual probability values
    Factor(){};
    // Initialize from a CPT (node index, network pointer, variable names)
    Factor(int node_idx, network* net);
};

// --- Utilities & Learning Functions (Must be implemented in starter.cpp) ---
network read_network(const char *filename);
void write_network(const char *filename, network &BayesNet);
vector<vector<string>> read_dataset(const string &filename);


// --- CORE INFERENCE OPERATIONS ---
Factor product(const Factor& f1, const Factor& f2);
Factor sum_out(const Factor& f, int var_index);
Factor restrict_factor(const Factor& f, int var_index, int value_index);


// --- CORE EM FUNCTIONS ---
void initialize_random_cpts(network &BayesNet);
float calculate_log_likelihood(const network &BayesNet, const vector<vector<string>> &data);
void perform_em_iteration(network &BayesNet, const vector<vector<string>> &data);

// MAIN LEARNING ENTRY POINT (Multi-Start EM Wrapper)
void learn_parameters(network &BayesNet, const vector<vector<string>> &data);

#endif // STARTER_H


