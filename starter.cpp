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

using namespace std;

// --------------------------- Helper trim ---------------------------
string trim(const string &str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (string::npos == first)
    {
        return "";
    }
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

// --------------------------- Graph Node ---------------------------
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
    Graph_Node()
    {
        Node_Name = "";
        nvalues = 0;
    }

    Graph_Node(string name, int n, vector<string> vals)
    {
        Node_Name = name;
        nvalues = n;
        values = vals;
    }

    string get_name()
    {
        return Node_Name;
    }

    vector<int> get_children()
    {
        return Children;
    }

    vector<string> get_Parents()
    {
        return Parents;
    }

    vector<float> get_CPT()
    {
        return CPT;
    }

    int get_nvalues()
    {
        return nvalues;
    }

    vector<string> get_values()
    {
        return values;
    }

    void set_CPT(vector<float> new_CPT)
    {
        CPT.clear();
        CPT = new_CPT;
    }

    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents = Parent_Nodes;
    }

    int add_child(int new_child_index)
    {
        for (int i = 0; i < (int)Children.size(); i++)
        {
            if (Children[i] == new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }

    void set_name(const string &nm)
    {
        Node_Name = nm;
    }

    void set_values(const vector<string> &vals)
    {
        values = vals;
        nvalues = vals.size();
    }
};

// --------------------------- Network ---------------------------
class network
{
    list<Graph_Node> Pres_Graph;

public:
    int addNode(Graph_Node node)
    {
        Pres_Graph.push_back(node);
        return 0;
    }

    list<Graph_Node>::iterator getNode(int i)
    {
        int count = 0;
        list<Graph_Node>::iterator listIt;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (count++ == i)
                break;
        }
        return listIt;
    }

    int netSize()
    {
        return Pres_Graph.size();
    }

    int get_index(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (listIt->get_name().compare(val_name) == 0)
                return count;
            count++;
        }
        return -1;
    }

    list<Graph_Node>::iterator get_nth_node(int n)
    {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (count == n)
                return listIt;
            count++;
        }
        return listIt;
    }

    list<Graph_Node>::iterator search_node(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (listIt->get_name().compare(val_name) == 0)
                return listIt;
        }
        // if not found, return end iterator (but here we return last)
        // To be safe, print message and return end-like iterator
        // We will return Pres_Graph.end() if not found
        return Pres_Graph.end();
    }

    // helper to allow iteration in write_network
    list<Graph_Node> &get_nodes_list()
    {
        return Pres_Graph;
    }
};

// --------------------------- Read .bif network ---------------------------
network read_network(const char *filename)
{
    network BayesNet;
    string line;
    ifstream myfile(filename);

    if (!myfile.is_open())
    {
        cout << "Error: Could not open file " << filename << endl;
        return BayesNet;
    }

    while (getline(myfile, line))
    {
        line = trim(line);

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        stringstream ss(line);
        string token;
        ss >> token;

        // Parse variable declarations
        if (token == "variable")
        {
            string var_name;
            ss >> var_name;
            var_name = trim(var_name);

            // read lines until we find "type" line containing values between { }
            string type_line;
            while (getline(myfile, type_line))
            {
                type_line = trim(type_line);
                if (type_line.find("type") != string::npos)
                    break;
            }

            // extract values inside { ... }
            vector<string> values;
            size_t open_br = type_line.find('{');
            size_t close_br = type_line.find('}');
            if (open_br != string::npos && close_br != string::npos && close_br > open_br)
            {
                string inside = type_line.substr(open_br + 1, close_br - open_br - 1);
                stringstream vals(inside);
                string tokenv;
                while (getline(vals, tokenv, ','))
                {
                    tokenv = trim(tokenv);
                    if (tokenv.size() > 0)
                    {
                        // remove possible surrounding quotes
                        if (tokenv.front() == '"' && tokenv.back() == '"')
                        {
                            tokenv = tokenv.substr(1, tokenv.size() - 2);
                        }
                        values.push_back(tokenv);
                    }
                }
            }

            Graph_Node new_node(var_name, (int)values.size(), values);
            BayesNet.addNode(new_node);

            // consume the closing brace of variable block if present
            // (some BIFs have "}" after)
            // no-op - continue reading
        }
        // Parse probability tables
        else if (token == "probability")
        {
            // read full header until '{'
            string header = line;
            while (header.find('{') == string::npos)
            {
                string nextline;
                if (!getline(myfile, nextline))
                    break;
                header += " " + trim(nextline);
            }

            // find node name and parents inside parentheses
            size_t lp = header.find('(');
            size_t rp = header.find(')');
            if (lp == string::npos || rp == string::npos || rp <= lp)
                continue;
            string inside = header.substr(lp + 1, rp - lp - 1);
            inside = trim(inside);

            // split by '|' if parents exist
            string node_name;
            vector<string> parents;
            size_t pipe = inside.find('|');
            if (pipe == string::npos)
            {
                // no parents
                node_name = trim(inside);
            }
            else
            {
                node_name = trim(inside.substr(0, pipe));
                string pstr = inside.substr(pipe + 1);
                stringstream ps(pstr);
                string tmp;
                while (ps >> tmp)
                {
                    // remove commas
                    if (!tmp.empty() && tmp.back() == ',')
                        tmp.pop_back();
                    tmp = trim(tmp);
                    if (!tmp.empty())
                        parents.push_back(tmp);
                }
            }

            // find node in BayesNet
            int index = BayesNet.get_index(node_name);
            if (index == -1)
            {
                // node not yet present (possible ordering differences) -> create placeholder
                Graph_Node placeholder(node_name, 0, vector<string>());
                BayesNet.addNode(placeholder);
                index = BayesNet.get_index(node_name);
            }

            list<Graph_Node>::iterator listIt = BayesNet.get_nth_node(index);
            // attach parents to node and register child's index in parents' children lists
            listIt->set_Parents(parents);
            for (auto &pname : parents)
            {
                int pidx = BayesNet.get_index(pname);
                if (pidx == -1)
                {
                    // parent node might appear later - create placeholder
                    Graph_Node placeholder(pname, 0, vector<string>());
                    BayesNet.addNode(placeholder);
                    pidx = BayesNet.get_index(pname);
                }
                list<Graph_Node>::iterator parentIt = BayesNet.get_nth_node(pidx);
                parentIt->add_child(index);
            }

            // Now read CPT lines until "};"
            vector<float> cpt;
            string cpline;
            while (getline(myfile, cpline))
            {
                string t = trim(cpline);
                if (t.empty())
                    continue;
                if (t == "};" || t == "};")
                    break;

                // if line contains "table", take tokens after it
                size_t tablepos = t.find("table");
                string probs_part;
                if (tablepos != string::npos)
                {
                    probs_part = t.substr(tablepos + 5);
                }
                else
                {
                    // if line has ')', take part after ')'
                    size_t closep = t.find(')');
                    if (closep != string::npos)
                    {
                        probs_part = t.substr(closep + 1);
                    }
                    else
                    {
                        probs_part = t;
                    }
                }

                // Extract numeric tokens (could be -1.000000, 0.1234, etc.)
                string tokenp;
                stringstream ssp(probs_part);
                while (ssp >> tokenp)
                {
                    // strip trailing commas and semicolons
                    while (!tokenp.empty() && (tokenp.back() == ',' || tokenp.back() == ';'))
                        tokenp.pop_back();

                    // remove any surrounding quotes
                    if (!tokenp.empty() && tokenp.front() == '"' && tokenp.back() == '"')
                    {
                        tokenp = tokenp.substr(1, tokenp.size() - 2);
                    }

                    // token might include parentheses - skip them
                    if (tokenp == "(" || tokenp == ")")
                        continue;

                    // check if token is numeric (- or digits or .)
                    if (!tokenp.empty() && (isdigit(tokenp[0]) || tokenp[0] == '-' || tokenp[0] == '.'))
                    {
                        float val = atof(tokenp.c_str());
                        cpt.push_back(val);
                    }
                }
            }

            listIt->set_CPT(cpt);
        }
    }

    myfile.close();
    return BayesNet;
}

// --------------------------- Write network back to file ---------------------------
void write_network(const char *filename, network &BayesNet)
{
    ofstream outfile(filename);

    if (!outfile.is_open())
    {
        cout << "Error: Could not open file " << filename << " for writing" << endl;
        return;
    }

    outfile << "// Bayesian Network" << endl
            << endl;

    // Write all nodes first
    int idx = 0;
    list<Graph_Node> &nodes = BayesNet.get_nodes_list();

    // Ensure ordering: write in stored order
    for (auto it = nodes.begin(); it != nodes.end(); ++it)
    {
        outfile << "variable " << it->get_name() << " {" << endl;
        outfile << "  type discrete [ " << it->get_nvalues() << " ] = { ";
        vector<string> vals = it->get_values();
        for (int j = 0; j < (int)vals.size(); j++)
        {
            outfile << "\"" << vals[j] << "\"";
            if (j < (int)vals.size() - 1)
                outfile << ", ";
        }
        outfile << " };" << endl;
        outfile << "}" << endl
                << endl;
        idx++;
    }

    outfile << fixed << setprecision(4);

    // Write probability tables
    for (auto it = nodes.begin(); it != nodes.end(); ++it)
    {
        vector<string> parents = it->get_Parents();
        vector<string> values = it->get_values();
        vector<float> cpt = it->get_CPT();

        outfile << "probability ( " << it->get_name();
        if (!parents.empty())
        {
            outfile << " | ";
            for (int j = 0; j < (int)parents.size(); j++)
            {
                outfile << parents[j];
                if (j < (int)parents.size() - 1)
                    outfile << ", ";
            }
        }
        outfile << " ) {" << endl;

        // compute parent radices
        vector<int> radices;
        for (auto &pname : parents)
        {
            auto pIt = BayesNet.search_node(pname);
            if (pIt != BayesNet.get_nodes_list().end())
                radices.push_back(pIt->get_nvalues());
            else
                radices.push_back(1);
        }

        int parent_combinations = 1;
        for (int r : radices)
            parent_combinations *= r;
        if (parent_combinations == 0)
            parent_combinations = 1;

        int cpt_index = 0;

        if (parents.empty())
        {
            outfile << "    table ";
            for (int k = 0; k < (int)values.size(); k++)
            {
                if (cpt_index < (int)cpt.size())
                    outfile << cpt[cpt_index++];
                else
                    outfile << -1.0;
                if (k < (int)values.size() - 1)
                    outfile << ", ";
            }
            outfile << ";" << endl;
        }
        else
        {
            for (int comb = 0; comb < parent_combinations; comb++)
            {
                // compute mixed radix digits
                vector<int> idxs(parents.size(), 0);
                int tmp = comb;
                for (int p = (int)parents.size() - 1; p >= 0; p--)
                {
                    if (radices[p] > 0)
                    {
                        idxs[p] = tmp % radices[p];
                        tmp /= radices[p];
                    }
                    else
                        idxs[p] = 0;
                }

                // print parent values
                outfile << "    ( ";
                for (int p = 0; p < (int)parents.size(); p++)
                {
                    auto pIt = BayesNet.search_node(parents[p]);
                    if (pIt != BayesNet.get_nodes_list().end())
                    {
                        vector<string> pvals = pIt->get_values();
                        int vidx = idxs[p] % max(1, (int)pvals.size());
                        outfile << pvals[vidx];
                    }
                    else
                    {
                        outfile << "UNK";
                    }
                    if (p < (int)parents.size() - 1)
                        outfile << ", ";
                }
                outfile << " ) ";

                for (int k = 0; k < (int)values.size(); k++)
                {
                    if (cpt_index < (int)cpt.size())
                        outfile << cpt[cpt_index++];
                    else
                        outfile << -1.0;
                    if (k < (int)values.size() - 1)
                        outfile << ", ";
                }
                outfile << ";" << endl;
            }
        }

        outfile << "};" << endl
                << endl;
    }

    outfile.close();
    cout << "Network written to file: " << filename << endl;
}

// --------------------------- Read dataset ---------------------------
// Expect each record as space-separated tokens (assignment said records.dat has space separated)
// Remove quotes if present. Each row must have same number of columns as network nodes.
vector<vector<string>> read_dataset(const string &filename)
{
    ifstream fin(filename);
    vector<vector<string>> data;
    if (!fin.is_open())
    {
        cout << "Error: could not open data file " << filename << endl;
        return data;
    }

    string line;
    while (getline(fin, line))
    {
        line = trim(line);
        if (line.empty())
            continue;

        vector<string> row;
        string token;
        stringstream ss(line);

        // split by whitespace
        while (ss >> token)
        {
            // remove surrounding quotes and trailing commas
            if (!token.empty() && token.front() == '"')
                token.erase(token.begin());
            if (!token.empty() && token.back() == '"')
                token.pop_back();
            if (!token.empty() && token.back() == ',')
                token.pop_back();
            token = trim(token);
            row.push_back(token);
        }
        if (!row.empty())
            data.push_back(row);
    }
    fin.close();
    return data;
}

// --------------------------- Utility: map value to index ---------------------------
int get_value_index(const vector<string> &values, const string &val)
{
    for (int i = 0; i < (int)values.size(); i++)
    {
        if (values[i] == val)
            return i;
    }
    return -1; // not found, or "?"
}

// --------------------------- Learning: MLE ---------------------------
void learn_parameters(network &BayesNet, const vector<vector<string>> &data)
{
    int N = BayesNet.netSize();
    int num_records = data.size();
    cout << "Learning parameters from " << num_records << " records..." << endl;

    // For each node
    for (int node_idx = 0; node_idx < N; node_idx++)
    {
        auto it = BayesNet.get_nth_node(node_idx);
        string node_name = it->get_name();
        vector<string> node_vals = it->get_values();
        vector<string> parents = it->get_Parents();

        // prepare parent value lists
        vector<vector<string>> parent_vals;
        for (auto &p : parents)
        {
            auto pIt = BayesNet.search_node(p);
            if (pIt != BayesNet.get_nodes_list().end())
                parent_vals.push_back(pIt->get_values());
            else
                parent_vals.push_back(vector<string>(1, "UNK"));
        }

        // Counting structures:
        // key = concatenated parent values (with '|'), counts[key][node_value] = freq
        map<string, map<string, int>> counts;
        map<string, int> parent_counts;

        // iterate records
        for (const auto &row : data)
        {
            if ((int)row.size() != N)
                continue; // skip malformed rows

            string this_val = row[node_idx];
            // if node value missing, skip this row for this node
            if (this_val == "?")
                continue;

            // build parent key (skip row if any parent is '?')
            string key = "";
            bool skip = false;
            for (auto &p : parents)
            {
                int pidx = BayesNet.get_index(p);
                if (pidx < 0 || pidx >= (int)row.size())
                {
                    skip = true;
                    break;
                }
                string pv = row[pidx];
                if (pv == "?")
                {
                    skip = true;
                    break;
                }
                key += pv;
                key.push_back('|');
            }
            if (skip)
                continue;

            counts[key][this_val]++;
            parent_counts[key]++;
        }

        vector<float> new_cpt;

        if (parents.empty())
        {
            // unconditional prior
            map<string, int> val_counts;
            int total = 0;
            for (const auto &row : data)
            {
                if ((int)row.size() <= node_idx)
                    continue;
                string v = row[node_idx];
                if (v == "?")
                    continue;
                val_counts[v]++;
                total++;
            }
            for (auto &val : node_vals)
            {
                float p = 0.0f;
                if (total == 0)
                    p = 1.0f / (float)node_vals.size();
                else
                    p = (float)val_counts[val] / (float)total;
                // round to 4 decimals
                p = round(p * 10000.0f) / 10000.0f;
                new_cpt.push_back(p);
            }
        }
        else
        {
            // conditional: iterate all possible parent combinations in lexicographic mixed-radix order
            int parent_combos = 1;
            vector<int> radices;
            for (auto &pv : parent_vals)
            {
                radices.push_back((int)pv.size());
                parent_combos *= max(1, (int)pv.size());
            }

            for (int comb = 0; comb < parent_combos; comb++)
            {
                // build key for this comb
                vector<int> digits(parents.size(), 0);
                int tmp = comb;
                for (int p = (int)parents.size() - 1; p >= 0; p--)
                {
                    int r = radices[p];
                    if (r > 0)
                    {
                        digits[p] = tmp % r;
                        tmp /= r;
                    }
                    else
                        digits[p] = 0;
                }

                string key = "";
                for (int p = 0; p < (int)parents.size(); p++)
                {
                    key += parent_vals[p][digits[p]];
                    key.push_back('|');
                }

                int denom = parent_counts.count(key) ? parent_counts[key] : 0;

                for (auto &val : node_vals)
                {
                    int numer = 0;
                    if (counts.count(key) && counts[key].count(val))
                        numer = counts[key][val];

                    float p = 0.0f;
                    if (denom == 0)
                    {
                        // no data for this parent combination -> assign uniform probability
                        p = 1.0f / (float)node_vals.size();
                    }
                    else
                    {
                        p = (float)numer / (float)denom;
                    }
                    p = round(p * 10000.0f) / 10000.0f;
                    new_cpt.push_back(p);
                }
            }
        }

        it->set_CPT(new_cpt);
    }

    cout << "Parameter learning completed." << endl;
}

// --------------------------- MAIN ---------------------------
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Hardcoded input files
    string bif_file = "hailfinder.bif";
    string data_file = "records.dat";

    // 1. Read network
    cout << "Reading network from: " << bif_file << endl;
    network BayesNet = read_network(bif_file.c_str());
    cout << "Network loaded. Nodes: " << BayesNet.netSize() << endl;

    // 2. Read dataset
    vector<vector<string>> data = read_dataset(data_file);
    cout << "Records loaded: " << data.size() << endl;

    // 3. Learn parameters (handle missing data as required)
    learn_parameters(BayesNet, data);

    // 4. Write solved network to required output file
    write_network("solved_hailfinder.bif", BayesNet);

    cout << "Parameter learning completed." << endl;

    return 0;
}


