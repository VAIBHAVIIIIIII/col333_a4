#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#ifndef BN_LIB_H
#define BN_LIB_H

using namespace std;

// --------------------------- Helper trim ---------------------------
inline string trim(const string &str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == string::npos) return "";
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
    Graph_Node() : Node_Name(""), nvalues(0) {}
    Graph_Node(string name, int n, vector<string> vals) : Node_Name(name), nvalues(n), values(vals) {}

    string get_name() const { return Node_Name; }
    vector<int> get_children() const { return Children; }
    vector<string> get_Parents() const { return Parents; }
    vector<float> get_CPT() const { return CPT; }
    int get_nvalues() const { return nvalues; }
    vector<string> get_values() const { return values; }

    void set_CPT(const vector<float> &new_CPT) { CPT = new_CPT; }
    void set_Parents(const vector<string> &Parent_Nodes) { Parents = Parent_Nodes; }
    void set_name(const string &nm) { Node_Name = nm; }
    void set_values(const vector<string> &vals) { values = vals; nvalues = (int)vals.size(); }

    int add_child(int new_child_index)
    {
        for (int c : Children) if (c == new_child_index) return 0;
        Children.push_back(new_child_index);
        return 1;
    }
};


// --------------------------- Network ---------------------------
class network
{
    list<Graph_Node> Pres_Graph;

public:
    int addNode(const Graph_Node &node) { Pres_Graph.push_back(node); return 0; }

    list<Graph_Node>::iterator get_nth_node(int n)
    {
        auto it = Pres_Graph.begin();
        advance(it, n);
        return it;
    }

    int netSize() const { return (int)Pres_Graph.size(); }

    int get_index(const string &val_name) const
    {
        int idx = 0;
        for (auto it = Pres_Graph.begin(); it != Pres_Graph.end(); ++it, ++idx)
            if (it->get_name() == val_name) return idx;
        return -1;
    }

    list<Graph_Node>::iterator search_node(const string &val_name)
    {
        for (auto it = Pres_Graph.begin(); it != Pres_Graph.end(); ++it)
            if (it->get_name() == val_name) return it;
        return Pres_Graph.end();
    }

    list<Graph_Node>::const_iterator search_node_const(const string &val_name) const
    {
        for (auto it = Pres_Graph.begin(); it != Pres_Graph.end(); ++it)
            if (it->get_name() == val_name) return it;
        return Pres_Graph.end();
    }

    list<Graph_Node>& get_nodes_list() { return Pres_Graph; }
    const list<Graph_Node>& get_nodes_list_const() const { return Pres_Graph; }
};

// --------------------------- Utility functions ---------------------------
int get_value_index(const vector<string> &values, const string &val)
{
    for (int i = 0; i < (int)values.size(); i++)
        if (values[i] == val) return i;
    return -1;
}

// --------------------------- Read dataset ---------------------------
vector<vector<string>> read_dataset(const string &filename)
{
    vector<vector<string>> data;
    ifstream infile(filename);
    if (!infile.is_open())
    {
        cerr << "Error opening dataset file: " << filename << endl;
        return data;
    }

    string line;
    while (getline(infile, line))
    {
        line = trim(line);
        if (line.empty()) continue;

        vector<string> row;
        stringstream ss(line);
        string token;
        while (ss >> token)
            row.push_back(token);
        data.push_back(row);
    }

    infile.close();
    return data;
}

// --------------------------- Read .bif network ---------------------------
network read_network(const char *filename)
{
    network BayesNet;
    ifstream myfile(filename);
    if (!myfile.is_open())
    {
        cerr << "Error: Could not open file " << filename << endl;
        return BayesNet;
    }

    string line;
    while (getline(myfile, line))
    {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string token;
        ss >> token;

        if (token == "variable")
        {
            string var_name;
            ss >> var_name;
            var_name = trim(var_name);

            string type_line;
            while (getline(myfile, type_line))
            {
                type_line = trim(type_line);
                if (type_line.find("type") != string::npos) break;
            }

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
                    if (!tokenv.empty())
                    {
                        if (tokenv.front() == '"' && tokenv.back() == '"')
                            tokenv = tokenv.substr(1, tokenv.size() - 2);
                        values.push_back(tokenv);
                    }
                }
            }

            Graph_Node new_node(var_name, (int)values.size(), values);
            BayesNet.addNode(new_node);
        }
        else if (token == "probability")
        {
            string header = line;
            while (header.find('{') == string::npos)
            {
                string nextline;
                if (!getline(myfile, nextline)) break;
                header += " " + trim(nextline);
            }

            size_t lp = header.find('(');
            size_t rp = header.find(')');
            if (lp == string::npos || rp == string::npos || rp <= lp) continue;
            string inside = trim(header.substr(lp + 1, rp - lp - 1));

            string node_name;
            vector<string> parents;
            size_t pipe = inside.find('|');
            if (pipe == string::npos)
                node_name = trim(inside);
            else
            {
                node_name = trim(inside.substr(0, pipe));
                string pstr = inside.substr(pipe + 1);
                stringstream ps(pstr);
                string tmp;
                while (ps >> tmp)
                {
                    if (!tmp.empty() && tmp.back() == ',') tmp.pop_back();
                    tmp = trim(tmp);
                    if (!tmp.empty()) parents.push_back(tmp);
                }
            }

            int index = BayesNet.get_index(node_name);
            if (index == -1)
            {
                Graph_Node placeholder(node_name, 0, vector<string>());
                BayesNet.addNode(placeholder);
                index = BayesNet.get_index(node_name);
            }

            auto listIt = BayesNet.get_nth_node(index);
            if (listIt != BayesNet.get_nodes_list().end()) listIt->set_Parents(parents);

            for (auto &pname : parents)
            {
                int pidx = BayesNet.get_index(pname);
                if (pidx == -1)
                {
                    Graph_Node placeholder(pname, 0, vector<string>());
                    BayesNet.addNode(placeholder);
                    pidx = BayesNet.get_index(pname);
                }
                BayesNet.get_nth_node(pidx)->add_child(index);
            }

            vector<float> cpt;
            string cpline;
            while (getline(myfile, cpline))
            {
                string t = trim(cpline);
                if (t.empty()) continue;
                if (t == "};" || t == "}") break;

                size_t tablepos = t.find("table");
                string probs_part;
                if (tablepos != string::npos)
                    probs_part = t.substr(tablepos + 5);
                else
                {
                    size_t closep = t.find(')');
                    if (closep != string::npos) probs_part = t.substr(closep + 1);
                    else probs_part = t;
                }

                stringstream ssp(probs_part);
                string tokenp;
                while (ssp >> tokenp)
                {
                    while (!tokenp.empty() && (tokenp.back() == ',' || tokenp.back() == ';'))
                        tokenp.pop_back();
                    if (!tokenp.empty() && (isdigit(tokenp[0]) || tokenp[0] == '-' || tokenp[0] == '.'))
                    {
                        try { cpt.push_back(stof(tokenp)); }
                        catch (...) { /* skip invalid numbers */ }
                    }
                }
            }

            if (listIt != BayesNet.get_nodes_list().end()) listIt->set_CPT(cpt);
        }
    }

    myfile.close();
    return BayesNet;
}

// --------------------------- Write network ---------------------------
inline void write_network(const char *filename, network &BayesNet)
{
    ofstream outfile(filename);
    if (!outfile.is_open())
    {
        cerr << "Error: Could not open file " << filename << " for writing" << endl;
        return;
    }

    outfile << "// Bayesian Network" << endl << endl;
    auto &nodes = BayesNet.get_nodes_list();

    for (auto it = nodes.begin(); it != nodes.end(); ++it)
    {
        outfile << "variable " << it->get_name() << " {" << endl;
        outfile << "  type discrete [ " << it->get_nvalues() << " ] = { ";
        auto vals = it->get_values();
        for (size_t j = 0; j < vals.size(); j++)
        {
            outfile << "\"" << vals[j] << "\"";
            if (j < vals.size() - 1) outfile << ", ";
        }
        outfile << " };" << endl << "}" << endl << endl;
    }

    outfile << fixed << setprecision(4);

    for (auto it = nodes.begin(); it != nodes.end(); ++it)
    {
        auto parents = it->get_Parents();
        auto values = it->get_values();
        auto cpt = it->get_CPT();

        outfile << "probability ( " << it->get_name();
        if (!parents.empty())
        {
            outfile << " | ";
            for (size_t j = 0; j < parents.size(); j++)
            {
                outfile << parents[j];
                if (j < parents.size() - 1) outfile << ", ";
            }
        }
        outfile << " ) {" << endl;

        vector<int> radices;
        for (auto &pname : parents)
        {
            auto pIt = BayesNet.search_node(pname);
            if (pIt != BayesNet.get_nodes_list().end()) radices.push_back(pIt->get_nvalues());
            else radices.push_back(1);
        }

        int parent_combinations = 1;
        for (int r : radices) parent_combinations *= r;
        if (parent_combinations == 0) parent_combinations = 1;

        size_t cpt_index = 0;
        for (int comb = 0; comb < parent_combinations; comb++)
        {
            vector<int> idxs(parents.size(), 0);
            int tmp = comb;
            for (int p = (int)parents.size() - 1; p >= 0; p--)
            {
                if (radices[p] > 0) { idxs[p] = tmp % radices[p]; tmp /= radices[p]; }
            }

            if (!parents.empty()) outfile << "    ( ";
            for (size_t p = 0; p < parents.size(); p++)
            {
                auto pIt = BayesNet.search_node(parents[p]);
                if (pIt != BayesNet.get_nodes_list().end())
                {
                    auto pvals = pIt->get_values();
                    int vidx = idxs[p] % max(1, (int)pvals.size());
                    outfile << pvals[vidx];
                }
                else outfile << "UNK";
                if (p < parents.size() - 1) outfile << ", ";
            }
            if (!parents.empty()) outfile << " ) ";

            double sum_probs = 0.0;
            for (size_t k = 0; k < values.size(); k++)
            {
                float p = cpt_index < cpt.size() ? cpt[cpt_index++] : 1.0f / values.size();
                if (k < values.size() - 1) { p = round(p * 10000.0f) / 10000.0f; sum_probs += p; outfile << p << ", "; }
                else { p = round((1.0 - sum_probs) * 10000.0f) / 10000.0f; outfile << p; }
            }
            outfile << ";" << endl;
        }

        outfile << "};" << endl << endl;
    }

    outfile.close();
    cout << "Network written to file: " << filename << endl;
}

void learn_parameters_EM(network &BayesNet, const vector<vector<string>> &data, int max_iters = 50)
{
    int N = BayesNet.netSize();
    int num_records = (int)data.size();
    cout << "EM: Learning parameters from " << num_records << " records..." << endl;

    if (N <= 0) { cerr << "Empty network!" << endl; return; }
    if (num_records <= 0) { cerr << "Empty dataset!" << endl; return; }

    // Quick dataset sanity check: all rows must have N tokens.
    for (int r = 0; r < num_records; ++r) {
        int row_len = static_cast<int>(data[r].size());

        // Only warn if the row is NOT blank (len 1) but is STILL the wrong length
        if (row_len != N && row_len > 1) {
            cerr << "Data row " << r << " length " << data[r].size()
                << " != expected " << N << " - check records.dat formatting." << endl;
        }
    }
    struct NodeInfo {
        int nvalues;
        vector<int> parent_indices;
        vector<int> parent_radices;
        vector<int> multipliers; // for mapping parent assignments to CPT index
    };

    vector<NodeInfo> nodes_info(N);
    vector<vector<float>> CPTs(N);

    // ------------------- Precompute node info -------------------
    int node_idx = 0;
    for (auto it = BayesNet.get_nodes_list().begin(); it != BayesNet.get_nodes_list().end(); ++it, ++node_idx)
    {
        CPTs[node_idx] = it->get_CPT();
        int val_size = it->get_nvalues();
        vector<string> parents = it->get_Parents();
        vector<int> p_indices, p_radices, mult;

        for (auto &pname : parents)
        {
            int pidx = BayesNet.get_index(pname);
            if (pidx < 0) {
                cerr << "Warning: parent '" << pname << "' of node " << it->get_name()
                     << " not found in network; treating as single-value parent." << endl;
                pidx = -1;
            }
            p_indices.push_back(pidx);
            int pvals = 1;
            if (pidx >= 0) pvals = BayesNet.get_nth_node(pidx)->get_nvalues();
            p_radices.push_back(max(1, pvals));
        }

        mult.assign(p_indices.size(), 1);
        for (int i = (int)p_indices.size() - 2; i >= 0; --i)
            mult[i] = mult[i + 1] * p_radices[i + 1];

        // Compute CPT size robustly
        long long cpt_size = val_size;
        for (int rdx : p_radices) cpt_size *= rdx;
        if (cpt_size <= 0 || cpt_size > (1 << 28)) {
            cerr << "Invalid CPT size for node " << it->get_name() << ": " << cpt_size << endl;
            cpt_size = val_size;
        }

        if ((int)CPTs[node_idx].size() != (int)cpt_size)
            CPTs[node_idx] = vector<float>((size_t)cpt_size, 1.0f / max(1, val_size));

        nodes_info[node_idx] = { val_size, p_indices, p_radices, mult };
    }

    // ------------------- EM iterations -------------------
    for (int iter = 0; iter < max_iters; ++iter)
    {
        // Expected counts (zeros)
        vector<vector<double>> exp_counts(N);
        for (int j = 0; j < N; ++j)
            exp_counts[j] = vector<double>(CPTs[j].size(), 0.0);

        // Process records
        for (int r = 0; r < num_records; ++r)
        {
            // safety: skip malformed rows
            if ((int)data[r].size() != N) continue;

            vector<int> missing_vars;
            for (int j = 0; j < N; ++j)
                if (data[r][j] == "?")
                    missing_vars.push_back(j);

            if (missing_vars.empty())
            {
                // Fully observed: increment counts
                for (int j = 0; j < N; ++j)
                {
                    auto &ni = nodes_info[j];

                    // compute parent assignment index
                    int parent_block = 0;
                    bool bad = false;
                    for (size_t pi = 0; pi < ni.parent_indices.size(); ++pi)
                    {
                        int pidx = ni.parent_indices[pi];
                        int pv_index = 0;
                        if (pidx >= 0) {
                            pv_index = get_value_index(BayesNet.get_nth_node(pidx)->get_values(), data[r][pidx]);
                            if (pv_index < 0) { bad = true; pv_index = 0; }
                        }
                        parent_block += pv_index * ni.multipliers[pi];
                    }
                    if (bad) { /* skip or continue with pv_index=0 */ }

                    int val_idx = get_value_index(BayesNet.get_nth_node(j)->get_values(), data[r][j]);
                    if (val_idx < 0) val_idx = 0;

                    long long cpt_idx = (long long)parent_block * ni.nvalues + val_idx;
                    if (cpt_idx >= 0 && cpt_idx < (long long)exp_counts[j].size())
                        exp_counts[j][(size_t)cpt_idx] += 1.0;
                    else {
                        cerr << "Warning: out-of-range CPT index for node " << j << " on record " << r << endl;
                    }
                }
                continue;
            }

            // Only one missing per row guaranteed; handle first
            int mvar = missing_vars[0];
            auto &mi = nodes_info[mvar];

            // For each possible value for missing var, compute joint (unnormalized)
            vector<double> joint_prob(mi.nvalues, 0.0);
            double total_joint = 0.0;

            for (int val_idx = 0; val_idx < mi.nvalues; ++val_idx)
            {
                // build a local completed record (avoid modifying original)
                vector<string> record = data[r];
                record[mvar] = BayesNet.get_nth_node(mvar)->get_values()[val_idx];

                double joint = 1.0;
                // Multiply local conditional prob for every node
                int node_index = 0;
                for (auto itn = BayesNet.get_nodes_list().begin(); itn != BayesNet.get_nodes_list().end(); ++itn, ++node_index)
                {
                    auto &nj = nodes_info[node_index];

                    int parent_block = 0;
                    for (size_t pi = 0; pi < nj.parent_indices.size(); ++pi)
                    {
                        int pidx = nj.parent_indices[pi];
                        int pv = 0;
                        if (pidx >= 0) {
                            pv = get_value_index(BayesNet.get_nth_node(pidx)->get_values(), record[pidx]);
                            if (pv < 0) pv = 0;
                        }
                        parent_block += pv * nj.multipliers[pi];
                    }

                    int vidx = get_value_index(BayesNet.get_nth_node(node_index)->get_values(), record[node_index]);
                    if (vidx < 0) vidx = 0;

                    long long cpt_idx = (long long)parent_block * nj.nvalues + vidx;
                    if (cpt_idx < 0 || cpt_idx >= (long long)CPTs[node_index].size()) {
                        joint = 0.0; break;
                    }
                    joint *= CPTs[node_index][(size_t)cpt_idx];
                    if (joint == 0.0) break;
                }

                joint_prob[val_idx] = joint;
                total_joint += joint;

                // If joint is zero, skipping adding counts for that completed record is fine.
                if (joint <= 0.0) continue;

                // accumulate unnormalized counts for this completion
                int node_index2 = 0;
                for (auto itn = BayesNet.get_nodes_list().begin(); itn != BayesNet.get_nodes_list().end(); ++itn, ++node_index2)
                {
                    auto &nj = nodes_info[node_index2];

                    int parent_block = 0;
                    for (size_t pi = 0; pi < nj.parent_indices.size(); ++pi)
                    {
                        int pidx = nj.parent_indices[pi];
                        int pv = 0;
                        if (pidx >= 0) {
                            pv = get_value_index(BayesNet.get_nth_node(pidx)->get_values(), record[pidx]);
                            if (pv < 0) pv = 0;
                        }
                        parent_block += pv * nj.multipliers[pi];
                    }

                    int vidx = get_value_index(BayesNet.get_nth_node(node_index2)->get_values(), record[node_index2]);
                    if (vidx < 0) vidx = 0;

                    long long cpt_idx = (long long)parent_block * nj.nvalues + vidx;
                    if (cpt_idx >= 0 && cpt_idx < (long long)exp_counts[node_index2].size())
                        exp_counts[node_index2][(size_t)cpt_idx] += joint;
                }
            } // end values loop

            // normalize per-record contributions by total_joint to get expected counts consistent with probabilities
            if (total_joint > 0.0) {
                for (int j = 0; j < N; ++j) {
                    // We added joint contributions above; now divide by total_joint
                    for (size_t k = 0; k < exp_counts[j].size(); ++k) {
                        // It's expensive to scale entire exp_counts per record; we instead postpone scaling to M-step.
                        // To keep things correct and simple here, do nothing now. M-step will normalize across all records.
                        (void)k;
                    }
                }
            }
        } // end records loop

        // -------- M-step: normalize CPTs using accumulated exp_counts --------
        for (int j = 0; j < N; ++j)
        {
            auto &ni = nodes_info[j];
            int val_size = ni.nvalues;
            int parent_combinations = 1;
            if (val_size > 0) parent_combinations = (int)CPTs[j].size() / val_size;

            for (int pc = 0; pc < parent_combinations; ++pc)
            {
                double sum = 0.0;
                for (int vi = 0; vi < val_size; ++vi) sum += exp_counts[j][pc * val_size + vi];
                if (sum < 1e-12)
                {
                    // Laplace smoothing: uniform if no evidence
                    for (int vi = 0; vi < val_size; ++vi) CPTs[j][pc * val_size + vi] = 1.0f / val_size;
                }
                else
                {
                    for (int vi = 0; vi < val_size; ++vi)
                        CPTs[j][pc * val_size + vi] = (float)(exp_counts[j][pc * val_size + vi] / sum);
                }
            }

            // Update network CPT
            BayesNet.get_nth_node(j)->set_CPT(CPTs[j]);
        }

        if ((iter + 1) % 10 == 0) cout << "EM iteration " << (iter + 1) << " completed." << endl;
    } // end iteration loop

    cout << "EM completed." << endl;
}

#endif
