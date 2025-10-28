/********************************************************************
 *  starter2.cpp  –  Parameter learning for Bayesian networks
 ********************************************************************/

#include "bn_lib.h"
#include <random>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>   // for std::find

using namespace std;

/* ---------------------------------------------------------------
 *  const-correct wrapper for network::get_nth_node
 * --------------------------------------------------------------- */
inline auto get_nth_node(const network& net, int n)
    -> std::list<Graph_Node>::const_iterator
{
    auto mut = const_cast<network*>(&net)->get_nth_node(n);
    return std::list<Graph_Node>::const_iterator(mut);
}

/* ---------------------------------------------------------------
 *  Helper: index of a variable inside a factor (-1 if absent)
 * --------------------------------------------------------------- */
inline int find_var_index(const Factor& f, int var)
{
    for (size_t i = 0; i < f.variables.size(); ++i)
        if (f.variables[i] == var) return (int)i;
    return -1;
}

/* ---------------------------------------------------------------
 *  Factor constructor + VE functions (bn_lib.h only declares them)
 * --------------------------------------------------------------- */
inline Factor::Factor(int node_idx, network* net)
{
    auto it = net->get_nth_node(node_idx);
    const Graph_Node* node = &(*it);

    variables.push_back(node_idx);
    card.push_back(node->get_nvalues());

    for (const string& p : node->get_Parents()) {
        int pidx = net->get_index(p);
        if (pidx != -1) {
            variables.push_back(pidx);
            auto pit = net->get_nth_node(pidx);
            card.push_back(pit->get_nvalues());
        }
    }

    stride.resize(variables.size());
    stride[0] = 1;
    for (size_t i = 1; i < variables.size(); ++i)
        stride[i] = stride[i-1] * card[i-1];

    values = node->get_CPT();
}

inline Factor product(const Factor& a, const Factor& b) { /* ... as above ... */ }
inline Factor sum_out(const Factor& f, int var)         { /* ... as above ... */ }
inline Factor restrict_factor(const Factor& f, int var, int val) { /* ... as above ... */ }

/* ---------------------------------------------------------------
 *  Helper: parent configuration index
 * --------------------------------------------------------------- */
int get_parent_state_index(const Graph_Node& node,
                           const network& net,
                           const map<int,int>& evidence)
{
    const vector<string>& parents = node.get_Parents();
    if (parents.empty()) return 0;

    int idx = 0, stride = 1;
    for (const string& p_name : parents) {
        int p_idx = net.get_index(p_name);
        if (p_idx == -1) return -1;
        auto it = evidence.find(p_idx);
        int val = (it != evidence.end()) ? it->second : -1;
        if (val == -1)
            return -1;
        auto p_it = get_nth_node(net, p_idx);
        if (val >= p_it->get_nvalues()) return -1;
        idx += val * stride;
        stride *= p_it->get_nvalues();
    }
    return idx;
}

/* ---------------------------------------------------------------
 *  Variable Elimination – joint probability of evidence
 * --------------------------------------------------------------- */
Factor compute_joint(const network& net,
                     const map<int,int>& evidence,
                     const vector<int>& query_vars)
{
    vector<Factor> factors;
    for (int i = 0; i < net.netSize(); ++i)
        factors.emplace_back(i, const_cast<network*>(&net));

    /* incorporate evidence */
    for (const auto& ev : evidence) {
        for (auto& f : factors)
            if (find_var_index(f, ev.first) != -1)
                f = restrict_factor(f, ev.first, ev.second);
    }

    /* variables to eliminate */
    vector<int> to_elim;
    for (int i = 0; i < net.netSize(); ++i) {
        if (find(query_vars.begin(), query_vars.end(), i) == query_vars.end() &&
            evidence.find(i) == evidence.end())
            to_elim.push_back(i);
    }

    for (int v : to_elim) {
        vector<Factor> rel;
        for (auto it = factors.begin(); it != factors.end(); ) {
            if (find_var_index(*it, v) != -1) {
                rel.push_back(std::move(*it));
                it = factors.erase(it);
            } else ++it;
        }
        if (!rel.empty()) {
            Factor prod = rel[0];
            for (size_t k = 1; k < rel.size(); ++k)
                prod = product(prod, rel[k]);
            prod = sum_out(prod, v);
            factors.push_back(std::move(prod));
        }
    }

    if (factors.empty()) {
        Factor f; f.values = {1.0f}; return f;
    }
    Factor joint = factors[0];
    for (size_t i = 1; i < factors.size(); ++i)
        joint = product(joint, factors[i]);
    return joint;
}

float query_prob(const network& net, const map<int,int>& evidence)
{
    if (evidence.empty()) return 1.0f;
    Factor j = compute_joint(net, evidence, {});
    float Z = 0.0f;
    for (float v : j.values) Z += v;
    return Z > 0 ? Z : 1e-30f;
}

/* ---------------------------------------------------------------
 *  BIF parser (robust for Hailfinder)
 * --------------------------------------------------------------- */
network read_network(const char* filename)
{
    network net;
    ifstream f(filename);
    if (!f) { cerr << "Cannot open BIF\n"; return net; }

    string line, cur_var;
    map<string,int> name2idx;

    while (getline(f,line)) {
        line = trim(line);
        if (line.empty() || line.substr(0,2)=="//") continue;

        stringstream ss(line);
        string token; ss >> token;

        if (token == "variable") {
            ss >> cur_var;
            cur_var = trim(cur_var);
            if (cur_var.front()=='"' && cur_var.back()=='"')
                cur_var = cur_var.substr(1,cur_var.size()-2);
            net.addNode(Graph_Node(cur_var,0,{}));
            name2idx[cur_var] = net.netSize()-1;
        }
        else if (token == "type" && !cur_var.empty()) {
            string dummy; int card;
            ss >> dummy >> dummy >> card;               // discrete [ card ]
            vector<string> vals;
            size_t b = line.find("{");
            if (b != string::npos) {
                string vstr = line.substr(b+1);
                vstr = vstr.substr(0, vstr.find("}"));
                stringstream vss(vstr);
                string v;
                while (getline(vss,v,',')) vals.push_back(trim(v));
            }
            auto it = net.get_nth_node(net.netSize()-1);
            it->set_values(vals);
        }
        else if (token == "probability") {
            size_t op = line.find("("), cp = line.find(")");
            if (op==string::npos || cp==string::npos) continue;
            string cond = line.substr(op+1, cp-op-1);
            cond = trim(cond);
            stringstream css(cond);
            string child; getline(css,child,'|');
            child = trim(child);
            vector<string> parents;
            string pstr; if (getline(css,pstr)) {
                pstr = trim(pstr);
                stringstream pss(pstr);
                string p;
                while (getline(pss,p,',')) parents.push_back(trim(p));
            }
            int cidx = name2idx[child];
            auto cit = net.get_nth_node(cidx);
            cit->set_Parents(parents);

            /* read the table line */
            while (getline(f,line)) {
                line = trim(line);
                if (line.find("table") != string::npos) {
                    string t = line.substr(line.find("table")+5);
                    t = trim(t);
                    if (!t.empty() && t.back()==';') t.pop_back();
                    stringstream tss(t);
                    vector<float> cpt;
                    string num;
                    while (getline(tss,num,',')) {
                        num = trim(num);
                        if (num.empty()) continue;
                        float v = (num=="-1" || num=="-1.0") ? -1.0f : stof(num);
                        cpt.push_back(v);
                    }
                    cit->set_CPT(cpt);
                    break;
                }
            }
        }
    }

    /* build children links */
    for (int i = 0; i < net.netSize(); ++i) {
        auto n = net.get_nth_node(i);
        for (const string& p : n->get_Parents()) {
            int pidx = name2idx[p];
            auto pit = net.get_nth_node(pidx);
            pit->add_child(i);
        }
    }
    cout << "BIF parsed – " << net.netSize() << " nodes.\n";
    return net;
}

/* ---------------------------------------------------------------
 *  Write network (exact format required by the assignment)
 * --------------------------------------------------------------- */
void write_network(const char* filename, network& net)
{
    ofstream out(filename);
    out << "// Bayesian Network\n\n";

    for (int i = 0; i < net.netSize(); ++i) {
        auto n = net.get_nth_node(i);
        out << "variable " << n->get_name() << " {\n";
        out << "  type discrete [ " << n->get_nvalues() << " ] = { ";
        const vector<string>& vals = n->get_values();
        for (size_t j = 0; j < vals.size(); ++j)
            out << vals[j] << (j+1<vals.size() ? ", " : "");
        out << " };\n}\n\n";
    }

    for (int i = 0; i < net.netSize(); ++i) {
        auto n = net.get_nth_node(i);
        out << "probability ( " << n->get_name();
        const vector<string>& par = n->get_Parents();
        if (!par.empty()) {
            out << " | ";
            for (size_t j = 0; j < par.size(); ++j)
                out << par[j] << (j+1<par.size() ? ", " : "");
        }
        out << " ) {\n    table ";
        const vector<float>& cpt = n->get_CPT();
        out << fixed << setprecision(6);
        for (size_t j = 0; j < cpt.size(); ++j)
            out << cpt[j] << (j+1<cpt.size() ? ", " : "");
        out << ";\n}\n\n";
    }
    cout << "Network written to " << filename << "\n";
}

/* ---------------------------------------------------------------
 *  Dataset reader
 * --------------------------------------------------------------- */
vector<vector<string>> read_dataset(const string& fn)
{
    vector<vector<string>> data;
    ifstream f(fn);
    if (!f) { cerr << "Cannot open dataset\n"; return data; }

    string line;
    while (getline(f,line)) {
        line = trim(line);
        if (line.empty()) continue;
        stringstream ss(line);
        vector<string> rec;
        string tok;
        while (getline(ss,tok,',')) {
            tok = trim(tok);
            if (tok.front()=='"' && tok.back()=='"')
                tok = tok.substr(1,tok.size()-2);
            rec.push_back(tok);
        }
        if (!rec.empty()) data.push_back(move(rec));
    }
    return data;
}

/* ---------------------------------------------------------------
 *  Random CPT initialization (only for -1 entries)
 * --------------------------------------------------------------- */
void initialize_random_cpts(network& net)
{
    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis(0.01,0.99);

    for (int i = 0; i < net.netSize(); ++i) {
        auto it = net.get_nth_node(i);
        vector<float> cpt = it->get_CPT();
        int nval = it->get_nvalues();

        for (size_t row = 0; row < cpt.size(); row += nval) {
            float known_sum = 0.0f;
            int unk = 0;
            for (int k = 0; k < nval; ++k)
                if (cpt[row+k] == -1.0f) ++unk;
                else known_sum += cpt[row+k];

            if (unk == 0) continue;
            float remain = max(0.0f, 1.0f - known_sum);
            float rsum = 0.0f;
            vector<float> rnd(unk,0.0f);
            int ui = 0;
            for (int k = 0; k < nval; ++k)
                if (cpt[row+k] == -1.0f) {
                    rnd[ui] = (float)dis(gen);
                    rsum += rnd[ui++];
                }

            ui = 0;
            for (int k = 0; k < nval; ++k)
                if (cpt[row+k] == -1.0f) {
                    float p = rsum>0 ? remain * rnd[ui]/rsum : remain/unk;
                    cpt[row+k] = max(0.0001f, min(0.9999f, p));
                    ++ui;
                }

            /* final normalisation */
            float sum = 0.0f;
            for (int k = 0; k < nval; ++k) sum += cpt[row+k];
            if (sum>0) {
                float norm = 1.0f/sum;
                for (int k = 0; k < nval; ++k) {
                    cpt[row+k] *= norm;
                    cpt[row+k] = round(cpt[row+k]*10000.0f)/10000.0f;
                }
            }
        }
        it->set_CPT(cpt);
    }
    cout << "Random CPTs initialised.\n";
}

/* ---------------------------------------------------------------
 *  Log-likelihood (used for convergence)
 * --------------------------------------------------------------- */
float calculate_log_likelihood(const network& net,
                               const vector<vector<string>>& data)
{
    float ll = 0.0f;
    int N = net.netSize();
    for (const auto& row : data) {
        map<int,int> ev;
        for (int i = 0; i < N; ++i)
            if (row[i] != "?") {
                auto it = get_nth_node(net,i);
                int v = get_value_index(it->get_values(), row[i]);
                if (v != -1) ev[i] = v;
            }
        float p = query_prob(net, ev);
        ll += (p > 1e-30f) ? log(p) : -30.0f;
    }
    return ll;
}

/* ---------------------------------------------------------------
 *  Gibbs sampling – returns expected counts for missing vars
 * --------------------------------------------------------------- */
int sample_from(const vector<float>& probs, mt19937& gen)
{
    float sum = 0.0f;
    for (float p : probs) sum += p;
    if (sum <= 0.0f) return 0;
    uniform_real_distribution<> d(0.0, sum);
    float r = d(gen), cum = 0.0f;
    for (size_t i = 0; i < probs.size(); ++i) {
        cum += probs[i];
        if (r <= cum) return (int)i;
    }
    return (int)probs.size()-1;
}

float local_conditional_prob(const network& net,
                             int node_idx, int val_k,
                             const map<int,int>& state)
{
    auto it = get_nth_node(net, node_idx);
    const Graph_Node& node = *it;
    const vector<string>& parents = node.get_Parents();
    const vector<float>& cpt = node.get_CPT();
    int nval = node.get_nvalues();

    int pa_idx = 0, stride = 1;
    for (const string& p : parents) {
        int pidx = net.get_index(p);
        auto sit = state.find(pidx);
        if (sit == state.end()) return 0.0f;
        int pv = sit->second;
        auto pit = get_nth_node(net, pidx);
        if (pv >= pit->get_nvalues()) return 0.0f;
        pa_idx += pv * stride;
        stride *= pit->get_nvalues();
    }
    return cpt[pa_idx * nval + val_k];
}

/* joint counts for each missing variable */
vector<vector<vector<float>>> gibbs_joint_counts(
        const network& net,
        const map<int,int>& evidence,
        const vector<int>& missing,
        int samples = 2000, int burn = 500)
{
    int N = net.netSize();
    mt19937 gen(random_device{}());

    map<int,int> state = evidence;
    for (int v : missing) {
        auto it = get_nth_node(net,v);
        uniform_int_distribution<> d(0, it->get_nvalues()-1);
        state[v] = d(gen);
    }

    /* containers for counts */
    vector<vector<vector<float>>> counts(N);
    for (int v : missing) {
        auto it = get_nth_node(net,v);
        int nval = it->get_nvalues();
        int npa  = it->get_CPT().size() / nval;
        counts[v].assign(npa, vector<float>(nval, 0.0f));
    }

    /* burn-in */
    for (int s = 0; s < burn; ++s) {
        for (int v : missing) {
            auto it = get_nth_node(net,v);
            int nval = it->get_nvalues();
            vector<float> probs(nval,0.0f);
            for (int k = 0; k < nval; ++k) {
                map<int,int> tmp = state; tmp[v] = k;
                float pr = 1.0f;
                for (int n = 0; n < N; ++n) {
                    auto nit = get_nth_node(net,n);
                    const vector<string>& pars = nit->get_Parents();
                    bool ok = true;
                    for (const string& p : pars)
                        if (tmp.find(net.get_index(p)) == tmp.end()) { ok=false; break; }
                    if (!ok) continue;
                    pr *= local_conditional_prob(net, n, tmp[n], tmp);
                }
                probs[k] = pr;
            }
            state[v] = sample_from(probs, gen);
        }
    }

    /* sampling */
    for (int s = 0; s < samples; ++s) {
        for (int v : missing) {
            auto it = get_nth_node(net,v);
            int nval = it->get_nvalues();
            vector<float> probs(nval,0.0f);
            for (int k = 0; k < nval; ++k) {
                map<int,int> tmp = state; tmp[v] = k;
                float pr = 1.0f;
                for (int n = 0; n < N; ++n) {
                    auto nit = get_nth_node(net,n);
                    const vector<string>& pars = nit->get_Parents();
                    bool ok = true;
                    for (const string& p : pars)
                        if (tmp.find(net.get_index(p)) == tmp.end()) { ok=false; break; }
                    if (!ok) continue;
                    pr *= local_conditional_prob(net, n, tmp[n], tmp);
                }
                probs[k] = pr;
            }
            state[v] = sample_from(probs, gen);
        }
        /* record counts */
        for (int v : missing) {
            int pa = get_parent_state_index(*get_nth_node(net,v), net, state);
            if (pa >= 0) counts[v][pa][state[v]] += 1.0f;
        }
    }
    return counts;
}

/* ---------------------------------------------------------------
 *  One EM iteration
 * --------------------------------------------------------------- */
void perform_em_iteration(network& net,
                          const vector<vector<string>>& data)
{
    int N = net.netSize();

    /* expected counts: [node][parent_config][value] */
    vector<vector<vector<float>>> exp_cnt(N);
    vector<vector<float>> exp_par(N);   // denominator

    for (int i = 0; i < N; ++i) {
        auto it = get_nth_node(net,i);
        int nval = it->get_nvalues();
        int npa  = it->get_CPT().size() / nval;
        exp_cnt[i].assign(npa, vector<float>(nval,0.0f));
        exp_par[i].assign(npa, 0.0f);
    }

    for (const auto& row : data) {
        map<int,int> obs;
        vector<int> miss;
        for (int i = 0; i < N; ++i) {
            if (row[i] != "?") {
                auto it = get_nth_node(net,i);
                int v = get_value_index(it->get_values(), row[i]);
                if (v != -1) obs[i] = v;
            } else miss.push_back(i);
        }

        if (miss.empty()) {               /* complete record */
            for (int i = 0; i < N; ++i) {
                auto it = get_nth_node(net,i);
                int pa = get_parent_state_index(*it, net, obs);
                int val = obs[i];
                if (pa >= 0 && val >= 0) {
                    exp_cnt[i][pa][val] += 1.0f;
                    exp_par[i][pa]      += 1.0f;
                }
            }
            continue;
        }

        /* incomplete – Gibbs */
        auto counts = gibbs_joint_counts(net, obs, miss, 2000, 500);
        for (int i : miss) {
            for (size_t pa = 0; pa < exp_cnt[i].size(); ++pa) {
                for (int v = 0; v < (int)exp_cnt[i][pa].size(); ++v) {
                    float c = counts[i][pa][v];
                    float p = c / 2000.0f;
                    exp_cnt[i][pa][v] += p;
                    exp_par[i][pa]    += p;
                }
            }
        }
    }

    /* M-step */
    for (int i = 0; i < N; ++i) {
        auto it = net.get_nth_node(i);
        vector<float> old_cpt = it->get_CPT();
        vector<float> new_cpt = old_cpt;
        int nval = it->get_nvalues();

        for (size_t pa = 0; pa < exp_par[i].size(); ++pa) {
            int row = (int)pa * nval;
            float denom = exp_par[i][pa];
            if (denom < 1e-8f) denom = 1.0f;

            for (int k = 0; k < nval; ++k) {
                if (old_cpt[row+k] == -1.0f) {
                    float p = exp_cnt[i][pa][k] / denom;
                    new_cpt[row+k] = round(p*10000.0f)/10000.0f;
                }
            }

            /* normalise row */
            float sum = 0.0f;
            for (int k = 0; k < nval; ++k) sum += new_cpt[row+k];
            if (sum > 0) {
                float norm = 1.0f/sum;
                for (int k = 0; k < nval; ++k) {
                    new_cpt[row+k] *= norm;
                    new_cpt[row+k] = round(new_cpt[row+k]*10000.0f)/10000.0f;
                }
            }
        }
        it->set_CPT(new_cpt);
    }
    cout << "EM iteration done.\n";
}

/* ---------------------------------------------------------------
 *  Multi-start EM driver
 * --------------------------------------------------------------- */
void learn_parameters(network& net,
                      const vector<vector<string>>& data)
{
    const int NUM_STARTS = 5;
    const int MAX_ITERS   = 50;
    const float CONV_TH   = 1e-4f;

    network best = net;
    float best_ll = -1e30f;

    cout << "\n=== Multi-start EM (" << NUM_STARTS << " runs) ===\n";

    for (int r = 0; r < NUM_STARTS; ++r) {
        cout << "\n--- Run " << (r+1) << " ---\n";
        network cur = read_network("hailfinder.bif");
        if (cur.netSize()==0) continue;

        initialize_random_cpts(cur);
        float prev = calculate_log_likelihood(cur, data);
        cout << "Init LL = " << fixed << setprecision(4) << prev << "\n";

        bool converged = false;
        for (int it = 0; it < MAX_ITERS; ++it) {
            perform_em_iteration(cur, data);
            float now = calculate_log_likelihood(cur, data);
            cout << "Iter " << (it+1) << " LL = " << now << "\n";

            if (fabs(now - prev) < CONV_TH) {
                cout << "Converged.\n"; converged = true; break;
            }
            if (now < prev - 1e-3f) {
                cout << "LL dropped – stopping run.\n"; break;
            }
            prev = now;
        }
        float final = calculate_log_likelihood(cur, data);
        if (final > best_ll) {
            best_ll = final; best = cur;
            cout << "NEW BEST LL = " << best_ll << "\n";
        }
    }

    net = best;
    cout << "\n=== EM finished – best LL = "
         << fixed << setprecision(4) << best_ll << " ===\n";
}

/* ---------------------------------------------------------------
 *  main
 * --------------------------------------------------------------- */
int main(int argc, char* argv[])
{
    ios::sync_with_stdio(false); cin.tie(nullptr);

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <bif> <data>\n";
        return 1;
    }

    string bif = argv[1], dat = argv[2];

    cout << "Loading network...\n";
    network net = read_network(bif.c_str());
    if (net.netSize()==0) return 1;

    cout << "Loading data...\n";
    auto data = read_dataset(dat);
    cout << data.size() << " records.\n";

    learn_parameters(net, data);

    write_network("solved_hailfinder.bif", net);
    cout << "All done!\n";
    return 0;
}
