#include <iostream>
#include <string>
#include <vector>
#include "bn_lib.h"

using namespace std;

// --------------------------- MAIN ---------------------------

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string bif_file, data_file;

    if (argc >= 3)
    {
        // Read file names from command-line arguments
        bif_file = argv[1];
        data_file = argv[2];
    }
    else
    {
        // Default values (fallback)
        bif_file = "hailfinder.bif";
        data_file = "records.dat";
    }

    // 1. Read network
    cout << "Reading network from: " << bif_file << endl;
    network BayesNet = read_network(bif_file.c_str());
    cout << "Network loaded. Nodes: " << BayesNet.netSize() << endl;

    // 2. Read dataset
    vector<vector<string>> data = read_dataset(data_file);
    cout << "Records loaded: " << data.size() << endl;

    // 3. Learn parameters using EM (handles missing values)
    int em_iterations = 50; // You can adjust as needed
    cout << "Learning parameters using EM (" << em_iterations << " iterations)..." << endl;
    learn_parameters_EM(BayesNet, data, em_iterations);

    // 4. Write solved network to file
    string solved_file = "solved.bif";
    write_network(solved_file.c_str(), BayesNet);
    cout << "EM-based parameter learning completed." << endl;
    cout << "Solved network written to: " << solved_file << endl;

    return 0;
}
