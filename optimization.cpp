#include <bits/stdc++.h>
#include <fstream>
#include <experimental/filesystem>

#include "Variables.h"
#include "Node.h"
#include "InitialSolution.h"
#include "LocalSearch.h"
#include "Util.h"

#define sqr(a) ((a) * (a))
#define EPS 0.000001

using namespace std;
namespace filesys = std::experimental::filesystem;

typedef vector <int> vi;
typedef vector <Node> vNode;

bool isNumber(string s) {
    for (int i = 0; i < s.size(); ++i) if (s[i] < '0' || s[i] > '9') return false;
    return true;
}

vi getInt(stringstream lineStream) {
    vi result;
    string cell;
    while (lineStream >> cell) if (isNumber(cell)) result.push_back(atoi(cell.c_str()));
    return result;
}

void readFile(const string filepath, vector<Node> &X) {
    ifstream fi(filepath.c_str());

    N = -1;
    string line, headline;
    vector <int> result;
    
    for (;;) {
        getline(fi, line, '\n');        

        stringstream lineStream(line);        
        lineStream >> headline;
        
        if (headline == "EOF") break;
        if (headline == "DIMENSION") {            
            result = getInt(stringstream(line));
            N = result[0];            
            X.resize(N + 10);
        }
        if (headline == "CAPACITY") {            
            result = getInt(stringstream(line));
            CAPACITY = result[0];
        }        
        if (headline == "NODE_COORD_SECTION") {               
            for (int i = 0; i < N; ++i) {
                getline(fi, line, '\n');
                stringstream ss(line);  
                ss >> X[i+1].id;
                ss >> X[i+1].x >> X[i+1].y;
            }                     
        }
        if (headline == "DEMAND_SECTION") {               
            for (int i = 0; i < N; ++i) {
                getline(fi, line, '\n');    
                stringstream ss(line);  
                ss >> X[i+1].id;
                ss >> X[i+1].demand;
            }
        }
    }

    fi.close();
}

void outFile(double min_cost, vi rout, string fileout) {   
    cerr << fileout << endl;
    // for (auto r: rout) cerr << r << "\n";
    ofstream of(fileout);
        
    of << setprecision(0) << fixed << min_cost << endl;

    // calculate_split number of trip, it is the number of 0 - depot on the path minus 1
    int n_trip = -1;
    for (auto r: rout) n_trip += (r == 0);
    of << n_trip;

    // print the trips
    int trip_th = 0;  // the order of trip
    for (int i = 0; i+1 < rout.size(); ++i) {
        if (!rout[i]) { // depot - start a new trip 
            trip_th ++;
            of << "\nRout #" << trip_th << ":";            
        } 
        else {
            of << " " << rout[i];
        }
    }
 
    of.close();
}

vNode convert_rout_to_nodes(const vi rout) {
    vNode result(1, X[0]);
    for (auto r: rout) if (r) result.push_back(X[r]);
    return result;
}

pair <double, vi> iterated_local_search() {    
    vector <Node> P = (new InitialSolution(X))->execute();
    
    pair <double, vi> result = (new LocalSearch(P))->execute();

    for (int iter = 0; iter < ITERATE_TIME; iter++) {
        const int NUMBER_OF_SWAPS = 1000;

        // kick solution by finding its subsequence' s next permutations
        int index[N+1];
        for (int i = 0; i <= N; ++i) index[i] = P[i].id;
        for (int i = 0; i < NUMBER_OF_SWAPS; ++i) {
            int u = rand() % N + 1;
            int v = rand() % N + 1;
            if (u > v) swap(u, v);
            next_permutation(index + u, index + v);
        }

        // next permutation by index
        for (int i = 0; i <= N; ++i) P[i] = X[index[i]];

        pair <double, vi> new_result = (new LocalSearch(P))->execute();        
        P = convert_rout_to_nodes(new_result.second);

        if (DEBUG) cerr << "old -> new: " << result.first << " " << new_result.first << endl;
        if (new_result.first < result.first) result = new_result;        
    }

    // result = splitting(P);
    return result;
}

int main(int argc, char **argv) {
    srand(time(NULL));

    string filename = "X-n153-k22.vrp";
    DEBUG = atoi(argv[1]);
    if (argc > 2) filename = argv[2];
    cerr << filename << endl;
    // fast read write on file
    ios_base::sync_with_stdio(false);
    
    // get directory containing .vrp files
    const string fileDir = string(filesys::current_path()) + "/X/";
    
    // for all file in a directory
    for (auto &p : filesys::directory_iterator(fileDir)) {}

    // testing on reading a .vrp file
    X.clear();
    readFile(fileDir + filename, X);

    X[0] = Node(1, -1, 0);  // X[0] is the depot

    // solution obtained by ils algorithm                
    pair <double, vi> answer = iterated_local_search();    
    cerr << "Answert is: " << answer.first << endl;

    // cerr << "Verify answer:\n" << verify(answer.first, answer.second);
    outFile(answer.first, answer.second, fileDir + "MyAnswer.sol");
    
    cerr << "file: " << filename << endl;    
    return 0;
}
