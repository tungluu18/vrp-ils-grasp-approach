#include <bits/stdc++.h>
#include <fstream>
#include <experimental/filesystem>

#include "control_IO.h"

#define sqr(a) ((a) * (a))
#define EPS 0.000001

#define DEBUG 0

using namespace std;
namespace filesys = std::experimental::filesystem;

struct Node
{
    int x, y, demand, id;
    Node() {}
    Node(int _x, int _y, int _demand) : x(_x), y(_y), demand(_demand) {}
};

typedef vector < int > vi;
typedef vector < Node > vNode;
// Global Variables
int N;
int CAPACITY;
vector < Node > X;

bool isNumber(string s)
{
    for (int i = 0; i < s.size(); ++i) if (s[i] < '0' || s[i] > '9') return false;
    return true;
}

vi getInt(stringstream lineStream)
{
    vi result;
    string cell;
    while (lineStream >> cell) if (isNumber(cell)) result.push_back(atoi(cell.c_str()));
    return result;
}

void readFile(const string filepath, vector<Node> &X)
{
    ifstream fi(filepath.c_str());

    N = -1;
    string line, headline;
    vector <int> result;
    
    for (;;)
    {
        getline(fi, line, '\n');        

        stringstream lineStream(line);        
        lineStream >> headline;
        
        if (headline == "EOF") break;
        if (headline == "DIMENSION")
        {            
            result = getInt(stringstream(line));
            N = result[0];            
            X.resize(N + 10);
        }
        if (headline == "CAPACITY")
        {            
            result = getInt(stringstream(line));
            CAPACITY = result[0];
        }        
        if (headline == "NODE_COORD_SECTION")        
        {               
            for (int i = 0; i < N; ++i) 
            {
                getline(fi, line, '\n');
                stringstream ss(line);  
                ss >> X[i+1].id;
                ss >> X[i+1].x >> X[i+1].y;
            }                     
        }
        if (headline == "DEMAND_SECTION")
        {               
            for (int i = 0; i < N; ++i) 
            {
                getline(fi, line, '\n');    
                stringstream ss(line);  
                ss >> X[i+1].id;
                ss >> X[i+1].demand;
            }
        }
        // if (headline == "DEPOT_SECTION")
        // {
        //     getline(fi, line, '\n');    
        //     stringstream ss(line);  
        //     ss >> X[0].x;

        //     getline(fi, line, '\n');    
        //     ss << line;  
        //     ss >> X[0].y;
        // }
    }

    fi.close();
}

void outFile(double min_cost, vi rout, string fileout)
{   
    // for (auto r: rout) cerr << r << "\n";
    ofstream of(fileout);
        
    of << setprecision(0) << fixed << min_cost << endl;

    // calculate_split number of trip, it is the number of 0 - depot on the path minus 1
    int n_trip = -1;
    for (auto r: rout) n_trip += (r == 0);
    of << n_trip;

    // print the trips
    int trip_th = 0;  // the order of trip
    for (int i = 0; i+1 < rout.size(); ++i)
    {
        if (!rout[i]) // depot - start a new trip
        {
            trip_th ++;
            of << "\nRout #" << trip_th << ":";            
        } 
        else 
        {
            of << " " << rout[i];
        }
    }
 
    of.close();
}

double distance(Node P, Node Q) 
{
    return sqrt(sqr(P.x - Q.x) + sqr(P.y - Q.y));
}

double distance_index(int i, int j)
{
    return distance(X[i], X[j]);
}

int find_nearest_node(Node current_node, const vector <Node> &P, const vector <int> &chosen, int capAllowance)
{       
    double minDistance = 1e10;  
    int result = -1;
    
    // find the node has not been visited before and its demand is not upper limit
    for (int i = 1; i <= N; ++i) if (!chosen[i] && P[i].demand <= capAllowance) 
    {        
        if (distance(current_node, P[i]) < minDistance) 
        {
            // the node is nearer to the current_node
            minDistance = distance(current_node, P[i]);
            result = i;
        }
    }
    return result;
}

vector < Node > initialize_solution(const vector < Node > &P) 
{        
    double cost = 0.0;
    vector <Node> res;
    vector <int> chosen(N + 5, 0);  // whether a node has been visited
    res.push_back(P[0]);  // push depot into vector as the first element

    while (res.size() < N+1)
    {           
        // always start the tour from depot
        Node current_node = P[0];  
        int cap = 0;
        for(;;)
        {
            int x = find_nearest_node(current_node, P, chosen, CAPACITY - cap);
            if (x == -1) 
            {
                // cannot find new node fulfilling capacity allowance
                cost += distance(current_node, P[0]);  // go back to depot from current_node
                break;
            }
            // next node to go is P[x]
            res.push_back(P[x]);        
            chosen[x] = 1;      
            cap += P[x].demand;
            cost += distance(current_node, P[x]);  // go from current_node to P[x]
            current_node = P[x];
        }
    }

    return res;
}

vNode convert_rout_to_nodes(const vi rout)
{
    vNode result(1, X[0]);
    for (auto r: rout) if (r) result.push_back(X[r]);
    return result;
}

pair < double, vi > calculate_split(const vector < Node > &P) 
{
    vector < double > dp(N+5);
    vector < double > trace(N+5);
    vi rout;
    dp[0] = 0.0;
    for (int i = 1; i <= N; ++i) 
    {
        dp[i] = 1e10;
        int cap = 0;
        double dist = 0.0;
        // if (DEBUG) cerr << "Debugging on " << i << "\n";
        for (int j = i; j > 0; --j) 
        {
            cap += P[j].demand;
            if (j < i) dist += distance(P[j], P[j+1]);          
            if (cap > CAPACITY) break;
            if (dist + distance(P[0], P[i]) + distance(P[0], P[j]) + dp[j-1] < dp[i]) 
            {
                dp[i] = dp[j-1] + distance(P[0], P[i]) + distance(P[0], P[j]) + dist;
                trace[i] = j;
                // if (DEBUG) cerr << cap << "\n";
            }
        }
    }
    
    // tracing    
    rout.push_back(0);
    if (DEBUG) cerr << "Capacity debug on each rout:\n";
    for (int u = N; u > 0; u = trace[u] - 1) 
    {                    
        for (int i = trace[u]; i <= u; ++i) rout.push_back(P[i].id);        
        rout.push_back(0);
        // if (DEBUG) cerr << trace[u] << " " << u << "\n";
        // if (DEBUG) cerr << cap << endl;
    }

    return {dp[N], rout};
}

vi calculate_cap_rout(const vi &rout) 
{
    vi result(N+1);
    int start = 0, cap = 0;
    for (int i = 1; i <= N; ++i)
    {
        if (!rout[i])
        {
            for (int j = start+1; j < i; ++j) result[j] = cap;
            cap = 0;
            start = i;
        }
        else 
            cap += X[rout[i]].demand;
    }
    return result;
}

pair < double, vi > swap_movement(double current_min_cost, const vi &rout)
{    
    double new_min_cost = current_min_cost;
    vi in_capacity = calculate_cap_rout(rout);
    int u = -1, v = -1;
    
    for (int i = 0; i < rout.size(); ++i) if (rout[i]) 
    {
        for (int j = i+1; j < rout.size(); ++j) if (rout[j]) 
        {
            // try to swap rout[i] and rout[j]
            double new_cost = current_min_cost;

            if (i+1 == j) 
            {
                new_cost -= distance_index(rout[i-1], rout[i]) + distance_index(rout[i+1], rout[i+2]);
                new_cost += distance_index(rout[i-1], rout[i+1]) + distance_index(rout[i], rout[i+2]);
            }
            else 
            {
                // check constraint on capacity
                if (in_capacity[i] - X[rout[i]].demand + X[rout[j]].demand > CAPACITY) continue;
                if (in_capacity[j] - X[rout[j]].demand + X[rout[i]].demand > CAPACITY) continue;

                // remove rout[i] and rout[j] from routing path
                new_cost -= distance_index(rout[i-1], rout[i]) + distance_index(rout[i], rout[i+1]);
                new_cost -= distance_index(rout[j-1], rout[j]) + distance_index(rout[j], rout[j+1]);
                
                // add rout[j] between rout[i-1] and rout[i+1]
                new_cost += distance_index(rout[i-1], rout[j]) + distance_index(rout[j], rout[i+1]);

                // add rout[i] between rout[j-1] and rout[j+1]
                new_cost += distance_index(rout[j-1], rout[i]) + distance_index(rout[i], rout[j+1]);
            }
            
            if (new_cost < new_min_cost) 
            {
                new_min_cost = new_cost;
                u = i, v = j;
            }
        }
    }

    if (u == -1) // cannot find a pair to swap that reduces the cost    
        return {current_min_cost, rout};    

    vi new_rout = rout;  
    swap(new_rout[u], new_rout[v]);
    vNode P = convert_rout_to_nodes(new_rout);
    return calculate_split(P);
}

pair < double, vi > shift_movement(double current_min_cost, const vi &rout)
{    
    double new_min_cost = current_min_cost;
    vi in_capacity = calculate_cap_rout(rout);
    int u = -1, v = -1;
    
    for (int i = 0; i < rout.size(); ++i) if (rout[i]) 
    {
        for (int j = 0; j+1 < rout.size(); ++j) if (j+1 != i) 
        {
            // check constraint on capacity
            if (in_capacity[j] + X[rout[i]].demand > CAPACITY) continue;

            // try to relocate rout[i] after rout[j]
            double new_cost = current_min_cost;
            
            // remove rout[i]
            new_cost -= distance_index(rout[i-1], rout[i]) + distance_index(rout[i], rout[i+1]);

            // insert rout[i] between rout[j] and rout[j+1]
            new_cost -= distance_index(rout[j], rout[j+1]);
            new_cost += distance_index(rout[j], rout[i]) + distance_index(rout[i], rout[j+1]);
            
            if (new_cost < new_min_cost)
            {
                new_min_cost = new_cost;
                u = i, v = j;
            }
        }
    }

    if (u == -1)
        return {current_min_cost, rout};

    vi new_rout = rout;    
    new_rout.insert(new_rout.begin() + v, rout[u]);  // insert rout[u] after rout[v]
    if (v < u) 
        new_rout.erase(new_rout.begin() + u + 1);
    else 
        new_rout.erase(new_rout.begin() + u);
    vNode P = convert_rout_to_nodes(new_rout);
    return calculate_split(P);
}

pair < double, vi > _2_opt_movement(double current_min_cost, const vi &rout)
{
    double new_min_cost = current_min_cost;
    vi in_capacity = calculate_cap_rout(rout);
    int u = -1, v = -1;

    for (int i = 1; i < rout.size(); ++i) if (rout[i] && rout[i+1])
    {
        for (int j = i+2; j < rout.size(); ++j) if (rout[j] && rout[j+1]) 
        {
            // check constraint on capacity
            if (in_capacity[i] - X[rout[i+1]].demand + X[rout[j+1]].demand > CAPACITY) continue;
            if (in_capacity[j] - X[rout[j+1]].demand + X[rout[i+1]].demand > CAPACITY) continue;

            double new_cost = current_min_cost;
            new_cost -= distance_index(rout[i], rout[i+1]) + distance_index(rout[i+1], rout[i+2]);
            new_cost -= distance_index(rout[j], rout[j+1]) + distance_index(rout[j+1], rout[j+2]);

            new_cost += distance_index(rout[i], rout[j+1]) + distance_index(rout[j+1], rout[i+2]);
            new_cost += distance_index(rout[j], rout[i+1]) + distance_index(rout[i+1], rout[j+2]);            

            if (new_cost < new_min_cost)
            {
                new_min_cost = new_cost;
                u = i, v = j;
            }
        }
    }

    if (u == -1)
        return {current_min_cost, rout};
    
    vi new_rout = rout;
    swap(new_rout[u+1], new_rout[v+1]);
    vNode P = convert_rout_to_nodes(new_rout);
    return calculate_split(P);
}

pair < double, vi > local_search(vector < Node > P)
{                    
    pair < double, vi > result = calculate_split(P);
    for (int loop = 0; loop < 200; ++loop)
    {
        pair < double, vi > new_result;
        new_result = shift_movement(result.first, result.second);
        new_result = swap_movement(new_result.first, new_result.second);
        new_result = shift_movement(new_result.first, new_result.second);
        new_result = swap_movement(new_result.first, new_result.second);
        new_result = _2_opt_movement(new_result.first, new_result.second);

        P = convert_rout_to_nodes(new_result.second);        
        result = calculate_split(P);        
    }        
    return result;
}

pair < double, vi > iterated_local_search()
{
    vector < Node > P = initialize_solution(X);
    pair < double, vi > result = local_search(P);

    srand(time(NULL));
    const int ITERATE_TIME = 100;
    for (int iter = 0; iter < ITERATE_TIME; iter++)
    {
        const int NUMBER_OF_SWAPS = 500;

        // kick solution by finding its subsequence' s next permutations
        int index[N+1];
        for (int i = 0; i <= N; ++i) index[i] = P[i].id;
        for (int i = 0; i < NUMBER_OF_SWAPS; ++i)
        {
            int u = rand() % N + 1;
            int v = rand() % N + 1;
            if (u > v) swap(u, v);
            next_permutation(index + u, index + v);
        }

        // next permutation by index
        for (int i = 0; i <= N; ++i) P[i] = X[index[i]];

        pair < double, vi > new_result = local_search(P);
        P = convert_rout_to_nodes(new_result.second);
        cerr << result.first << " " << new_result.first << "\n";
        if (new_result.first < result.first) result = new_result;        
    }

    return result;
}

string verify(double min_cost, vi rout)
{
    double total_cost = 0;
    int cap = 0;
    int previous = 0;
    for (auto &r: rout)
    {
        total_cost += distance_index(previous, r);
        if (!r)
        {
            cerr << "capacity is: " << cap << endl;
            if (cap > CAPACITY) return "False\n";
            cap = 0;        
        } 
        else 
        {
            cap += X[r].demand;
        }        
        previous = r;
    }
    cerr << total_cost << endl;
    if (fabs(total_cost-min_cost) < 1) 
        return "True\n";
    else 
        return "False\n";
}

int main(int argc, char **argv)
{
    string filename = "X-n153-k22.vrp";
    if (argc > 1) filename = argv[1];

    // fast read write on file
    ios_base::sync_with_stdio(false);
    
    // get directory containing .vrp files
    const string fileDir = string(filesys::current_path()) + "/X/";
    
    // for all file in a directory
    for (auto &p : filesys::directory_iterator(fileDir))
    {
        // cerr << p << endl;        
    }

    // testing on reading a .vrp file
    X.clear();
    readFile(fileDir + filename, X);

    X[0] = Node(1, -1, 0);  // X[0] is the depot

    // solution obtained by ils algorithm                
    pair < double, vi > answer = iterated_local_search();    
    cerr << answer.first << endl;
    
    // cerr << verify(answer.first, answer.second);
    outFile(answer.first, answer.second, fileDir + "MyAnswer.sol");
    return 0;
}
