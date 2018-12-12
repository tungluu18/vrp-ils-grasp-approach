#ifndef InitialSolution_H
#define InitialSolution_H

#include "Variables.h"
#include "Node.h"

#include <bits/stdc++.h>

using namespace std;

/**
 * @input is the list of customers, capacity 
 * @output: a permutation of customers that presents a route solution initially
 */
class InitialSolution {
private:    
    vector <Node> P;
    double density[1005][1005];

public:
    InitialSolution() {
        P.clear();
    }

    InitialSolution(const vector <Node> &_P) {
        P.clear();
        for (auto p: _P) P.push_back(p);
    }
    
    vector <Node> execute() {                
        double cost = 0.0;
        vector <Node> res;
        vector <int> chosen(N + 5, 0);  // whether a node has been visited
        res.push_back(P[0]);  // push depot into vector as the first element

        while (res.size() < N+1) {           
            // always start the tour from depot
            Node current_node = P[0];  
            int cap = 0;
            for(;;) {
                int x = find_nearest_node(current_node, chosen, CAPACITY - cap);
                if (x == -1) {
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

    int find_nearest_node(Node current_node, const vector <int> &chosen, int capAllowance) {
        double minDistance = 1e10;  
        int result = -1;
        
        // find the node has not been visited before and its demand is not upper limit
        for (int i = 1; i <= N; ++i) if (!chosen[i] && P[i].demand <= capAllowance) {        
            if (distance(current_node, P[i]) < minDistance) {
                // the node is nearer to the current_node
                minDistance = distance(current_node, P[i]);
                result = i;
            }
        }
        return result;
    }    
  
    void calculate_density() {
        const int p = 4;    // p has a big value, the customer with lower demand will be seleted
        const int k = 1;    // k has a big value, the customer whose distance is small will be selected        
        for (int i = 0; i < N+1; ++i) density[i][i] = 0.0;
        for (int j = 1; j <= N; ++j) {
            density[0][j] = pow(abs(CAPACITY - P[j].demand), p)
                /   (pow(distance(P[0], P[j]), k));
            density[j][0] = density[0][j];
        }
  
        for (int i = 1; i < N; ++i) {
            for (int j = i+1; j <= N; ++j) {
                density[i][j] = pow(abs(CAPACITY - P[i].demand - P[j].demand), p)
                    /   (pow(distance(P[i], P[j]), k) * density[0][i] * density[0][j]);
                density[j][i] = density[i][j];
            }
        }

        cerr << setprecision(10) << fixed << density[5][89] << endl;
    }

    vector <Node> run2() {        
        calculate_density();

        vector <int> indexRes;
        indexRes.push_back(0);
        vector <int> chosen(N+1, 0);  // whether a node has been visited
        chosen[0] = 1;
        for (int i = 0; i < N; ++i) {
            // find the next customer to visit is the customer that has 
            // the biggest density with the current customer
            int nextNode = -1;
            int curNode = indexRes.back();
            for (int j = 1; j <= N; ++j) if (!chosen[j]) {
                if (nextNode == -1 || density[curNode][nextNode] < density[curNode][j]) {
                    nextNode = j;
                }
            }
            indexRes.push_back(nextNode);
            chosen[nextNode] = 1;
        }

        vector <Node> res;
        for (auto i: indexRes) res.push_back(P[i]);
        return res;
    }
};

#endif