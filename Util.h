#ifndef Util_H
#define Util_H

#include "Variables.h"
#include "Node.h"

#include <bits/stdc++.h>

using namespace std;

typedef vector <Node> vNode;
typedef vector <int> vi;

/**
 * @input: a permutation of customers
 * @output:
 *  -- first: the min cost to visit all customer follow that order
 *  -- second: the routes
 */
pair <double, vi> splitting(const vector < Node > &P) {
    vector < double > dp(N+5);
    vector < double > trace(N+5);
    vi rout;
    dp[0] = 0.0;
    for (int i = 1; i <= N; ++i) {
        dp[i] = 1e10;
        int cap = 0;
        double dist = 0.0;
        // if (DEBUG) cerr << "Debugging on " << i << "\n";
        for (int j = i; j > 0; --j) {
            cap += P[j].demand;
            if (j < i) dist += distance(P[j], P[j+1]);          
            if (cap > CAPACITY) break;
            if (dist + distance(P[0], P[i]) + distance(P[0], P[j]) + dp[j-1] < dp[i]) {
                dp[i] = dp[j-1] + distance(P[0], P[i]) + distance(P[0], P[j]) + dist;
                trace[i] = j;
                // if (DEBUG) cerr << cap << "\n";
            }
        }
    }
    
    // tracing routes   
    rout.push_back(0);    
    for (int u = N; u > 0; u = trace[u]-1) {                    
        for (int i = trace[u]; i <= u; ++i) rout.push_back(P[i].id);
        rout.push_back(0);        
    }

    return {dp[N], rout};
}

#define  dist(_i, _j) distance(X[_i], X[_j])

string verify(double min_cost, vi rout) {
    double total_cost = 0;
    int cap = 0;
    int previous = 0;
    for (auto &r: rout) {
        total_cost += dist(previous, r);
        if (!r) {
            cerr << "capacity is: " << cap << endl;
            if (cap > CAPACITY) return "False\n";
            cap = 0;        
        } else {
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

#endif