#ifndef InterMovement_H
#define InterMovement_H

#include <bits/stdc++.h>

#include "Variables.h"
#include "Node.h"
#include "Util.h"
using namespace std;

typedef vector <Node> vNode;
typedef vector <int> vi;

#define  dist(_i, _j) distance(X[_i], X[_j])
/**
 * @output: the demand on each route corresponding to each customer
 */
vi demand_on_rout(const vi &rout) {
    vi result(rout.size());
    int start = 0, cap = 0;
    for (int i = 1; i < rout.size(); ++i) {
        if (!rout[i]) {
            for (int j = start+1; j < i; ++j) result[j] = cap;
            cap = 0;
            start = i;
        } else  cap += X[rout[i]].demand;        
    }
    return result;
}

namespace inter_movement {    
    pair <double, vi> exchange_1_0(const double current_min_cost, const vi &rout) {
        double new_min_cost = current_min_cost;
        vi demanding = demand_on_rout(rout);

        // (u -> v) is the optimal movement
        int u = -1, v = -1;

        for (int i = 0; i < rout.size(); ++i) if (rout[i]) {            
            // move rout[i] to the left of rout[j]
            for (int j = 1; j < rout.size(); ++j) {
                if (i == j || i+1 == j) continue;
                if (rout[j]) {
                    if (demanding[j] + X[rout[i]].demand > CAPACITY) continue;
                } else  {
                    if (demanding[j-1] + X[rout[i]].demand > CAPACITY) continue;
                }
                double new_cost = current_min_cost;
                new_cost -= dist(rout[i-1], rout[i]) + dist(rout[i], rout[i+1]);
                new_cost -= dist(rout[j-1], rout[j]);
                new_cost += dist(rout[i-1], rout[i+1]);
                new_cost += dist(rout[j-1], rout[i]) + dist(rout[i], rout[j]);

                if (new_cost < new_min_cost) {
                    new_min_cost = new_cost;
                    u = i, v = j;
                }
            }
        }

        if (u == -1) return {current_min_cost, rout};

        vi new_rout = rout;
        new_rout.insert(new_rout.begin()+v, rout[u]);  // insert rout[v] to the left of rout[v]
        if (u < v)
            new_rout.erase(new_rout.begin()+u);
        else new_rout.erase(new_rout.begin()+u+1);
        
        // remove duplicate elements
        vector <int>::iterator it = unique(new_rout.begin(), new_rout.end());
        new_rout.resize(distance(new_rout.begin(), it));

        return {new_min_cost, new_rout};
    }
};

#endif