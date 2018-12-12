#ifndef IntraMovement_H
#define IntraMovement_H

#include <bits/stdc++.h>

#include "Variables.h"
#include "Node.h"
#include "Util.h"
using namespace std;

typedef vector <Node> vNode;
typedef vector <int> vi;

#define  dist(_i, _j) distance(X[_i], X[_j])

namespace intra_movement {
    pair <double, vi> inverse(const double current_min_cost, const vi &rout) {
        double new_min_cost = current_min_cost;

        // segment [u,v] is optimal to inverse
        int u = -1, v = -1;
        for (int i = 1; i < rout.size(); ++i) if (rout[i]) {
            for (int j = i+1; j < rout.size(); ++j) {
                if (!rout[j]) break;
                double new_cost = current_min_cost;
                new_cost -= dist(rout[i-1], rout[i]) + dist(rout[j], rout[j+1]);
                new_cost += dist(rout[i-1], rout[j]) + dist(rout[i], rout[j+1]);
                if (new_cost < new_min_cost) {
                    new_min_cost = new_cost;
                    u = i, v = j;
                }
            }
        }
        
        if (u == -1) return {current_min_cost, rout};
        vi new_rout = rout;
        reverse(new_rout.begin()+u, new_rout.begin()+v+1);
        return {new_min_cost, new_rout};
    }

    pair <double, vi> move(const double current_min_cost, const vi &rout) {
        double new_min_cost = current_min_cost;

        // (u -> v) is the optimal movement
        int u = -1, v = -1;

        for (int i = 0; i < rout.size(); ++i) if (rout[i]) {
            int i_left, i_right;
            for (i_left = i; rout[i_left]; --i_left);
            for (i_right = i; rout[i_right]; ++i_right);

            // move rout[i] to the left of rout[j]
            for (int j = i_left+2; j < i_right; ++j) {
                if (i == j || i+1 == j) continue;
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

        if (u == -1) 
            return {current_min_cost, rout};
        // if (DEBUG) cerr << u << " " << v << endl;
        vi new_rout = rout;
        new_rout.insert(new_rout.begin()+v, rout[u]);  // insert rout[v] to the left of rout[v]
        if (u < v)
            new_rout.erase(new_rout.begin()+u);
        else 
            new_rout.erase(new_rout.begin()+u+1);
        return {new_min_cost, new_rout};
    }

    pair <double, vi> swap(const double current_min_cost, const vi &rout) {
        double new_min_cost = current_min_cost;

        // (u, v) is the optimal pair to swap
        int u = -1, v = -1;
        
        for (int i = 0; i < rout.size(); ++i) if (rout[i]) {
            for (int j = i+1; j < rout.size(); ++j) {                
                if (!rout[j]) break;   // only consider other customers in the same route
                
                // swap rout[i] and rout[j]
                double new_cost = current_min_cost;

                if (i+1 == j) {
                    new_cost -= dist(rout[i-1], rout[i]) + dist(rout[i+1], rout[i+2]);
                    new_cost += dist(rout[i-1], rout[i+1]) + dist(rout[i], rout[i+2]);
                } else {
                    new_cost -= dist(rout[i-1], rout[i]) + dist(rout[i], rout[i+1]);
                    new_cost -= dist(rout[j-1], rout[j]) + dist(rout[j], rout[j+1]);
                    new_cost += dist(rout[i-1], rout[j]) + dist(rout[j], rout[i+1]);
                    new_cost += dist(rout[j-1], rout[i]) + dist(rout[i], rout[j+1]);
                }
                
                if (new_cost < new_min_cost) {
                    new_min_cost = new_cost;
                    u = i, v = j;
                }
            }
        }

        if (u == -1) // cannot find a pair to swap that reduces the cost    
            return {current_min_cost, rout};    

        vi new_rout = rout;  
        std::swap(new_rout[u], new_rout[v]);
        return {new_min_cost, new_rout};
    }
};

#endif