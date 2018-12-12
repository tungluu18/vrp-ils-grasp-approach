#ifndef LocalSearch_H
#define LocalSearch_H

#include "Variables.h"
#include "Node.h"

#include "Util.h" // splitting algorithm 

#include "IntraMovement.h"
#include "InterMovement.h"

using namespace std;

typedef vector <Node> vNode;
typedef vector <int> vi;

class LocalSearch {
private:
    vNode P;    

public:
    LocalSearch() {}

    LocalSearch(const vNode &_P) {P = _P;}

    pair <double, vi> execute() {
        pair <double, vi> result = splitting(P);
        // if (DEBUG) cerr << "Initially local search :" << result.first << endl;
        const int LOOP = 100;
        while (true) {
            pair <double, vi> new_result = result; 
            new_result = intra_movement::swap(new_result.first, new_result.second);
            new_result = intra_movement::move(new_result.first, new_result.second);
            new_result = intra_movement::inverse(new_result.first, new_result.second);
            new_result = inter_movement::exchange_1_0(new_result.first, new_result.second);
            new_result = inter_movement::exchange_1_0(new_result.first, new_result.second);
            if (result.first == new_result.first) break;
            result = new_result;
        }
        return result;
    }
};
#endif
