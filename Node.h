#ifndef Node_H
#define Node_H

#include <bits/stdc++.h>

#define sqr(a) ((a) * (a))

class Node {
public:
    int x, y, demand, id;
    Node() {}
    Node(int _x, int _y, int _demand) : x(_x), y(_y), demand(_demand) {}    
};

double distance(Node P, Node Q) {
    return sqrt(sqr(P.x - Q.x) + sqr(P.y - Q.y));
}

#endif