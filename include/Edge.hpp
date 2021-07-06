#ifndef _EDGE_H_
#define _EDGE_H_

#include "Vertex.hpp"
#include "Junction.hpp"

using namespace std;

class Vertex;

class Weight;

class Junction;

class Edge {
protected:
    // int mId;    // edge id

    double mCredibility;

    bool mIsVisited;
    // bool mIsPreferred;

    Weight *mWeight;
    Vertex *mSource;
    Vertex *mTarget;
    Junction *mJunction;

public:
    // constructor and destructor
    Edge(Vertex *aSource, Vertex *aTarget, Weight *aWeight, double aCredibility);

    ~Edge();

    // getter and setter
    // int id();
    string getInfo();

    double getCredibility();

    void setCredibility(double aCredibility);

    bool isVisited();

    // bool isPreferred();
    void setVisited();

    void resetVisited();
    // void setPreferred();

    Weight *getWeight();

    Vertex *getSource();

    Vertex *getTarget();

    Junction *getJunction();

    void setSource(Vertex *aSource);

    void setTarget(Vertex *aTarget);

    void setJunction(Junction *aJunction);

    // functionality
    bool hasCopy();

    bool doesConnectSameDir();

    void recover(double copy = 1);

    void traverse();

    Edge *getComplementEdge();

    // print func
    void print();
};

#endif
