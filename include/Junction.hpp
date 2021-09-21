#ifndef _JUNCTION_H_
#define _JUNCTION_H_

#include <iostream>

#include "Segment.hpp"
#include "Edge.hpp"
#include "Weight.hpp"

class Segment;

class Edge;

using namespace std;

class Junction {
protected:
    char mSourceDir;
    char mTargetDir;

    double mCredibility;

    bool mIsInferred;
    bool mHasLowerBoundLimit;     // used in ILP processing, if true then lower bound is 1, otherwise 0

    Weight *mWeight;
    Segment *mSource;
    Segment *mTarget;

    Edge *mEdgeA;         //
    Edge *mEdgeB;

public:
    Junction(Segment *aSource, Segment *aTarget, char aSourceDir, char aTargetDir, double aCoverage,
             double aCredibility, double aCopy, bool aInferred, bool aIsBounded, bool aIsSourceSinkJunction);

    ~Junction();

    vector<string> getInfo();

    double getCredibility();

    double setCredibility(double aCredibility);

    bool isInferred();

    bool hasLowerBoundLimit();

    bool hasCopy();

    void setInferred();

    void resetInferred();

    void setHasLowerBoundLimit();

    void resetHasLowerBoundLimit();

    void checkLowerBound();

    void restoreCopy();

    void backupCopy();

    Weight *getWeight();

    Segment *getSource();

    Segment *getTarget();

    Edge *getEdgeA();

    Edge *getEdgeB();

    char getSourceDir();
    char getTargetDir();

    // functionality
    void insertEdgesToVertices();

    void print();
};

#endif
