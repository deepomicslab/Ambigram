#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <vector>

#include "Edge.hpp"
#include "Segment.hpp"
#include "Weight.hpp"

using namespace std;

class Segment;

class Edge;

class Weight;

class Vertex {
protected:
    int mId;                    // vertex id
    char mDir;                  // +/-

    double mCredibility;        // credibility, used in linear programming

    bool mIsOrphan;              // record whether originally orphan
    bool mHasCheckedOrphan;
    bool mIsVisited;

    Edge *mShortestPrevEdge;

    Segment *mSegment;         // the host segment
    Weight *mWeight;           // weight object including copy number

    vector<Edge *> *mEdgesAsSource;    // edges for this vertex as source
    vector<Edge *> *mEdgesAsTarget;    // edges for this vertex as target

public:
    // constructor and destructor
    Vertex(int aId, char aDir, Weight *aWeight, double aCredibility);

    ~Vertex();

    // getter and setter
    int getId();

    int getStart();

    int getEnd();

    char getDir();

    string getInfo();

    inline bool operator==(const Vertex &other) const {
        if (mId == other.mId && mDir == other.mDir) return true;
        return false;
    };

    double getInCoverage();

    double getOutCoverage();

    double getCredibility();

    void setCredibility(double aCredibility);

    bool isVisited();

    bool hasCopy();

    bool hasEdgeAsSource();

    bool hasEdgeAsTarget();

    bool isOrphan(bool aIsOriginal = true);

    void checkOrphan();          // using only time RIGHT AFTER READING GRAPH
    void setVisited();

    void resetVisited();

    // traversal action and recover (by 1 copy number)
    void traverse();

    void recover();

    Weight *getWeight();

    Segment *getSegment();

    Edge *getShortestPrevEdge();

    void setSegment(Segment *aSeg);

    void setShortestPrevEdge(Edge *aEdge);

    vector<Edge *> *getEdgesAsSource();

    vector<Edge *> *getEdgesAsTarget();

    void insertEdgeAsSource(Edge *aEdge);

    void insertEdgeAsTarget(Edge *aEdge);

    Edge *findEdgeAsSource(Vertex *aTargetVertex);

    Vertex *getComplementVertex();

    // print func
    void print();

    bool isReverse(Vertex *aVertex);
};

#endif
