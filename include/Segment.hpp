#ifndef _SEGMENT_H_
#define _SEGMENT_H_

#include "Vertex.hpp"
#include "Weight.hpp"

using namespace std;

class Vertex;

class Weight;

class Segment {
protected:
    int mId;    // segment id
    string mChrom;
    int mStart;
    int mEnd;
//        seg partition, on of source id
    int partition;
public:
    int getPartition() const;

    void setPartition(int partition);

protected:
    double mCredibility;

    bool mHasLowerBoundLimit;  // used in ILP processing, if true then lower bound is 1, otherwise 0
    bool mIsOrphan;            // record whether segment is orphan originally
    bool mHasCheckedOrphan;

    Weight *mWeight;
    Vertex *mPositiveVertex;
    Vertex *mNegativeVertex;

public:
    // constructor and destructor
    Segment(int aId, string aChrom, int aStart, int aEnd, double aCoverage, double aCredibility, double aCopy);

    ~Segment();

    // getter and setter
    int getId();

    string getChrom();

    int getStart();

    int getEnd();

    double getCredibility();

    void setCredibility(double aCredibility);

    bool isOrphan(bool aIsOriginal = true);

    bool isDeadEnd();

    bool hasCopy();

    bool hasLowerBoundLimit();

    void checkOrphan();

    void checkLowerBound();

    void setHasLowerBoundLimit();

    void resetHasLowerBoundLimit();

    void restoreCopy();

    void backupCopy();

    Weight *getWeight();

    Vertex *getPositiveVertex();

    Vertex *getNegativeVertex();

    void setPositiveVertex(Vertex *aVertex);

    void setNegativeVertex(Vertex *aVertex);

    // print func
    void print();
};

#endif
