#include <iostream>

#include "Segment.hpp"
#include "Weight.hpp"

using namespace std;

Segment::Segment(int aId, int chr, string aChrom, int aStart, int aEnd, double aCoverage, double aCredibility, double aCopy) {
    mId = aId;
    chrId = chr;
    mChrom = aChrom;
    mStart = aStart;
    mEnd = aEnd;
    mCredibility = aCredibility;

    mHasLowerBoundLimit = true;
    mHasCheckedOrphan = false;
    mIsOrphan = false;

    mWeight = new Weight(aCoverage);
    mWeight->setCopyNum(aCopy);
    mPositiveVertex = new Vertex(aId, '+', mWeight, aCredibility);
    mNegativeVertex = new Vertex(aId, '-', mWeight, aCredibility);
    mPositiveVertex->setSegment(this);
    mNegativeVertex->setSegment(this);
}
Segment::Segment(int aId, int chr, Segment* seg) {
    mId = aId;
    chrId = chr;
    mChrom = seg->getChrom();
    mStart = seg->getStart();
    mEnd = seg->getEnd();
    mCredibility = seg->getCredibility();

    mHasLowerBoundLimit = true;
    mHasCheckedOrphan = false;
    mIsOrphan = false;

    mWeight = new Weight(seg->getWeight()->getCoverage());
    mWeight->setCopyNum(seg->getWeight()->getCopyNum());
    mPositiveVertex = new Vertex(aId, '+', mWeight, mCredibility);
    mNegativeVertex = new Vertex(aId, '-', mWeight, mCredibility);
    mPositiveVertex->setSegment(this);
    mNegativeVertex->setSegment(this);
}

Segment::~Segment() { ; }

int Segment::getId() { return mId; }
int Segment::getChrId() { return chrId; }

string Segment::getChrom() { return mChrom; }

int Segment::getStart() { return mStart; }

int Segment::getEnd() { return mEnd; }

double Segment::getCredibility() { return mCredibility; }

void Segment::setCredibility(double aCredibility) { mCredibility = aCredibility; }

bool Segment::isOrphan(bool aIsOriginal) {
    // pos orphan <=> neg orphan
    if (aIsOriginal) return mIsOrphan;
    else return mPositiveVertex->isOrphan(false) && mNegativeVertex->isOrphan(false);
}

bool Segment::isDeadEnd() {
    return (!mPositiveVertex->hasEdgeAsSource()) || (!mPositiveVertex->hasEdgeAsTarget());
}

bool Segment::hasLowerBoundLimit() { return mHasLowerBoundLimit; }

bool Segment::hasCopy() { return mPositiveVertex->hasCopy(); }

void Segment::checkOrphan() {
    // if (!mHasCheckedOrphan) {
    mPositiveVertex->checkOrphan();
    mNegativeVertex->checkOrphan();
    mIsOrphan = mPositiveVertex->isOrphan() && mNegativeVertex->isOrphan();
    // mHasLowerBoundLimit = !mIsOrphan;
    // mHasCheckedOrphan = true;
    // }
}

void Segment::checkLowerBound() { mHasLowerBoundLimit = mWeight->getCopyNum() >= 1; }

void Segment::setHasLowerBoundLimit() { mHasLowerBoundLimit = true; }

void Segment::resetHasLowerBoundLimit() { mHasLowerBoundLimit = false; }

void Segment::restoreCopy() { mWeight->restore(); }

void Segment::backupCopy() { mWeight->backup(); }

Weight *Segment::getWeight() { return mWeight; }

Vertex *Segment::getPositiveVertex() { return mPositiveVertex; }

Vertex *Segment::getNegativeVertex() { return mNegativeVertex; }

void Segment::setPositiveVertex(Vertex *aVertex) { mPositiveVertex = aVertex; }

void Segment::setNegativeVertex(Vertex *aVertex) { mNegativeVertex = aVertex; }


void Segment::print() {
    // TODO
}

int Segment::getPartition() const {
    return partition;
}

void Segment::setPartition(int partition) {
    Segment::partition = partition;
}
