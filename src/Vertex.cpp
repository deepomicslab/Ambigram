#include <iostream>

#include "Vertex.hpp"

using namespace std;

Vertex::Vertex(int aId, char aDir, Weight *aWeight, double aCredibility) {
    mId = aId;
    mDir = aDir;

    mWeight = aWeight;
    mCredibility = aCredibility;

    mIsOrphan = false;
    mHasCheckedOrphan = false;
    mIsVisited = false;

    mShortestPrevEdge = NULL;
    mSegment = NULL;

    mEdgesAsSource = new vector<Edge *>();
    mEdgesAsTarget = new vector<Edge *>();
}

int Vertex::getId() { return mId; }

int Vertex::getStart() { return (mDir == '+') ? mSegment->getStart() : mSegment->getEnd(); }

int Vertex::getEnd() { return (mDir == '+') ? mSegment->getEnd() : mSegment->getStart(); }

char Vertex::getDir() { return mDir; }

string Vertex::getInfo() { return to_string(mId) + mDir; }

double Vertex::getInCoverage() {
    double inCoverage = 0;
    for (Edge *e : *mEdgesAsTarget) {
        inCoverage += e->getWeight()->getCoverage();
    }
    return inCoverage;
}

double Vertex::getOutCoverage() {
    double outCoverage = 0;
    for (Edge *e : *mEdgesAsSource) {
        outCoverage += e->getWeight()->getCoverage();
    }
    return outCoverage;
}

double Vertex::getCredibility() { return mCredibility; }

void Vertex::setCredibility(double aCredibility) { mCredibility = aCredibility; }

bool Vertex::isVisited() { return mIsVisited; }

bool Vertex::hasCopy() { return mWeight->getCopyNum() >= 1; }

bool Vertex::hasEdgeAsSource() { return mEdgesAsSource->size() > 0; }

bool Vertex::hasEdgeAsTarget() { return mEdgesAsTarget->size() > 0; }

bool Vertex::isOrphan(bool aIsOriginal) {
    // if (aIsOriginal) return mIsOrphan;
    // else return (!this->hasEdgeAsSource()) && (!this->hasEdgeAsTarget());
    return (!this->hasEdgeAsSource()) && (!this->hasEdgeAsTarget());
}

void Vertex::checkOrphan() {
    // if (!mHasCheckedOrphan) {
    mIsOrphan = this->isOrphan(false);
    // mHasCheckedOrphan = true;
    // }
}

void Vertex::setVisited() { mIsVisited = true; }

void Vertex::resetVisited() {
    mIsVisited = false;
    // mShortestPrevEdge = NULL;
}

void Vertex::traverse() { mWeight->decreaseCopyNum(); }

void Vertex::recover() { mWeight->increaseCopyNum(); }

Weight *Vertex::getWeight() { return mWeight; }

Edge *Vertex::getShortestPrevEdge() { return mShortestPrevEdge; }

Segment *Vertex::getSegment() { return mSegment; }

void Vertex::setShortestPrevEdge(Edge *aEdge) { mShortestPrevEdge = aEdge; }

void Vertex::setSegment(Segment *aSeg) { mSegment = aSeg; }

vector<Edge *> *Vertex::getEdgesAsSource() { return mEdgesAsSource; }

vector<Edge *> *Vertex::getEdgesAsTarget() { return mEdgesAsTarget; }

void Vertex::insertEdgeAsSource(Edge *aEdge) { mEdgesAsSource->push_back(aEdge); }

void Vertex::insertEdgeAsTarget(Edge *aEdge) { mEdgesAsTarget->push_back(aEdge); }

Edge *Vertex::findEdgeAsSource(Vertex *aTargetVertex) {
    for (Edge *e : *mEdgesAsSource) {
        if (e->getTarget() == aTargetVertex) {
            return e;
        }
    }
    return NULL;
}

Vertex *Vertex::getComplementVertex() {
    return (mDir == '-') ? mSegment->getPositiveVertex() : mSegment->getNegativeVertex();
}

void print() {
    // TODO
}
