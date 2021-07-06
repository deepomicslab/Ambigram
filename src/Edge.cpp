#include <iostream>

#include "Edge.hpp"

using namespace std;

Edge::Edge(Vertex *aSource, Vertex *aTarget, Weight *aWeight, double aCredibility) {
    mCredibility = aCredibility;

    mIsVisited = false;

    mWeight = aWeight;
    mSource = aSource;
    mTarget = aTarget;
}

Edge::~Edge() {
    vector<Edge *>::iterator removedIter;
    for (vector<Edge *>::iterator i = mSource->getEdgesAsSource()->begin();
         i != mSource->getEdgesAsSource()->end(); i++) {
        if (*i == this) {
            removedIter = i;
            break;
        }
    }
    mSource->getEdgesAsSource()->erase(removedIter);

    for (vector<Edge *>::iterator i = mTarget->getEdgesAsTarget()->begin();
         i != mSource->getEdgesAsTarget()->end(); i++) {
        if (*i == this) {
            removedIter = i;
            break;
        }
    }
    mTarget->getEdgesAsTarget()->erase(removedIter);
}

string Edge::getInfo() {
    return mSource->getInfo() + "=>" + mTarget->getInfo();
}

double Edge::getCredibility() { return mCredibility; }

void Edge::setCredibility(double aCredibility) { mCredibility = aCredibility; }

bool Edge::isVisited() { return mIsVisited; }

void Edge::setVisited() { mIsVisited = true; }

void Edge::resetVisited() { mIsVisited = false; }

Weight *Edge::getWeight() { return mWeight; }

Vertex *Edge::getSource() { return mSource; }

Vertex *Edge::getTarget() { return mTarget; }

Junction *Edge::getJunction() { return mJunction; }

void Edge::setSource(Vertex *aSource) { mSource = aSource; }

void Edge::setTarget(Vertex *aTarget) { mSource = aTarget; }

void Edge::setJunction(Junction *aJunction) { mJunction = aJunction; }

bool Edge::hasCopy() { return mWeight->getCopyNum() >= 1; }

bool Edge::doesConnectSameDir() { return mSource->getDir() == mTarget->getDir(); }

void Edge::recover(double copy) {
    mWeight->increaseCopyNum(copy);
    // mSource->getWeight()->increaseCopyNum(copy);
}

void Edge::traverse() {
    mWeight->decreaseCopyNum();
    // mSource->getWeight()->decreaseCopyNum();
}

Edge *Edge::getComplementEdge() {
    return (this == mJunction->getEdgeA()) ? mJunction->getEdgeB() : mJunction->getEdgeA();
}

void Edge::print() {
    // TODO
    cout << mSource->getId() << mSource->getDir() << "=>" << mTarget->getId() << mTarget->getDir() << endl;
}
