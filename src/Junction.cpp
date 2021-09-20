#include <iostream>

#include "Junction.hpp"

using namespace std;

Junction::Junction(Segment *aSource, Segment *aTarget, char aSourceDir, char aTargetDir, double aCoverage,
                   double aCredibility, double aCopy, bool aInferred, bool aIsBounded, bool aIsSourceSinkJunction) {
    mSourceDir = aSourceDir;
    mTargetDir = aTargetDir;

    mCredibility = aCredibility;

    mIsInferred = aInferred;
    mHasLowerBoundLimit = aIsBounded;
    // aSource->setHasLowerBoundLimit();
    // aTarget->setHasLowerBoundLimit();

    mWeight = new Weight(aCoverage);
    mWeight->setCopyNum(aCopy);
    aIsSourceSinkJunction ? mWeight->setInferred() : mWeight->resetInferred();

    mSource = aSource;
    mTarget = aTarget;

    if (aSourceDir == '+' && aTargetDir == '+') {
        mEdgeA = new Edge(mSource->getPositiveVertex(), mTarget->getPositiveVertex(), mWeight, mCredibility);
        mEdgeB = new Edge(mTarget->getNegativeVertex(), mSource->getNegativeVertex(), mWeight, mCredibility);
    } else if (aSourceDir == '-' && aTargetDir == '-') {
        mEdgeA = new Edge(mSource->getNegativeVertex(), mTarget->getNegativeVertex(), mWeight, mCredibility);
        mEdgeB = new Edge(mTarget->getPositiveVertex(), mSource->getPositiveVertex(), mWeight, mCredibility);
    } else if (aSourceDir == '+' && aTargetDir == '-') {
        mEdgeA = new Edge(mSource->getPositiveVertex(), mTarget->getNegativeVertex(), mWeight, mCredibility);
        mEdgeB = new Edge(mTarget->getPositiveVertex(), mSource->getNegativeVertex(), mWeight, mCredibility);
    } else if (aSourceDir == '-' && aTargetDir == '+') {
        mEdgeA = new Edge(mSource->getNegativeVertex(), mTarget->getPositiveVertex(), mWeight, mCredibility);
        mEdgeB = new Edge(mTarget->getNegativeVertex(), mSource->getPositiveVertex(), mWeight, mCredibility);
    }

    mEdgeA->setJunction(this);
    mEdgeB->setJunction(this);
}

Junction::~Junction() {
    delete mEdgeA;
    delete mEdgeB;
    delete mWeight;
}

vector<string> Junction::getInfo() {
    vector<string> info;
    info.push_back(mEdgeA->getInfo());
    info.push_back(mEdgeB->getInfo());
    return info;
}

double Junction::getCredibility() { return mCredibility; }

double Junction::setCredibility(double aCredibility) { mCredibility = aCredibility; }

bool Junction::isInferred() { return mIsInferred; }

bool Junction::hasLowerBoundLimit() { return mHasLowerBoundLimit; }

bool Junction::hasCopy() { return mEdgeA->hasCopy(); }

void Junction::setInferred() { mIsInferred = true; }

void Junction::resetInferred() { mIsInferred = false; }

void Junction::setHasLowerBoundLimit() { mHasLowerBoundLimit = true; }

void Junction::resetHasLowerBoundLimit() { mHasLowerBoundLimit = false; }

void Junction::checkLowerBound() { mHasLowerBoundLimit = mWeight->getCopyNum() >= 0 && !mIsInferred; }

void Junction::restoreCopy() { mWeight->restore(); }

void Junction::backupCopy() { mWeight->backup(); }

Weight *Junction::getWeight() { return mWeight; }

Segment *Junction::getSource() { return mSource; }

Segment *Junction::getTarget() { return mTarget; }

Edge *Junction::getEdgeA() { return mEdgeA; }

Edge *Junction::getEdgeB() { return mEdgeB; }

// functionality
void Junction::insertEdgesToVertices() {
    if (mSourceDir == '+' && mTargetDir == '+') {
        mSource->getPositiveVertex()->insertEdgeAsSource(mEdgeA);
        mTarget->getPositiveVertex()->insertEdgeAsTarget(mEdgeA);
        mSource->getNegativeVertex()->insertEdgeAsTarget(mEdgeB);
        mTarget->getNegativeVertex()->insertEdgeAsSource(mEdgeB);
    } else if (mSourceDir == '-' && mTargetDir == '-') {
        mSource->getNegativeVertex()->insertEdgeAsSource(mEdgeA);
        mTarget->getNegativeVertex()->insertEdgeAsTarget(mEdgeA);
        mSource->getPositiveVertex()->insertEdgeAsTarget(mEdgeB);
        mTarget->getPositiveVertex()->insertEdgeAsSource(mEdgeB);
    } else if (mSourceDir == '+' && mTargetDir == '-') {
        mSource->getPositiveVertex()->insertEdgeAsSource(mEdgeA);
        mTarget->getNegativeVertex()->insertEdgeAsTarget(mEdgeA);
        if (mSource != mTarget) {
            mSource->getNegativeVertex()->insertEdgeAsTarget(mEdgeB);
            mTarget->getPositiveVertex()->insertEdgeAsSource(mEdgeB);
        }
    } else if (mSourceDir == '-' && mTargetDir == '+') {
        mSource->getNegativeVertex()->insertEdgeAsSource(mEdgeA);
        mTarget->getPositiveVertex()->insertEdgeAsTarget(mEdgeA);
        if (mSource != mTarget) {
            mSource->getPositiveVertex()->insertEdgeAsTarget(mEdgeB);
            mTarget->getNegativeVertex()->insertEdgeAsSource(mEdgeB);
        }
    }
}

void Junction::print() {
    // TODO
}
