#include "Weight.hpp"

Weight::Weight(double aCoverage) {
    mCoverage = aCoverage;
    mCorrectedCoverage = mCoverage;
    // mCoverageOriginal = aCoverage;
    // mCoverageAdjusted = aCoverage;
    mCopyNum = 0;
    mCopyNumOriginal = 0;
    mCopyNumBackup = 0;

    mIsInferred = false;
}

Weight::~Weight() { ; }

double Weight::getCoverage() { return mCoverage; }

double Weight::getCorrectedCoverage() { return mCorrectedCoverage; }
void Weight::setCorrectedCoverage(double v) { mCorrectedCoverage = v; }


// double Weight::getOriginalCoverage() { return mCoverageOriginal; }
// double Weight::getAdjustedCoverage() { return mCoverageAdjusted; }
double Weight::getCopyNum() { return mCopyNum; }

double Weight::getCopyNumBackup() { return mCopyNumBackup; }

void Weight::setCoverage(double aCoverage) { mCoverage = aCoverage; }

// void Weight::setOriginalCoverage(double aCoverage) { mCoverageOriginal = aCoverage; }
// void Weight::setAdjustedCoverage(double aCoverage) { mCoverageAdjusted = aCoverage; }
void Weight::setCopyNum(double aCopyNum) { mCopyNum = mCopyNumOriginal = mCopyNumBackup = aCopyNum; }

void Weight::backup() { mCopyNumBackup = mCopyNum; }

void Weight::restore() { mCopyNum = mCopyNumBackup; }

void Weight::increaseCopyNum(double aIncrement) { mCopyNum += aIncrement; }

void Weight::decreaseCopyNum(double aDecrement) { mCopyNum -= aDecrement; }

bool Weight::isInferred() { return mIsInferred; }

void Weight::setInferred() { mIsInferred = true; }

void Weight::resetInferred() { mIsInferred = false; }
