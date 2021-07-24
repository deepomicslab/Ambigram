#ifndef _COVERAGE_H_
#define _COVERAGE_H_

using namespace std;

class Weight {
protected:
    double mCoverage;           // coverage as adjusted by LP
    // double mCoverageOriginal;   // coverage from input
    // double mCoverageAdjusted;   // coverage as adjusted by purity
    double mCopyNum;             // the copy number, coverage / average coverage of one copy
    double mCopyNumOriginal;
    double mCopyNumBackup;
    double mCorrectedCoverage;

    bool mIsInferred;

public:
    Weight(double aCoverage);

    ~Weight();

    double getCoverage();
    double getCorrectedCoverage();
    void setCorrectedCoverage(double v);

    // double getOriginalCoverage();
    // double getAdjustedCoverage();
    double getCopyNum();

    double getCopyNumBackup();

    void setCoverage(double aCoverage);

    // void setOriginalCoverage(double aCoverage);
    // void setAdjustedCoverage(double aCoverage);
    void setCopyNum(double aCopyNum);

    void backup();

    void restore();

    void increaseCopyNum(double aIncrement = 1);

    void decreaseCopyNum(double aDecrement = 1);

    bool isInferred();

    void setInferred();

    void resetInferred();

    void print();
};

#endif
