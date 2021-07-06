#ifndef _SEGMENTDB_H_
#define _SEGMENTDB_H_

#include "SVprofile.hpp"
#include "LocusDB.hpp"

using namespace std;

struct locus;

class SVprofile;

struct seg {
    int id;
    string chr;
    int start;     // TODO the start and end are (pos - 1), need to modify all related codes
    int end;

    vector<locus *> loci;
};

class SegmentDB {
protected:
    string mChr;
    int mRegionStart;
    int mRegionEnd;

    vector<int> *mBps;
    vector<seg *> *mSegs;

public:
    SegmentDB(string aChr, int aRegionStart, int aRegionEnd);

    ~SegmentDB();

    vector<seg *> *getSegs();

    vector<int> *getBps();

    void readSegs(const char *aFilename);

    void updateBps(SVprofile *aSVprofile);

    void constructSegsFromBps();

    void writeSegs(const char *outFn);

    void print();
};

#endif
