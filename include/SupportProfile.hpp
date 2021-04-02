#ifndef _SUPPORT_PROFILE_H_
#define _SUPPORT_PROFILE_H_

#include "htslib/sam.h"
#include "LocusDB.hpp"

using namespace std;

struct readCount {
    int rr, ra, ar, aa;
};

struct support {
    vector<locus *> pairedLoci;
    vector<readCount *> pairedCounts;
};

class SupportProfile {
    protected:
        string mSampleName;
        vector< pair<locus *, support *> > * mCoveredLociSupports;

        LocusDB * mLociRef;
        int * mGT;
        int mNumHom;
        int mNumHet;
        int mNumUnknown;

    public:
        SupportProfile(string aSampleName);
        ~SupportProfile();

        string getSampleName();
        vector< pair<locus *, support *> > * getLociSupports();
        int * getGT();
        
        void setLociRef(LocusDB * db);
        void readGenotypes(const char * aVCF);
        void readSupport(const char * aSupportFn);

        void countSupport(const char * aBamFn);
        void writeSupport(const char * outFn);

        void getInSameSegSupports(locus * l, vector<locus *> & pLoci, vector<readCount *> & pCounts);
        
        int getBaseIdx(bam1_t * aln, int pos);

        void print();
        void printStatistics();
};

#endif
