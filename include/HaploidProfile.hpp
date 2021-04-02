#ifndef _HAPLOID_PROFILE_H_
#define _HAPLOID_PROFILE_H_

#include "SegmentDB.hpp"
#include "LocusDB.hpp"
#include "SupportProfile.hpp"

using namespace std;

struct SignedSeg;

struct Variant {
    locus * l;
    int type;
};

struct Strand {
    // SignedSeg * belongSeg;
    seg * belongSeg;
    char sign;
    vector<Variant *> variants;

    // Strand * compStrand;
};

// struct SignedSeg {
//     seg * oriSeg;
//     Strand * sp, * sn;
//     SignedSeg(seg * aSeg);
// };

class HaploidProfile {
    protected:
        string mSampleName;

        SegmentDB * mSegRef;
        bool * mSegNormal;

        vector<Strand *> * mHaploid1;
        vector<Strand *> * mHaploid2;

        SupportProfile * mSp;

    public:
        HaploidProfile(string aSampleName);
        ~HaploidProfile();

        void readHaploids(const char * aFilename);
        void setSegRef(SegmentDB * aSegDB);
        void identifyNormal();
    
        void setSupportProfile(SupportProfile * sp);
        void placeVariantsInSeg(seg * aSeg);
        void placeVariants();

        void print();
        void printNormal();
};

#endif
