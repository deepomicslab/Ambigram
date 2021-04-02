#ifndef _LOCUS_DB_H_
#define _LOCUS_DB_H_

#include "SegmentDB.hpp"

using namespace std;

struct seg;

struct locus {
    int id;
    int pos;   // 0-based
    char ref;
    char alt;

    seg * belongSeg;
};

class SegmentDB;

class LocusDB {
    protected:
        string mChrom;
        int mStart;
        int mEnd;

        vector<locus *> * mLoci;
        vector<int> * mPos;

        SegmentDB * mSegRef;

        void readVCF(const char * aVCF);
        void readLegend(const char * aLegend);

    public:
        LocusDB(string aChrom, int aStart, int aEnd);
        ~LocusDB();

        string getChr();
        int getStart();
        int getEnd();
        vector<locus *> * getLoci();
        
        void read(const char * aFileName, const int mode);
        void findLociInRange(int startPos, int endPos, vector<locus *>::const_iterator & begin, vector<locus *>::const_iterator & end);
        
        void setSegRef(SegmentDB * aSegDB);
        void assignLocusToSeg();

        void print();

        static const int MODE_VCF = 0, MODE_LEGEND = 1;
};

#endif
