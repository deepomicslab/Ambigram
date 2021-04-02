#ifndef _SV_PROFILE_H_
#define _SV_PROFILE_H_

#include <tuple>
#include <vector>

#include "SegmentDB.hpp"

using namespace std;

typedef tuple<string, int, char, string, int, char> SVInfo;     // left_chr, left_pos<left_id>, left_strand, right_chr, right_pos<right_id>, right_strand, index
typedef tuple<int, string, int, string> SVSupport;   // left_clip, right_clip, left_cigar, right_cigar
// typedef tuple<string, int, int, int, double> seg;  // chr, idx, start, end, depth

class SegmentDB;

class SVprofile {
    protected:
        string mSample;
        vector<int> * mBps;
        vector<int> * _mBps;

        vector<SVInfo> * mAbnormalJuncInfo;
        vector<SVSupport> * mAbnormalJuncSupports;
        vector<SVInfo> * _mAbnormalJuncInfo;
        vector<SVSupport> * _mAbnormalJuncSupports;

        vector<SVInfo> * mNormalJuncInfo;
        vector<SVSupport> * mNormalJuncSupports;

        double mAvgDp;
        double * mSegDepth;
        SegmentDB * mSegs;

    public:
        SVprofile(const char * aRawSVFn, string aSample);
        ~SVprofile();

        vector<int> * getBps();
        vector<SVInfo> * getAbnormalInfo();
        vector<SVSupport> * getAbnormalSupports();

        void readRawSV(const char * aRawSVFn);
        void insertSVEntry(string leftChr, int leftPos, char leftStrand, int leftClip, string leftCigar, string rightChr, int rightPos, char rightStrand, int rightClip, string rightCigar);
        int getMatchNum(string);
        void filterAbnormal(string chrom, int start, int end, int aClipThres = 5, int aMatchThres = 19);
        void setSegDB(SegmentDB * aSegDB);
        void pos2id();
        void countSegDepth(const char * aDepthFn);
        void countNormal(const char * aBamFn, int aEndMatchThres = 5);

        void writeLocalHap(const char * outFn);
        void writeNormal(const char * outFn);
        void writeAbnormal(const char * outFn);

        void print();
        void printJunc();
};

#endif
