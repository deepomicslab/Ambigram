#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>

#include "htslib/sam.h"
#include "SVprofile.hpp"

using namespace std;
SVprofile::SVprofile(const char * aRawSVFn, string aSample) {
    mSample = aSample;
    mBps = new vector<int>();
    mAbnormalJuncInfo = new vector<SVInfo>();
    mAbnormalJuncSupports = new vector<SVSupport>();
    _mBps = new vector<int>();
    _mAbnormalJuncInfo = new vector<SVInfo>();
    _mAbnormalJuncSupports = new vector<SVSupport>();
    mNormalJuncInfo = new vector<SVInfo>();
    mNormalJuncSupports = new vector<SVSupport>();
    mAvgDp = 0;
    mSegDepth = NULL;

    this->readRawSV(aRawSVFn);
}

SVprofile::~SVprofile() {
    delete [] mBps;
    delete [] mAbnormalJuncInfo;
    delete [] mAbnormalJuncSupports;
    delete [] _mBps;
    delete [] _mAbnormalJuncInfo;
    delete [] _mAbnormalJuncSupports;
}

vector<int> * SVprofile::getBps() { return mBps; }
vector<SVInfo> * SVprofile::getAbnormalInfo() { return mAbnormalJuncInfo; }
vector<SVSupport> * SVprofile::getAbnormalSupports() { return mAbnormalJuncSupports; }

void SVprofile::readRawSV(const char * aRawSVFn) {
    ifstream rawSV(aRawSVFn);
    if (!rawSV) {
        cerr << "Cannot open file: " << aRawSVFn << endl;
        exit(1);
    }

    string line;
    getline(rawSV, line); // title line
    string temp;
    while(getline(rawSV, line)) {
        stringstream ss(line);
        int leftPos, leftClip, rightPos, rightClip;
        char leftStrand, rightStrand;
        string leftChr, rightChr, leftCigar, rightCigar;
        ss >> leftChr >> leftPos >> leftStrand >> leftClip
            >> rightChr >> rightPos >> rightStrand >> rightClip;
        for (int i = 0; i < 11; i++) ss >> temp;
        ss >> leftCigar >> rightCigar;
        this->insertSVEntry(leftChr, leftPos - 1, leftStrand, leftClip, leftCigar, rightChr, rightPos - 1, rightStrand, rightClip, rightCigar);
    }
    rawSV.close();
    
    *_mBps = *mBps;
    *_mAbnormalJuncInfo = *mAbnormalJuncInfo;
    *_mAbnormalJuncSupports = *mAbnormalJuncSupports;
}

void SVprofile::insertSVEntry(string leftChr, int leftPos, char leftStrand, int leftClip, string leftCigar, string rightChr, int rightPos, char rightStrand, int rightClip, string rightCigar) {
    SVInfo info = make_tuple(leftChr, leftPos, leftStrand, rightChr, rightPos, rightStrand);
    SVSupport support = make_tuple(leftClip, leftCigar, rightClip, rightCigar);
    vector<SVInfo>::iterator iterSVInfo = lower_bound(mAbnormalJuncInfo->begin(), mAbnormalJuncInfo->end(), info);
    vector<SVSupport>::iterator iterSVSupport = iterSVInfo - mAbnormalJuncInfo->begin() + mAbnormalJuncSupports->begin();
    if (iterSVInfo == mAbnormalJuncInfo->end() || *iterSVInfo != info) {
        mAbnormalJuncInfo->insert(iterSVInfo, info);
        mAbnormalJuncSupports->insert(iterSVSupport, support);
        
        vector<int>::iterator iterBps = lower_bound(mBps->begin(), mBps->end(), leftPos);
        if (iterBps == mBps->end() || *iterBps != leftPos) {
            mBps->insert(iterBps, leftPos);
        }
        iterBps = lower_bound(mBps->begin(), mBps->end(), rightPos);
        if (iterBps == mBps->end() || *iterBps != rightPos) {
            mBps->insert(iterBps, rightPos);
        }
    } else {
        get<0>(*iterSVSupport) += leftClip;
        get<2>(*iterSVSupport) += rightClip;
    }
}

int SVprofile::getMatchNum(string cigar) {
    int match = 0;
    size_t start = 0;
    size_t found = cigar.find_first_of("MIDNSHP=X", start);
    while (found != string::npos) {
        int num = stoi(cigar.substr(start, found - start));
        char type = cigar[found];
        if (type == 'M') {
            match += num;
        }
        start = found + 1;
        found = cigar.find_first_of("MIDNSHP=X", start);
    }
    return match;
}

void SVprofile::filterAbnormal(string chrom, int start, int end, int aClipThres, int aMatchThres) {
    mBps->clear();
    mAbnormalJuncInfo->clear();
    mAbnormalJuncSupports->clear();
    for (int i = 0; i < _mAbnormalJuncInfo->size(); i++) {
        SVSupport support = (*_mAbnormalJuncSupports)[i];
        int leftPos, rightPos, leftClip, rightClip;
        char leftStrand, rightStrand;
        string leftChr, leftCigar, rightChr, rightCigar;
        tie(leftChr, leftPos, leftStrand, rightChr, rightPos, rightStrand) = (*_mAbnormalJuncInfo)[i];
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*_mAbnormalJuncSupports)[i];
        if (leftChr == chrom && rightChr == chrom
                && leftPos >= start && leftPos <= end && rightPos >= start && rightPos <= end
                && this->getMatchNum(leftCigar) >= aMatchThres && this->getMatchNum(rightCigar) >=aMatchThres
                && leftClip + rightClip >= aClipThres) {
            this->insertSVEntry(leftChr, leftPos, leftStrand, leftClip, leftCigar, rightChr, rightPos, rightStrand, rightClip, rightCigar);
        }
    }
}

void SVprofile::setSegDB(SegmentDB * aSegDB) {
    mSegs = aSegDB;
}

void SVprofile::pos2id() {
    for (int i = 0; i < mAbnormalJuncInfo->size(); i++) {
        int idLeft = lower_bound(mSegs->getBps()->begin(), mSegs->getBps()->end(), get<1>((*mAbnormalJuncInfo)[i])) - mSegs->getBps()->begin();
        int idRight = lower_bound(mSegs->getBps()->begin(), mSegs->getBps()->end(), get<4>((*mAbnormalJuncInfo)[i])) - mSegs->getBps()->begin() + 1;
        get<1>((*mAbnormalJuncInfo)[i]) = idLeft;
        get<4>((*mAbnormalJuncInfo)[i]) = idRight;
    }
}

void SVprofile::countSegDepth(const char * aDepthFn) {
    if (mSegDepth != NULL) {
        delete [] mSegDepth;
    }
    mSegDepth = new double[mSegs->getSegs()->size()];

    gzFile depthFin = gzopen(aDepthFn, "rb");
    if (!depthFin) {
        cerr << "cannot open file: " << aDepthFn << endl;
        exit(1);
    }
    char line[8196];
    int segsIdx = 0;
    string segChr;
    int segId, segStart, segEnd;
    double segDepth;
    seg * s = (*(mSegs->getSegs()))[segsIdx];
    int totDepth = 0;
    while (gzgets(depthFin, line, 8196) != NULL) {
        stringstream ss(line);
        string chr;
        int pos;
        double depth;
        ss >> chr >> pos >> depth;
        if (chr != s->chr) continue;
        totDepth += depth;
        mAvgDp += depth;
        if (pos - 1 == s->end) {
            mSegDepth[segsIdx] = totDepth * 1.0 / (s->end - s->start + 1);
            totDepth = depth;
            segsIdx++;
            if (segsIdx >= mSegs->getSegs()->size()) break;
            s = (*(mSegs->getSegs()))[segsIdx];
        }
    }
    gzclose(depthFin);
    mAvgDp = mAvgDp / (mSegs->getSegs()->back()->end - mSegs->getSegs()->front()->start + 1);
}

void SVprofile::countNormal(const char * aBamFn, int aEndMatchThres) {
    mNormalJuncInfo->clear();
    mNormalJuncSupports->clear();
    samFile * bam = sam_open(aBamFn, "rb");
    hts_idx_t * bamIdx = sam_index_load(bam, aBamFn);
    bam_hdr_t * bamHeader = sam_hdr_read(bam);
    hts_itr_t * itr = NULL;
    bam1_t * aln = bam_init1();
    for (int i = 0; i < mSegs->getSegs()->size() - 1; i++) {
        seg * s = (*(mSegs->getSegs()))[i];
        itr = bam_itr_queryi(bamIdx, bam_name2id(bamHeader, s->chr.c_str()), s->end, s->end + 1);
        int support = 0;
        while (bam_itr_next(bam, itr, aln) >= 0) {
            uint32_t * cigar = bam_get_cigar(aln);
            int coveredLen = bam_cigar2rlen(aln->core.n_cigar, cigar);
            if (s->end - aln->core.pos + 1 >= aEndMatchThres & aln->core.pos + 1 + coveredLen - s->end >= aEndMatchThres) {
                support++;
            }
        }
        mNormalJuncInfo->push_back(SVInfo(s->chr, s->id, '+', s->chr, s->id + 1, '+'));
        mNormalJuncSupports->push_back(SVSupport(support, "", support, ""));
    }
}

void SVprofile::writeLocalHap(const char * outFn) {
    ofstream fout(outFn);
    fout << "SAMPLE " << mSample << endl
        << "AVG_DP " << mAvgDp << endl
        << "PURITY 1" << endl
        << "AVG_PLOIDY 2" << endl
        << "PLOIDY 2m1" << endl
        << "SOURCE H:" << mSegs->getSegs()->front()->id + 1 << endl
        << "SINK H:" << mSegs->getSegs()->back()->id + 1 << endl;
    for (seg * s : *(mSegs->getSegs())) {
        fout << "SEG H:" << s->id + 1 << " " << mSegDepth[s->id] << endl; 
    }
    for (int i = 0; i < mAbnormalJuncInfo->size(); i++) {
        string leftChr, rightChr;
        int leftId, rightId;
        char leftStrand, rightStrand;
        tie(leftChr, leftId, leftStrand, rightChr, rightId, rightStrand) = (*mAbnormalJuncInfo)[i];
        int leftClip, rightClip;
        string leftCigar, rightCigar;
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mAbnormalJuncSupports)[i];
        fout << "JUNC H:" << leftId + 1 << ":" << leftStrand << " H:" << rightId + 1 << ":" << rightStrand << " " << leftClip + rightClip << endl;
    }
    for (int i = 0; i < mNormalJuncInfo->size(); i++) {
        string leftChr, rightChr;
        int leftId, rightId;
        char leftStrand, rightStrand;
        tie(leftChr, leftId, leftStrand, rightChr, rightId, rightStrand) = (*mNormalJuncInfo)[i];
        int leftClip, rightClip;
        string leftCigar, rightCigar;
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mNormalJuncSupports)[i];
        fout << "JUNC H:" << leftId + 1 << ":" << leftStrand << " H:" << rightId + 1 << ":" << rightStrand << " " << leftClip << endl;
    }
    fout.close();
}

void SVprofile::writeAbnormal(const char * outFn) {
    ofstream fout(outFn);
    for (int i = 0; i < mAbnormalJuncInfo->size(); i++) {
        string leftChr, rightChr;
        int leftId, rightId;
        char leftStrand, rightStrand;
        tie(leftChr, leftId, leftStrand, rightChr, rightId, rightStrand) = (*mAbnormalJuncInfo)[i];
        int leftClip, rightClip;
        string leftCigar, rightCigar;
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mAbnormalJuncSupports)[i];
        fout << leftId + 1 << " " << leftStrand << " " << rightId + 1 << " " << rightStrand << " " << leftClip + rightClip << endl;
    }
    fout.close();
}

void SVprofile::writeNormal(const char * outFn) {
    ofstream fout(outFn);
    for (int i = 0; i < mNormalJuncInfo->size(); i++) {
        string leftChr, rightChr;
        int leftId, rightId;
        char leftStrand, rightStrand;
        tie(leftChr, leftId, leftStrand, rightChr, rightId, rightStrand) = (*mNormalJuncInfo)[i];
        int leftClip, rightClip;
        string leftCigar, rightCigar;
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mNormalJuncSupports)[i];
        fout << leftId + 1 << " " << leftStrand << " " << rightId + 1 << " " << rightStrand << " " << leftClip + rightClip << endl;
    }
    fout.close();
}

void SVprofile::print() {
    cout << "Before filter: " << _mAbnormalJuncInfo->size() << ", after filter: " << mAbnormalJuncInfo->size() << endl;
    for (vector<int>::iterator i = mBps->begin(); i != mBps->end(); i++) {
        cout << *i;
        if (i != mBps->begin() && (i - mBps->begin()) % 10 == 0) cout << endl; else cout << " ";
    }
    cout << endl;
    for (int i = 0; i < mAbnormalJuncInfo->size(); i++) {
        int leftPos, rightPos, leftClip, rightClip;
        char leftStrand, rightStrand;
        string leftChr, leftCigar, rightChr, rightCigar;
        tie(leftChr, leftPos, leftStrand, rightChr, rightPos, rightStrand) = (*mAbnormalJuncInfo)[i];
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mAbnormalJuncSupports)[i];
        cout << leftChr << " " << leftPos << " " << leftStrand << " " << leftClip << " " << leftCigar << " " << rightChr << " " << rightPos << " " << rightStrand << " " << rightClip << " " << rightCigar << endl;
    }
}

void SVprofile::printJunc() {
    for (int i = 0; i < mAbnormalJuncInfo->size(); i++) {
        string leftChr, rightChr;
        int leftId, rightId;
        char leftStrand, rightStrand;
        tie(leftChr, leftId, leftStrand, rightChr, rightId, rightStrand) = (*mAbnormalJuncInfo)[i];
        int leftClip, rightClip;
        string leftCigar, rightCigar;
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mAbnormalJuncSupports)[i];
        cout << leftId + 1 << " " << leftStrand << " " << rightId + 1 << " " << rightStrand << " " << leftClip + rightClip << endl; 
    }
    for (int i = 0; i < mNormalJuncInfo->size(); i++) {
        string leftChr, rightChr;
        int leftId, rightId;
        char leftStrand, rightStrand;
        tie(leftChr, leftId, leftStrand, rightChr, rightId, rightStrand) = (*mNormalJuncInfo)[i];
        int leftClip, rightClip;
        string leftCigar, rightCigar;
        tie(leftClip, leftCigar, rightClip, rightCigar) = (*mNormalJuncSupports)[i];
        cout << leftId + 1 << " " << leftStrand << " " << rightId + 1 << " " << rightStrand << " " << leftClip << endl; 
    }
}
