#include <iostream>
#include <algorithm>

#include "LocusDB.hpp"
#include "htslib/synced_bcf_reader.h"

using namespace std;

LocusDB::LocusDB(string aChrom, int aStart, int aEnd) {
    mChrom = aChrom;
    mStart = aStart;
    mEnd = aEnd;

    mLoci = new vector<locus *>();
    mPos = new vector<int>();
}

LocusDB::~LocusDB() {
    delete mLoci;
    delete mPos;
}

string LocusDB::getChr() { return mChrom; }
int LocusDB::getStart() { return mStart; }
int LocusDB::getEnd() { return mEnd; }
vector<locus *> * LocusDB::getLoci() { return mLoci; }

void LocusDB::read(const char * aFileName, const int mode) {
    if (mode == this->MODE_VCF) this->readVCF(aFileName);
    else if (mode == this->MODE_LEGEND) this->readLegend(aFileName);
}

void LocusDB::readVCF(const char * aVCF) {
    bcf_srs_t * sr = bcf_sr_init();
    bcf_sr_set_regions(sr, (mChrom + ':' + to_string(mStart) + '-' + to_string(mEnd)).c_str(), 0);
    bcf_sr_add_reader(sr, aVCF);
    // bcf_hdr_t * vcf_header = bcf_sr_get_header(sr, 0);
    bcf1_t * rec = bcf_init1();

    int id = 0;
    while (bcf_sr_next_line(sr)) {
        rec = bcf_sr_get_line(sr, 0);
        locus * l = new locus{id, rec->pos, rec->d.allele[0][0], rec->d.allele[1][0]};
        mLoci->push_back(l);
        mPos->push_back(rec->pos);
        id++;
        // cout << rec->pos << " " << rec->d.allele[0] << " " << rec->d.allele[1] << endl;
    }
    bcf_sr_destroy(sr);
}

void LocusDB::readLegend(const char * aLegend) {}

void LocusDB::findLociInRange(int startPos, int endPos, vector<locus *>::const_iterator & begin, vector<locus *>::const_iterator & end) {
    begin = lower_bound(mPos->begin(), mPos->end(), startPos) - mPos->begin() + mLoci->begin();
    end = upper_bound(mPos->begin(), mPos->end(), endPos) - mPos->begin() + mLoci->begin();
}

void LocusDB::setSegRef(SegmentDB * aSegDB) { mSegRef = aSegDB; }

void LocusDB::assignLocusToSeg() {
    int segIdx = 0;
    vector<seg *> * segs = mSegRef->getSegs();
    for (locus * l : *mLoci) {
        while (l->pos > (*segs)[segIdx]->end) {
            segIdx++;
        }
        l->belongSeg = (*segs)[segIdx];
        (*segs)[segIdx]->loci.push_back(l);
    }
}

void LocusDB::print() {
    for (locus * l : *mLoci) {
        cout << l->pos << " " << l->ref << " " << l->alt;
        if (l->belongSeg != NULL) cout << " " << l->belongSeg->id;
        cout << endl;
    }
}
