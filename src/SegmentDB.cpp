#include <iostream>
#include <sstream>
#include <fstream>

#include "SegmentDB.hpp"
#include "SVprofile.hpp"

using namespace std;

SegmentDB::SegmentDB(string aChr, int aRegionStart, int aRegionEnd) {
    mChr = aChr;
    mRegionStart = aRegionStart;
    mRegionEnd = aRegionEnd;

    mBps = new vector<int>();
    mSegs = new vector<seg *>();
}

SegmentDB::~SegmentDB() {
    delete[] mBps;
    delete[] mSegs;
}

vector<seg *> *SegmentDB::getSegs() { return mSegs; }

vector<int> *SegmentDB::getBps() { return mBps; }

void SegmentDB::readSegs(const char *aFilename) {
    ifstream fin(aFilename);
    string line;
    string chr;
    int id, start, end;
    while (getline(fin, line)) {
        stringstream ss(line);
        ss >> chr >> id >> start >> end;
        seg *s = new seg{id - 1, chr, start - 1, end - 1};
        mSegs->push_back(s);
        mBps->push_back(end - 1);
    }
    mBps->pop_back();
    fin.close();
}

void SegmentDB::updateBps(SVprofile *aSVprofile) {
    vector<int> *bps = aSVprofile->getBps();
    for (int i = 0; i < bps->size(); i++) {
        vector<int>::iterator iterBps = lower_bound(mBps->begin(), mBps->end(), (*bps)[i]);
        if (iterBps == mBps->end() || *iterBps != (*bps)[i]) {
            mBps->insert(iterBps, (*bps)[i]);
        }
    }
}

void SegmentDB::constructSegsFromBps() {
    mSegs->clear();
    mSegs->push_back(new seg{0, mChr, mRegionStart, mBps->front()});
    for (int i = 0; i < mBps->size() - 1; i++) {
        mSegs->push_back(new seg{i + 1, mChr, (*mBps)[i], (*mBps)[i + 1]});
    }
    mSegs->push_back(new seg{(int) mBps->size(), mChr, mBps->back(), mRegionEnd});
}

void SegmentDB::writeSegs(const char *outFn) {
    ofstream fout(outFn);
    for (seg *s : *mSegs) {
        fout << s->chr << " " << s->id + 1 << " " << s->start + 1 << " " << s->end + 1 << endl;
    }
    fout.close();
}

void SegmentDB::print() {
    for (seg *s : *mSegs) {
        cout << s->chr << " " << s->id << " " << s->start << " " << s->end << endl;
    }
}
