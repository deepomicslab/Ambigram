#include <iostream>
#include <algorithm>

#include "Record.hpp"
#include "Exceptions.hpp"

using namespace std;

Record::Record(string aChrom, int aPos, char aStrand) {
    mChrom = aChrom;
    mPos = aPos;
    mStrand = aStrand;

    mBackwardEntries = new vector<entry_t *>();
    mBackwardEntryUUID = new vector<string>();

    mForwardEntries = new vector<entry_t *>();
    mForwardEntryUUID = new vector<string>();

    mComplementRecord = NULL;
}

Record::~Record() {
    ;
}

string Record::getChrom() { return mChrom; }
int Record::getPos() { return mPos; }
char Record::getStrand() { return mStrand; }

vector<string> * Record::getBackwardEntryUUID() { return mBackwardEntryUUID; }
vector<string> * Record::getForwardEntryUUID() { return mForwardEntryUUID; }
vector<entry_t *> * Record::getBackwardEntries() { return mBackwardEntries; }
vector<entry_t *> * Record::getForwardEntries() { return mForwardEntries; }

void Record::insertBackwardEntry(string aChrom, int aPos, char aStrand, int aSupport, bool aIsComplement) {
    string uuid = (aStrand == '+') ? (aChrom + ':' + to_string(aPos)) : (aChrom + ':' + to_string(-aPos));
    vector<string>::iterator iterUUID = lower_bound(mBackwardEntryUUID->begin(), mBackwardEntryUUID->end(), uuid);
    vector<entry_t *>::iterator iterEntry = iterUUID - mBackwardEntryUUID->begin() + mBackwardEntries->begin();

    if (iterUUID == mBackwardEntryUUID->end() || *iterUUID != uuid) {
        // not exist
        mBackwardEntryUUID->insert(iterUUID, uuid);

        entry_t * e = new entry_t{aChrom, aPos, aStrand, aSupport, aIsComplement};
        mBackwardEntries->insert(iterEntry, e);
    } else {
        // exist
        (*iterEntry)->support += aSupport;
    }
}

void Record::insertForwardEntry(string aChrom, int aPos, char aStrand, int aSupport, bool aIsComplement) {
    string uuid = (aStrand == '+') ? (aChrom + ':' + to_string(aPos)) : (aChrom + ':' + to_string(-aPos));
    vector<string>::iterator iterUUID = lower_bound(mForwardEntryUUID->begin(), mForwardEntryUUID->end(), uuid);
    vector<entry_t *>::iterator iterEntry = iterUUID - mForwardEntryUUID->begin() + mForwardEntries->begin();

    if (iterUUID == mForwardEntryUUID->end() || *iterUUID != uuid) {
        // not exist
        mForwardEntryUUID->insert(iterUUID, uuid);

        entry_t * e = new entry_t{aChrom, aPos, aStrand, aSupport, aIsComplement};
        mForwardEntries->insert(iterEntry, e);
    } else {
        // exist
        (*iterEntry)->support += aSupport;
    }
}

void Record::sortEntry() {
    sort(mBackwardEntries->begin(), mBackwardEntries->end(), [](entry_t * e1, entry_t * e2) { return e1->support < e2->support; });
    sort(mForwardEntries->begin(), mForwardEntries->end(), [](entry_t * e1, entry_t * e2) { return e1->support < e2->support; });
}

entry_t * Record::findForwardEntry(string aChrom, int aPos, char aStrand) {
    int id = (aStrand == '+') ? aPos : -aPos;
    string uuid = aChrom + ':' + to_string(id);
    vector<string>::iterator iterUUID = lower_bound(mForwardEntryUUID->begin(), mForwardEntryUUID->end(), uuid);
    vector<entry_t *>::iterator iterEntry = iterUUID - mForwardEntryUUID->begin() + mForwardEntries->begin();

    if (iterUUID == mForwardEntryUUID->end() || *iterUUID != uuid) {
        return NULL;
    } else {
        return *iterEntry;
    }
}

void Record::print() {
    // TODO
}
