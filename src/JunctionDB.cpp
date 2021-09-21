#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <functional>
#include <algorithm>


#include "JunctionDB.hpp"
#include "Junction.hpp"

using namespace std;

// JunctionDB::JunctionDB(string chr, int start, int end) {
//     mChr = chr;
//     mStart = start;
//     mEnd = end;

//     // mEntries = new vector<junc_entry_t *>();
//     // mEntryId = new vector<int>();
//     mRecords = new vector<Record *>();
//     mRecordUUID = new vector<int>();
// }

JunctionDB::JunctionDB(const char *aFilename) {
    // mEntries = new vector<junc_entry_t *>();
    // mEntryId = new vector<int>();
    mRecords = new vector<Record *>();
    mRecordUUID = new vector<string>();

    this->readDB(aFilename);
}

JunctionDB::JunctionDB(vector<Junction *> &junctions) {
    mRecords = new vector<Record *>();
    mRecordUUID = new vector<string>();
    for (Junction *junc: junctions) {
        if (junc->getWeight() > 0) {
            string chrom5p = junc->getSource()->getChrom(),
                chrom3p = junc->getTarget()->getChrom();
            char strand5p = junc->getSourceDir(), 
                strand3p = junc->getTargetDir();
            int pos5p = 0, pos3p = 0;
            if (strand5p == '+' && strand3p == '+') {
                pos5p = junc->getSource()->getEnd();
                pos3p = junc->getTarget()->getStart();
            }
            else if (strand5p == '+' && strand3p == '-') {
                pos5p = junc->getSource()->getEnd();
                pos3p = junc->getTarget()->getEnd();
            }
            else if (strand5p == '-' && strand3p == '+') {
                pos5p = junc->getSource()->getStart();
                pos3p = junc->getTarget()->getStart();
            }
            else if (strand5p == '-' && strand3p == '-') {
                pos5p = junc->getSource()->getStart();
                pos3p = junc->getTarget()->getEnd();
            }
            char support = junc->getWeight()->getCoverage();
            this->insertRecord(chrom5p, pos5p, strand5p,
                            chrom3p, pos3p, strand3p,
                            support);
        }         
    }
}

JunctionDB::~JunctionDB() {
    // delete [] mEntries;
    // delete [] mEntryId;
    delete[] mRecords;
    delete[] mRecordUUID;
}

vector<Record *> *JunctionDB::getRecords() { return mRecords; }

void JunctionDB::readDB(const char *aFilename) {
    ifstream dbFile(aFilename);
    if (!dbFile) {
        cerr << "Cannot open file " << aFilename << endl;
        exit(1);
    }

    cout << "Reading database..." << endl;

    char line[8192];
    char *token;
    dbFile.getline(line, 8192);
    while (!dbFile.eof()) {
        dbFile.getline(line, 8192);
        token = strtok(line, "\t");
        if (token == NULL) continue;
        string chrom5p = token;
        int pos5p = atoi(strtok(NULL, "\t"));
        char strand5p = strtok(NULL, "\t")[0];
        string chrom3p = strtok(NULL, "\t");
        int pos3p = atoi(strtok(NULL, "\t"));
        char strand3p = strtok(NULL, "\t")[0];
        int support = atoi(strtok(NULL, "\t"));
        // cout << chrom5p << " " << pos5p << " " << strand5p << " " << chrom3p << " " << pos3p << " " << strand3p << " " << support << endl;
        if (support > 0) {
            // cout << support << endl;
            this->insertRecord(chrom5p, pos5p, strand5p,
                               chrom3p, pos3p, strand3p,
                               support);
        }
    }
    dbFile.close();
    cout << "Read done" << endl;
}

// void JunctionDB::insertRecord(int aPos_5p, char aStrand_5p, int aPos_3p, char aStrand_3p, int aSupport) {
//     vector<int>::iterator iterPos = lower_bound(mRecordPos->begin(), mRecordPos->end(), aPos_5p);
//     vector<Record *>::iterator iterRecord = iterPos - mRecordPos->begin() + mRecords->begin();
//     
//     if (iterPos == mRecordPos->end() || *iterPos != aPos_5p) {
//         // new entry
//         Record *new_record = new Record();
//         new_entry->pos = aPos_5p;
//         new_entry->from_plus_strand_pos = new vector<int>();
//         new_entry->from_plus_strand_sign = new vector<char>();
//         new_entry->from_plus_strand_support = new vector<int>();
//         new_entry->from_minus_strand_pos = new vector<int>();
//         new_entry->from_minus_strand_sign = new vector<char>();
//         new_entry->from_minus_strand_support = new vector<int>();
//         if (aStrand_5p == '+') {
//             new_entry->from_plus_strand_pos->push_back(aPos_3p);
//             new_entry->from_plus_strand_sign->push_back(aStrand_3p);
//             new_entry->from_plus_strand_support->push_back(aSupport);
//         } else {
//             new_entry->from_minus_strand_pos->push_back(aPos_3p);
//             new_entry->from_minus_strand_sign->push_back(aStrand_3p);
//             new_entry->from_minus_strand_support->push_back(aSupport);
//         }
// 
//         mEntries->insert(iterEntry, new_entry);
//         mEntryId->insert(iterId, aPos_5p);
//     } else {
//         // exists
//         if (aStrand_5p == '+') {
//             (*iterEntry)->from_plus_strand_pos->push_back(aPos_3p);
//             (*iterEntry)->from_plus_strand_sign->push_back(aStrand_3p);
//             (*iterEntry)->from_plus_strand_support->push_back(aSupport);
//         } else {
//             (*iterEntry)->from_minus_strand_pos->push_back(aPos_3p);
//             (*iterEntry)->from_minus_strand_sign->push_back(aStrand_3p);
//             (*iterEntry)->from_minus_strand_support->push_back(aSupport);
//         }
//     }
// }

void JunctionDB::insertRecord(string aChrom_5p, int aPos_5p, char aStrand_5p,
                              string aChrom_3p, int aPos_3p, char aStrand_3p,
                              int aSupport) {
    int id = (aStrand_5p == '+') ? aPos_5p : -aPos_5p;
    string uuid = aChrom_5p + ':' + to_string(id);
    vector<string>::iterator iterUUID = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    vector<Record *>::iterator iterRec = iterUUID - mRecordUUID->begin() + mRecords->begin();

    // for first segment
    if (iterUUID == mRecordUUID->end() || *iterUUID != uuid) {
        // new record
        Record *rec = new Record(aChrom_5p, aPos_5p, aStrand_5p);
        mRecordUUID->insert(iterUUID, uuid);
        mRecords->insert(iterRec, rec);
        rec->insertForwardEntry(aChrom_3p, aPos_3p, aStrand_3p, aSupport, false);
    } else {
        // exists
        (*iterRec)->insertForwardEntry(aChrom_3p, aPos_3p, aStrand_3p, aSupport, false);
    }
    id = -id;
    uuid = aChrom_5p + ':' + to_string(id);
    iterUUID = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    iterRec = iterUUID - mRecordUUID->begin() + mRecords->begin();
    if (iterUUID == mRecordUUID->end() || *iterUUID != uuid) {
        // new record
        Record *rec = new Record(aChrom_5p, aPos_5p, (aStrand_5p == '+') ? '-' : '+');
        mRecordUUID->insert(iterUUID, uuid);
        mRecords->insert(iterRec, rec);
        rec->insertBackwardEntry(aChrom_3p, aPos_3p, (aStrand_3p == '+') ? '-' : '+', aSupport, true);
    } else {
        // exists
        (*iterRec)->insertBackwardEntry(aChrom_3p, aPos_3p, (aStrand_3p == '+') ? '-' : '+', aSupport, true);
    }

    // for second segment
    id = (aStrand_3p == '+') ? aPos_3p : -aPos_3p;
    uuid = aChrom_3p + ':' + to_string(id);
    iterUUID = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    iterRec = iterUUID - mRecordUUID->begin() + mRecords->begin();
    if (iterUUID == mRecordUUID->end() || *iterUUID != uuid) {
        // new record
        Record *rec = new Record(aChrom_3p, aPos_3p, aStrand_3p);
        mRecordUUID->insert(iterUUID, uuid);
        mRecords->insert(iterRec, rec);
        rec->insertBackwardEntry(aChrom_5p, aPos_5p, aStrand_5p, aSupport, false);
    } else {
        (*iterRec)->insertBackwardEntry(aChrom_5p, aPos_5p, aStrand_5p, aSupport, false);
    }
    id = -id;
    uuid = aChrom_3p + ':' + to_string(id);
    iterUUID = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    iterRec = iterUUID - mRecordUUID->begin() + mRecords->begin();
    if (iterUUID == mRecordUUID->end() || *iterUUID != uuid) {
        // new record
        Record *rec = new Record(aChrom_3p, aPos_3p, (aStrand_3p == '+') ? '-' : '+');
        mRecordUUID->insert(iterUUID, uuid);
        mRecords->insert(iterRec, rec);
        rec->insertForwardEntry(aChrom_5p, aPos_5p, (aStrand_5p == '+') ? '-' : '+', aSupport, true);
    } else {
        (*iterRec)->insertForwardEntry(aChrom_5p, aPos_5p, (aStrand_5p == '+') ? '-' : '+', aSupport, true);
    }
}

void JunctionDB::sortRecordEntry() {
    for (Record *r : *mRecords) {
        r->sortEntry();
    }
}

Record *JunctionDB::findRecord(string aChrom, int aPos, char aStrand) {
    int id = (aStrand == '+') ? aPos : -aPos;
    string uuid = aChrom + ':' + to_string(id);
    // int uuid = (aStrand == '+') ? aPos : -aPos;
    vector<string>::iterator iterId = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    vector<Record *>::iterator iterRecord = iterId - mRecordUUID->begin() + mRecords->begin();
    if (iterId == mRecordUUID->end() || *iterId != uuid) {
        return NULL;
    } else {
        return *iterRecord;
    }
    // int uuid = (aStrand == '+') ? aId : -aId;
    // vector<int>::iterator iterUUID = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    // if (iterUUID == mRecordUUID->end() || *iterUUID != uuid) {
    //     return NULL;
    // } else {
    //     return *(iterUUID - mRecordUUID->begin() + mRecords->begin());
    // }
}

vector<Record *> *JunctionDB::findRecords(string aChrom, int aPos, char aStrand) {
    auto *recs = new vector<Record *>();
    int id = (aStrand == '+') ? aPos : -aPos;
    string uuid = aChrom + ':' + to_string(id);
    // int uuid = (aStrand == '+') ? aPos : -aPos;
    auto iterId = lower_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    auto upperId = upper_bound(mRecordUUID->begin(), mRecordUUID->end(), uuid);
    auto iterRecord = iterId - mRecordUUID->begin() + mRecords->begin();
    if (iterId == mRecordUUID->end() || *iterId != uuid) {
        return NULL;
    } else {
        for (auto it = iterRecord; it != mRecords->end(); it++) {
            recs->push_back(*it);
        }
        return recs;
    }
}
// void JunctionDB::updateRecordsFromFile(const char * aAbnormalFn) {
//     ifstream fin(aAbnormalFn);
//     if (!fin) {
//         cerr << "Cannot open file: " << aAbnormalFn << endl;
//         exit(1);
//     }
//     string line;
//     int leftId, rightId;
//     char leftStrand, rightStrand;
//     double support;
//     while (getline(fin, line)) {
//         cout << line << endl;
//         stringstream ss(line);
//         ss >> leftId >> leftStrand >> rightId >> rightStrand >> support;
//         this->insertRecord(leftId, leftStrand, rightId, rightStrand, support);
//     }
//     fin.close();
// }

// void JunctionDB::writeDB(const char * outFn) {
//     ofstream fout(outFn);
//     for (int i = 0; i < mRecordUUID->size(); i++) {
//         for (entry * e : *((*mRecords)[i]->getForwardEntries())) {
//             if (!e->isComplement) fout << (*mRecords)[i]->getId() << " " << (char)(*mRecords)[i]->getStrand() << " " << e->id << " " << e->strand << " " << e->support << endl;
//         }
//     }
//     fout.close();
// }

// void JunctionDB::print() {
//     for (junc_entry_t *r : *mEntries) {
//         for (int i = 0; i < r->from_plus_strand_pos->size(); i++) {
//             cout << r->pos << " + " << (*r->from_plus_strand_pos)[i] << " " << (*r->from_plus_strand_sign)[i] << " " << (*r->from_plus_strand_support)[i] << endl;
//         }
//         for (int i = 0; i < r->from_minus_strand_pos->size(); i++) {
//             cout << r->pos << " - " << (*r->from_minus_strand_pos)[i] << " " << (*r->from_minus_strand_sign)[i] << " " << (*r->from_minus_strand_support)[i] << endl;
//         }
//     }
// }

void JunctionDB::print() {
    cout << "Junction database: " << endl;
    for (Record *r : *mRecords) {
        vector<entry_t *> *eForward = r->getForwardEntries();
        vector<entry_t *> *eBackward = r->getBackwardEntries();
        for (entry_t *e : *eForward) {
            cout << r->getChrom() << ':' << r->getPos() << (char) (r->getStrand()) << " => " << e->chrom << ':'
                 << e->pos << e->strand << " " << e->support << endl;
        }
        cout << "forward done" << endl;
        for (entry_t *e : *eBackward) {
            cout << r->getChrom() << ':' << r->getPos() << (char) (r->getStrand()) << " <= " << e->chrom << ':'
                 << e->pos << e->strand << " " << e->support << endl;
        }
        cout << "backward done" << endl;
    }
}

// void JunctionDB::printBps() {
//     for (vector<int>::iterator i = mBps->begin(); i != mBps->end(); i++) {
//         cout << *i;
//         if (i != mBps->begin() && (i - mBps->begin()) % 10 == 0) cout << endl; else cout << " ";
//     }
//     cout << endl;
// }
