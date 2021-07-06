#ifndef _RECORD_H_
#define _RECORD_H_

#include <vector>
#include <tuple>

using namespace std;

typedef struct Entry {
    string chrom;
    int pos;
    char strand;
    int support;
    bool isComplement;
} entry_t;

// typedef struct JuncEntry {
//     int pos;
//     vector<int> *from_plus_strand_pos;
//     vector<char> *from_plus_strand_sign;
//     vector<int> *from_plus_strand_support;
//     vector<int> *from_minus_strand_pos;
//     vector<char> *from_minus_strand_sign;
//     vector<int> *from_minus_strand_support;
// } junc_entry_t;

class Record {
protected:
    string mChrom;
    int mPos;
    char mStrand;
    // vector<int> * mPositivePairId;
    // vector<entry *> * mBackwardEntries;
    // vector<int> * mBackwardEntryUUID;
    vector<entry_t *> *mBackwardEntries;
    vector<string> *mBackwardEntryUUID;

    // vector<int> * mNegativePairId;
    // vector<entry *> * mForwardEntries;
    // vector<int> * mForwardEntryUUID;
    vector<entry_t *> *mForwardEntries;
    vector<string> *mForwardEntryUUID;

    Record *mComplementRecord;

public:
    Record(string aChrom, int aPos, char aStrand);

    ~Record();

    string getChrom();

    int getPos();

    char getStrand();

    vector<string> *getBackwardEntryUUID();

    vector<string> *getForwardEntryUUID();

    vector<entry_t *> *getBackwardEntries();

    vector<entry_t *> *getForwardEntries();

    void insertBackwardEntry(string aChrom, int aPos, char aStrand, int aSupport, bool aIsComplement);

    void insertForwardEntry(string aChrom, int aPos, char aStrand, int aSupport, bool aIsComplement);

    void sortEntry();

    entry_t *findForwardEntry(string aChrom, int aPos, char aStrand);

    void print();
};

#endif
