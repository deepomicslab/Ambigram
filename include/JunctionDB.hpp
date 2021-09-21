#ifndef _JUNCTION_DB_H_
#define _JUNCTION_DB_H_

#include "Record.hpp"
#include "SVprofile.hpp"
#include "Junction.hpp"

using namespace std;

class Record;

class SVprofile;

class JunctionDB {
protected:
    // string mChr;
    // int mStart;
    // int mEnd;

    // vector<junc_entry_t *> *mEntries;
    // vector<int> *mEntryId;
    vector<Record *> *mRecords;
    vector<string> *mRecordUUID;

public:
    // JunctionDB(string chr, int start, int end);
    JunctionDB(const char *aFilename);
    JunctionDB(vector<Junction *> &junctions);

    ~JunctionDB();

    vector<Record *> *getRecords();

    void readDB(const char *aFilename);

    // void insertRecord(int aSourceId, char aSourceDir, int aTargetId, char aTargetDir, int aSupport);
    void insertRecord(string aChrom_5p, int aPos_5p, char aStrand_5p,
                      string aChrom_3p, int aPos_3p, char aStrand_3p,
                      int aSupport);

    void sortRecordEntry();

    Record *findRecord(string aChrom, int aPos, char aStrand);

    vector<Record *> *findRecords(string aChrom, int aPos, char aStrand);


    void updateRecordsFromFile(const char *aAbnormalFn);

    void writeDB(const char *outFn);

    void print();

    void printBps();
};

#endif
