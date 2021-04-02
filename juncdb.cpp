#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include "SVprofile.hpp"
#include "SegmentDB.hpp"
#include "JunctionDB.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
    if (strcmp(argv[1], "segs") == 0) {
        const char * rawSVListFn = argv[2];
        const char * chr = argv[3];
        int start = atoi(argv[4]);
        int end = atoi(argv[5]);
        const char * outSegsFn = argv[6];

        SegmentDB * segdb = new SegmentDB(chr, start - 1, end - 1);
        ifstream rawSVList(rawSVListFn);
        if (!rawSVList) {
            cerr << "Cannot open file: " << argv[2] << endl;
            exit(1);
        }
        string line;
        string samplename, path;
        while (getline(rawSVList, line)) {
            stringstream ss(line);
            ss >> samplename >> path;
            SVprofile * sv = new SVprofile(path.c_str(), samplename);
            sv->filterAbnormal(chr, start - 1, end - 1);
            segdb->updateBps(sv);
        }
        rawSVList.close();

        segdb->constructSegsFromBps();
        segdb->writeSegs(outSegsFn);
    } else if (strcmp(argv[1], "indv") == 0) {
        const char * juncdbSegsFn = argv[2];
        const char * indvRawSVFn = argv[3];
        const char * indvBamFn = argv[4];
        const char * indvDepthFn = argv[5];
        const char * chr = argv[6];
        int start = atoi(argv[7]);
        int end = atoi(argv[8]);
        const char * outLocalhapFn = argv[9];
        const char * outNormalJuncFn = argv[10];
        const char * outAbnormalJuncFn = argv[11];
        const char * samplename = argv[12];
        
        SegmentDB * segdb = new SegmentDB(chr, start, end);
        JunctionDB * juncdb = new JunctionDB(chr, start, end);
        cout << "Reading segments..." << endl;
        segdb->readSegs(juncdbSegsFn);
        cout << "Loading SV profile..." << endl;
        SVprofile * sv = new SVprofile(indvRawSVFn, samplename);
        cout << "Filtering SV profile..." << endl;
        sv->filterAbnormal(chr, start - 1, end - 1);
        cout << "Seting segments..." << endl;
        sv->setSegDB(segdb);
        cout << "Converting position to ID..." << endl;
        sv->pos2id();
        cout << "Couting segment depth..." << endl;
        sv->countSegDepth(indvDepthFn);
        cout << "Couting normal junction..." << endl;
        sv->countNormal(indvBamFn);

        sv->writeLocalHap(outLocalhapFn);
        sv->writeNormal(outNormalJuncFn);
        sv->writeAbnormal(outAbnormalJuncFn);
    } else if (strcmp(argv[1], "create") == 0) {
        const char * normalListFn = argv[2];
        const char * abnormalListFn = argv[3];
        const char * dbFn = argv[4];
        const char * chr = argv[5];
        int start = atoi(argv[6]);
        int end = atoi(argv[7]);

        JunctionDB * juncdb = new JunctionDB(chr, start, end);
        ifstream normalList(normalListFn);
        ifstream abnormalList(abnormalListFn);
        string line;
        while (getline(normalList, line)) {
            juncdb->updateRecordsFromFile(line.c_str());
        }
        while (getline(abnormalList, line)) {
            juncdb->updateRecordsFromFile(line.c_str());
        }
        abnormalList.close();
        juncdb->sortRecordEntry();
        juncdb->writeDB(dbFn);
    }
    return 0;
}
