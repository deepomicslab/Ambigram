#include <iostream>

#include "SegmentDB.hpp"
#include "LocusDB.hpp"
#include "SupportProfile.hpp"
#include "HaploidProfile.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    const char * samplename = argv[1];
    const char * chr = argv[2];
    int start = atoi(argv[3]);
    int end = atoi(argv[4]);
    const char * haploidFn = argv[5];
    const char * segdbFn = argv[6];
    const char * vcfFn = argv[7];
    const char * supportFn = argv[8];

    SegmentDB * segdb = new SegmentDB(chr, start, end);
    segdb->readSegs(segdbFn);

    LocusDB * locidb = new LocusDB(chr, start, end);
    locidb->read(vcfFn, LocusDB::MODE_VCF);
    locidb->setSegRef(segdb);
    locidb->assignLocusToSeg();

    SupportProfile * sp = new SupportProfile(samplename);
    sp->setLociRef(locidb);
    sp->readGenotypes(vcfFn);
    sp->readSupport(supportFn);
    // sp->print();
    // return 0;

    HaploidProfile * hp = new HaploidProfile(samplename);
    hp->setSegRef(segdb);
    hp->readHaploids(haploidFn);
    // hp->print();
    hp->identifyNormal();
    hp->setSupportProfile(sp);
    hp->placeVariants();
    sp->printStatistics();
    // hp->printNormal();
    return 0;
}
