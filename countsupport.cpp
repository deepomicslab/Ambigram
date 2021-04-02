#include <iostream>

#include "LocusDB.hpp"
#include "SupportProfile.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    const char * vcfFn = argv[1];
    const char * bamFn = argv[2];
    const char * samplename = argv[3];
    const char * outFn = argv[4];

    LocusDB * locidb = new LocusDB("chr6", 28460000, 33500000);
    locidb->read(vcfFn, LocusDB::MODE_VCF);

    SupportProfile * sp = new SupportProfile(samplename);
    sp->setLociRef(locidb);
    sp->readGenotypes(vcfFn);
    sp->countSupport(bamFn);
    // sp->print();
    sp->writeSupport(outFn);
    return 0;
}
