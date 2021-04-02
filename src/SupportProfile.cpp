#include <iostream>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <zlib.h>

#include "htslib/synced_bcf_reader.h"
#include "SupportProfile.hpp"

using namespace std; 

SupportProfile::SupportProfile(string aSampleName) {
    mSampleName = aSampleName;
    // mCoveredLoci = new vector< vector< pair<locus *, int> > *>();
    mCoveredLociSupports = new vector< pair<locus *, support *> >();
    // mSupportCount = new vector<int>();
    // mSupportIdentifiers = new vector<string>();

    mNumHet = mNumHom = mNumUnknown = 0;
}

SupportProfile::~SupportProfile() {
    // delete mCoveredLoci;
}

string SupportProfile::getSampleName() { return mSampleName; }
vector< pair<locus *, support *> > * SupportProfile::getLociSupports() { return mCoveredLociSupports; }
int * SupportProfile::getGT() { return mGT; }

void SupportProfile::setLociRef(LocusDB * db) {
    mLociRef = db;
    mGT = new int[db->getLoci()->size()]{-1};
}

void SupportProfile::readGenotypes(const char * aVCF) {
    bcf_srs_t * sr = bcf_sr_init();
    bcf_sr_set_regions(sr, (mLociRef->getChr() + ':' + to_string(mLociRef->getStart()) + '-' + to_string(mLociRef->getEnd())).c_str(), 0);
    bcf_sr_add_reader(sr, aVCF);
    bcf_hdr_t * header = bcf_sr_get_header(sr, 0);
    bcf_hdr_set_samples(header, mSampleName.c_str(), 0);

    bcf1_t * rec = bcf_init1();
    int * gt = NULL, ngt = 0;
    int idx = 0;
    while (bcf_sr_next_line(sr)) {
        rec = bcf_sr_get_line(sr, 0);
        ngt = bcf_get_genotypes(header, rec, &gt, &ngt);
        support * sp = new support();
        mGT[idx] = bcf_gt_allele(gt[0]) + bcf_gt_allele(gt[1]);
        if (mGT[idx] == 0 || mGT[idx] == 2) {
            mNumHom++;
        } else if (mGT[idx] == 1) {
            mNumHet++;
        } else {
            mNumUnknown++;
        }
        mCoveredLociSupports->push_back(make_pair((*(mLociRef->getLoci()))[idx], sp));
        idx++;
    }
    bcf_sr_destroy(sr);
}

void SupportProfile::readSupport(const char * aSupportFn) {
    gzFile fin = gzopen(aSupportFn, "rb");
    if (!fin) {
        cerr << "Cannot open file: " << aSupportFn << endl;
        exit(1);
    }

    char line[8192];
    while (gzgets(fin, line, 8192) != NULL) {
        stringstream ss(line);
        int id1, id2;
        int rr, ra, ar, aa;
        ss >> id1 >> id2 >> rr >> ra >> ar >> aa;
        
        support * sp1 = (*mCoveredLociSupports)[id1 - 1].second;
        locus * l2 = (*mCoveredLociSupports)[id2 - 1].first;
        sp1->pairedLoci.push_back(l2);
        sp1->pairedCounts.push_back(new readCount{rr, ra, ar, aa});
    }
    gzclose(fin);
}

void SupportProfile::countSupport(const char * aBamFn) {
    samFile * bam = sam_open(aBamFn, "rb");
    bam_hdr_t * header = sam_hdr_read(bam);
    bam1_t * aln = bam_init1();
    vector<bam1_t *> alnStack;
    int numDigitLociSize = log10((double)mLociRef->getLoci()->size()) + 1;

    int count = 0;
    while (sam_read1(bam, header, aln) >= 0) {
        count++;
        alnStack.push_back(bam_dup1(aln));
        while (sam_read1(bam, header, aln) >= 0) {
            count++;
            if (count % 50000 == 0) cout << count << endl;
            if (strcmp(bam_get_qname(aln), bam_get_qname(alnStack[0])) == 0) {
                alnStack.push_back(bam_dup1(aln));
            } else {
                break;
            }
        }
        
        vector< pair<locus *, int> > * coveredLoci = new vector < pair<locus *, int> >();
        for (bam1_t * aln_i : alnStack) {
            if ((int)aln_i->core.qual < 20 || (aln_i->core.flag & 0x900) != 0) continue;
            uint8_t * seqi = bam_get_seq(aln_i);
            vector<locus *>::const_iterator begin, end;
            mLociRef->findLociInRange(aln_i->core.pos, aln_i->core.pos + bam_cigar2rlen(aln_i->core.n_cigar, bam_get_cigar(aln_i)), begin, end);
            while (begin != end) {
                int gt = mGT[(*begin)->id];
                if (gt < 0 || gt % 2 == 0) {
                    begin++;
                    continue;
                }
                // cout << (*begin)->pos << " " << aln_i->core.pos << " " << (*begin)->ref << " " << bam_get_qname(aln_i) << " ";
                // uint32_t * cigar = bam_get_cigar(aln_i);
                // for (int i = 0; i < aln_i->core.n_cigar; i++) {
                //     int opchr = bam_cigar_opchr(cigar[i]);
                //     int oplen = bam_cigar_oplen(cigar[i]);
                //     cout << (char)opchr << oplen;
                // }
                // cout << " ";
                // for (int i = 0; i < aln_i->core.l_qseq; i++) {
                //     cout << seq_nt16_str[bam_seqi(seqi, i)];
                // }
                // cout << " ";
                int baseIdx = this->getBaseIdx(aln_i, (*begin)->pos);
                if (baseIdx < 0) {
                    begin++;
                    // cout << endl;
                    continue;
                }
                // cout << baseIdx << " ";
                // TODO low base quality may accur here, need to remove?
                char base = string(1, seq_nt16_str[bam_seqi(seqi, baseIdx)])[0];
                // cout << base << " " << (int)bam_get_qual(aln_i)[baseIdx] << " ";
                if (base == (*begin)->ref) {
                    coveredLoci->push_back(make_pair(*begin, 0));
                } else if (base == (*begin)->alt) {
                    coveredLoci->push_back(make_pair(*begin, 1));
                } else {
                    // cout << "NA";
                }
                // cout << endl;
                begin++;
            }
        }

        if (coveredLoci->size() > 1) {
            // sort(coveredLoci.begin(), coveredLoci.end(), [](pair<locus *, int> p1, pair<locus *, int> p2) { return p1.first->pos < p2.first->pos; });
            for (int i = 0; i < coveredLoci->size(); i++) {
                support * sp = (*mCoveredLociSupports)[(*coveredLoci)[i].first->id].second;
                for (int j = 0; j < coveredLoci->size(); j++) {
                    if (j == i) continue;
                    vector<locus *>::iterator iterPair = lower_bound(sp->pairedLoci.begin(), sp->pairedLoci.end(), (*coveredLoci)[j].first);
                    vector<readCount *>::iterator iterCount = iterPair - sp->pairedLoci.begin() + sp->pairedCounts.begin();
                    if (iterPair == sp->pairedLoci.end() || *iterPair != (*coveredLoci)[j].first) {
                        sp->pairedLoci.insert(iterPair, (*coveredLoci)[j].first);
                        int comb = (*coveredLoci)[i].second + (*coveredLoci)[j].second;
                        if (comb == 0) {
                            sp->pairedCounts.insert(iterCount, new readCount{1, 0, 0, 0});
                        } else if (comb == 1) {
                            if ((*coveredLoci)[i].second == 0) {
                                sp->pairedCounts.insert(iterCount, new readCount{0, 1, 0, 0});
                            } else {
                                sp->pairedCounts.insert(iterCount, new readCount{0, 0, 1, 0});
                            }
                        } else {
                            sp->pairedCounts.insert(iterCount, new readCount{0, 0, 0, 1});
                        }
                    } else {
                        int comb = (*coveredLoci)[i].second + (*coveredLoci)[j].second;
                        if (comb == 0) {
                            (*iterCount)->rr++;
                        } else if (comb == 1) {
                            if ((*coveredLoci)[i].second == 0) {
                                (*iterCount)->ra++;
                            } else {
                                (*iterCount)->ar++;
                            }
                        } else {
                            (*iterCount)->aa++;
                        }
                    }
                }
            }
        }

        alnStack.clear();
        alnStack.push_back(bam_dup1(aln));
    }
}

void SupportProfile::writeSupport(const char * outFn) {
    gzFile fout = gzopen(outFn, "wb");
    for (int i = 0; i < mCoveredLociSupports->size(); i++) {
        stringstream ss;
        locus * l = (*mCoveredLociSupports)[i].first;
        support * s = (*mCoveredLociSupports)[i].second;
        for (int j = 0; j < s->pairedLoci.size(); j++) {
            readCount * rc = s->pairedCounts[j];
            ss << l->id + 1 << " " << s->pairedLoci[j]->id + 1 << " " << rc->rr << " " << rc->ra << " " << rc->ar << " " << rc->aa << endl;
        }
        gzputs(fout, ss.str().c_str());
    }
    gzclose(fout);
}


void SupportProfile::getInSameSegSupports(locus * l, vector<locus *> & pLoci, vector<readCount *> & pCounts) {
    pLoci.clear();
    pCounts.clear();

    support * sp = (*mCoveredLociSupports)[l->id].second;
    for (int i = 0; i < sp->pairedLoci.size(); i++) {
        if (sp->pairedLoci[i]->belongSeg == l->belongSeg && l != sp->pairedLoci[i]) {
            pLoci.push_back(sp->pairedLoci[i]);
            pCounts.push_back(sp->pairedCounts[i]);
        }
    }
}

int SupportProfile::getBaseIdx(bam1_t * aln, int pos) {
    uint32_t * cigar = bam_get_cigar(aln);
    int alnStart = aln->core.pos;
    int op, oplen;
    int idx = 0;
    for (int i = 0; i < aln->core.n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH) {
            if (alnStart + oplen < pos) {
                alnStart += oplen;
                idx += oplen;
            } else {
                idx += pos - alnStart;
                if (idx >= bam_cigar2rlen(aln->core.n_cigar, cigar)) {
                    return -1;
                }
                return idx;
            }
        } else if (op == BAM_CDEL) {
            if (alnStart + oplen < pos) {
                alnStart += oplen;
            } else {
                return -1;
            }
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            idx += oplen;
        }
    }
    return -1;
}

void SupportProfile::print() {
    for (int i = 0; i < mCoveredLociSupports->size(); i++) {
        locus * l = (*mCoveredLociSupports)[i].first;
        support * s = (*mCoveredLociSupports)[i].second;
        for (int j = 0; j < s->pairedLoci.size(); j++) {
            readCount * rc = s->pairedCounts[j];
            cout << l->id + 1 << " " << s->pairedLoci[j]->id + 1 << " " << rc->rr << " " << rc->ra << " " << rc->ar << " " << rc->aa << endl;
        }
    }
}

void SupportProfile::printStatistics() {
    cout << mNumHom << " " << mNumHet << " " << mNumUnknown << " " << mLociRef->getLoci()->size() << endl;
}
