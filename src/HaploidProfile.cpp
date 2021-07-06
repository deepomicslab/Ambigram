#include <iostream>
#include <fstream>
#include <sstream>

#include "HaploidProfile.hpp"

using namespace std;

// Strand::Strand(SignedSeg * aSeg, char aSign) {
//     belongSeg = aSeg;
//     sign = aSign;
//     if (aSign == '+') {
//         aSeg->sp = this;
//     } else {
//         aSeg->sn = this;
//     }
// }

// SignedSeg::SignedSeg(seg * aSeg) {
//     oriSeg = aSeg;
//     sp = sn = NULL;
// }

HaploidProfile::HaploidProfile(string aSampleName) {
    mSampleName = aSampleName;
    mHaploid1 = new vector<Strand *>();
    mHaploid2 = new vector<Strand *>();
}

HaploidProfile::~HaploidProfile() {
    delete[] mHaploid1;
    delete[] mHaploid2;
}

void HaploidProfile::readHaploids(const char *aFilename) {
    ifstream fin(aFilename);
    if (!fin) {
        cerr << "Cannot open file: " << aFilename << endl;
        exit(1);
    }
    string line;
    stringstream ss;

    // first
    getline(fin, line);
    ss = stringstream(line);
    string node;

    while (ss >> node) {
        int id = stoi(node.substr(0, node.size() - 1));
        char sign = node.back();
        // SignedSeg * seg = new SignedSeg(id - 1, mChr, mStart - 1, mEnd - 1);
        Strand *strand;
        if (sign == '+') {
            strand = new Strand{(*(mSegRef->getSegs()))[id - 1], '+'};
        } else {
            strand = new Strand{(*(mSegRef->getSegs()))[id - 1], '-'};
        }
        mHaploid1->push_back(strand);
    }
    mHaploid1->pop_back();

    // second
    getline(fin, line);
    ss = stringstream(line);
    while (ss >> node) {
        int id = stoi(node.substr(0, node.size() - 1));
        char sign = node.back();
        // SignedSeg * seg = new SignedSeg(id - 1, mChr, mStart - 1, mEnd - 1);
        Strand *strand;
        if (sign == '+') {
            strand = new Strand{(*(mSegRef->getSegs()))[id - 1], '+'};
        } else {
            strand = new Strand{(*(mSegRef->getSegs()))[id - 1], '-'};
        }
        mHaploid2->push_back(strand);
    }
    mHaploid2->pop_back();
    fin.close();
}

void HaploidProfile::setSegRef(SegmentDB *db) {
    mSegRef = db;
    mSegNormal = new bool[db->getSegs()->size()]{false};
}

void HaploidProfile::identifyNormal() {
    int *segCount1 = new int[mSegRef->getSegs()->size()]{0};
    int *segCount2 = new int[mSegRef->getSegs()->size()]{0};
    for (Strand *s : *mHaploid1) {
        segCount1[s->belongSeg->id]++;
    }
    for (Strand *s : *mHaploid2) {
        segCount2[s->belongSeg->id]++;
    }

    for (int i = 0; i < mSegRef->getSegs()->size(); i++) {
        if (segCount1[i] == 1 && segCount2[i] == 1) {
            // normal seg
            mSegNormal[i] = true;
        }
    }
}

void HaploidProfile::setSupportProfile(SupportProfile *sp) {
    mSp = sp;
}

void HaploidProfile::placeVariantsInSeg(seg *aSeg) {
    // only for normal segments
    int nHom = 0;
    int nHet = 0;
    int nUnknown = 0;

    int solvedHet = 0;
    int noSupport = 0;
    vector<Variant *> hap;
    if (mSegNormal[aSeg->id]) {
        for (int i = 0; i < aSeg->loci.size(); i++) {
            locus *l = aSeg->loci[i];
            int gt = mSp->getGT()[l->id];
            if (gt != 1) {
                // hom, not related supports
                if (gt == 0) {
                    hap.push_back(new Variant{l, 0});
                    nHom++;
                } else if (gt == 2) {
                    hap.push_back(new Variant{l, 1});
                    nHom++;
                } else {
                    hap.push_back(new Variant{l, -2});
                    nUnknown++;
                }
            } else {
                nHet++;
                vector<locus *> sameSegLoci;
                vector<readCount *> sameSegCounts;
                mSp->getInSameSegSupports(l, sameSegLoci, sameSegCounts);
                if (sameSegLoci.size() == 0) {
                    noSupport++;
                    hap.push_back(new Variant{l, -2});
                    continue;
                }
                solvedHet++;
                int hap0support = 0;
                int hap1support = 0;
                for (int j = 0; j < sameSegLoci.size(); j++) {
                    if (sameSegLoci[j]->id < l->id) {
                        Variant *pv = hap[sameSegLoci[j]->id - hap.front()->l->id];
                        if (pv->type == 0) {
                            hap0support += sameSegCounts[j]->rr + sameSegCounts[j]->aa;
                            hap1support += sameSegCounts[j]->ra + sameSegCounts[j]->ar;
                        } else {
                            hap0support += sameSegCounts[j]->ra + sameSegCounts[j]->ar;
                            hap1support += sameSegCounts[j]->rr + sameSegCounts[j]->aa;
                        }
                    } else {
                        if (hap.size() == 0) {
                            hap.push_back(new Variant{l, 0});
                        }
                        break;
                    }
                    // cout << aSeg->id << " " << l->id << " " << sameSegLoci[j]->id << " " << sameSegCounts[j]->rr << " " << sameSegCounts[j]->ra << " " << sameSegCounts[j]->ar << " " << sameSegCounts[j]->aa << endl;
                }

                if (hap0support >= hap1support) {
                    hap.push_back(new Variant{l, 0});
                } else {
                    hap.push_back(new Variant{l, 1});
                }
            }
        }
        // if (hap.size() > 0) {
        //     cout << aSeg->id << " " << aSeg->loci.size() << endl;
        //     for (Variant * v : hap) {
        //         // if (mSp->getGT()[v->l->id] == 1) {
        //             cout << v->l->id << "::" << v->type << " ";
        //         //}
        //     }
        //     cout << endl;
        // }
        cout << aSeg->id << " " << nHom << " " << nHet << " " << nUnknown << " " << solvedHet << " " << solvedHet + nHom
             << " " << noSupport << " " << aSeg->loci.size() << " "
             << ((aSeg->loci.size() > 0) ? solvedHet * 1.0 / aSeg->loci.size() * 100 : 0) << " "
             << ((aSeg->loci.size() > 0) ? (nHom + solvedHet) * 1.0 / aSeg->loci.size() * 100 : 0) << " "
             << ((nHet > 0) ? solvedHet * 1.0 / nHet * 100 : 0) << endl;
    }
}

void HaploidProfile::placeVariants() {
    cout
            << "segId hom het unknown solvedHet hom+solvedHet noSupportHet total solvedHet/total hom_solvedHet/total solvedHet/het"
            << endl;
    for (int i = 0; i < mSegRef->getSegs()->size(); i++) {
        if (mSegNormal[i]) {
            // TODO get loci and supports in the same segment then use greedy to select;
            this->placeVariantsInSeg((*(mSegRef->getSegs()))[i]);
        }
    }
}

void HaploidProfile::print() {
    cout << "Haploid 1: " << endl;
    for (Strand *s : *mHaploid1) {
        cout << s->belongSeg->id << s->sign << " ";
    }
    cout << endl;
    cout << "Haploid 2: " << endl;
    for (Strand *s : *mHaploid2) {
        cout << s->belongSeg->id << s->sign << " ";
    }
    cout << endl;
}

void HaploidProfile::printNormal() {
    for (int i = 0; i < mSegRef->getSegs()->size(); i++) {
        seg *s = (*(mSegRef->getSegs()))[i];
        if (mSegNormal[i]) {
            if (s->loci.size() <= 1) continue;
            for (int j = 0; j < s->loci.size(); j++) {
                locus *l = (*(mSp->getLociSupports()))[s->loci[j]->id].first;
                support *sp = (*(mSp->getLociSupports()))[s->loci[j]->id].second;
                for (int k = 0; k < sp->pairedLoci.size(); k++) {
                    if (sp->pairedLoci[k]->belongSeg != l->belongSeg) continue;
                    readCount *rc = sp->pairedCounts[k];
                    cout << s->id << " " << s->loci[j]->id << " " << sp->pairedLoci[k]->id << " " << rc->rr << " "
                         << rc->ra << " " << rc->ar << " " << rc->aa << endl;
                }
            }
        }
    }
    cout << endl;
}
