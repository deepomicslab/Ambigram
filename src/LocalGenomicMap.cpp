#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <unistd.h>
#include <deque>

#include "LocalGenomicMap.hpp"
#include "Exceptions.hpp"
#include "coin/CbcModel.hpp"
#include "coin/OsiClpSolverInterface.hpp"

using namespace std;

LocalGenomicMap::LocalGenomicMap(Graph * aGraph) {
    mGraph = aGraph;

    mCircuits = new vector< VertexPath *>();
    mHaploids = new vector< VertexPath *>();
}

LocalGenomicMap::~LocalGenomicMap() {;}

Graph * LocalGenomicMap::getGraph() { return mGraph; }
void LocalGenomicMap::setGraph(Graph * aGraph) { mGraph = aGraph; }

vector<VertexPath *> *LocalGenomicMap::get_long_frags() {
    return mLongFrags;
}

void LocalGenomicMap::read_long_frags(const char *fn) {
    mLongFrags = new vector<VertexPath *>();
    ifstream in(fn);
    string line, value;
    stringstream ss;
    while (getline(in, line)) {
        ss = stringstream(line);
        VertexPath *p = new VertexPath();
        while (ss >> value) {
            int id = stoi(value.substr(0, value.length() - 1));
            char strand = value.back();
            if (strand == '+') {
                p->push_back(mGraph->getSegmentById(id)->getPositiveVertex());
            } else {
                p->push_back(mGraph->getSegmentById(id)->getNegativeVertex());
            }
        }
        if (p->front()->getDir() == '-') {
            if (p->back()->getDir() == '-') {
                mLongFrags->push_back(this->get_complement(p));
            } else {
                if (p->back()->getId() < p->front()->getId()) {
                    mLongFrags->push_back(this->get_complement(p));
                }
            }
            p->clear();
        } else {
            mLongFrags->push_back(p);
        }
    }
    sort(mLongFrags->begin(), mLongFrags->end(), [](VertexPath *p1, VertexPath *p2) {
        int size = min(p1->size(), p2->size());
        for (int i = 0; i < size; i++) {
            if (p1->at(i)->getId() != p2->at(i)->getId()) {
                return p1->at(i)->getId() < p2->at(i)->getId();
            }
        }
    });

    in.close();
    
    vector<VertexPath *> *long_frags_new = new vector<VertexPath *>();
    vector<VertexPath *> *temp;
    while (!this->equal_frags(mLongFrags, long_frags_new)) {
        temp = mLongFrags;
        mLongFrags = long_frags_new;
        long_frags_new = temp;
        mLongFrags->clear();
        for (VertexPath *frag : *long_frags_new) {
            this->merge_long_frags(frag);
        }
        sort(mLongFrags->begin(), mLongFrags->end(), [](VertexPath *p1, VertexPath *p2) {
            int size = min(p1->size(), p2->size());
            for (int i = 0; i < size; i++) {
                if (p1->at(i)->getId() != p2->at(i)->getId()) {
                    return p1->at(i)->getId() < p2->at(i)->getId();
                }
            }
        });
    }
    sort(mLongFrags->begin(), mLongFrags->end(), [](VertexPath *p1, VertexPath *p2) {
        return p1->size() > p2->size();        
    });
}

void LocalGenomicMap::merge_long_frags(VertexPath *aFrag) {
    VertexPath::iterator found_iter, iter;
    bool found = false;
    for (iter = aFrag->end(); iter != aFrag->begin(); iter--) {
        for (int i = 0; i < mLongFrags->size(); i++) {
            found_iter = find_end(mLongFrags->at(i)->begin(), mLongFrags->at(i)->end(), aFrag->begin(), iter);
            if (found_iter != mLongFrags->at(i)->end()) {
                found = true;
                VertexPath *new_f = new VertexPath(mLongFrags->at(i)->begin(), found_iter);
                new_f->insert(new_f->end(), aFrag->begin(), aFrag->end());
                found_iter = find_end(mLongFrags->at(i)->begin(), mLongFrags->at(i)->end(), new_f->begin(), new_f->end());
                if (found_iter == mLongFrags->at(i)->end()) {
                    // new_f in
                //     break;
                // } else {
                    // new_f not in
                    found_iter = find_end(new_f->begin(), new_f->end(), mLongFrags->at(i)->begin(), mLongFrags->at(i)->end());
                    if (found_iter != new_f->end()) {
                        // new_f include
                        mLongFrags->at(i)->clear();
                        mLongFrags->at(i)->insert(mLongFrags->at(i)->end(), new_f->begin(), new_f->end());
                    } else {
                        // new_f not include
                        mLongFrags->push_back(new_f);
                    }
                }
                break;
            }
        }
        if (found) {
            break;
        }
    }
    if (!found) {
        VertexPath *comp = this->get_complement(aFrag);
        for (iter = comp->end(); iter != comp->begin(); iter--) {
            for (int i = 0; i < mLongFrags->size(); i++) {
                found_iter = find_end(mLongFrags->at(i)->begin(), mLongFrags->at(i)->end(), comp->begin(), iter);
                if (found_iter != mLongFrags->at(i)->end()) {
                    found = true;
                    VertexPath *new_f = new VertexPath(mLongFrags->at(i)->begin(), found_iter);
                    new_f->insert(new_f->end(), comp->begin(), comp->end());
                    found_iter = find_end(mLongFrags->at(i)->begin(), mLongFrags->at(i)->end(), new_f->begin(), new_f->end());
                    if (found_iter == mLongFrags->at(i)->end()) {
                        // new_f in
                    //     break;
                    // } else {
                        // new_f not in
                        found_iter = find_end(new_f->begin(), new_f->end(), mLongFrags->at(i)->begin(), mLongFrags->at(i)->end());
                        if (found_iter != new_f->end()) {
                            // new_f include
                            mLongFrags->at(i)->clear();
                            mLongFrags->at(i)->insert(mLongFrags->at(i)->end(), new_f->begin(), new_f->end());
                        } else {
                            // new_f not include
                            mLongFrags->push_back(new_f);
                        }
                    }
                    break;
                }
            }
            if (found) {
                break;
            }
        }
    }
    if (!found) {
        mLongFrags->push_back(aFrag);
    }
}

bool LocalGenomicMap::equal_frags(vector<VertexPath *> *f1, vector<VertexPath *> *f2) {
    if (f1->size() != f2->size()) {
        return false;
    }
    for (int i = 0; i < f1->size(); i++) {
        if (*f1->at(i) != *f2->at(i)) {
            return false;
        }
    }
    return true;
}

VertexPath *LocalGenomicMap::get_complement(VertexPath *path) {
    VertexPath *comp = new VertexPath();
    for (VertexPath::reverse_iterator rit = path->rbegin(); rit != path->rend(); rit++) {
        comp->push_back((*rit)->getComplementVertex());
    }
    return comp;
}

vector<double> LocalGenomicMap::scaleILPCoef(vector<double> aCovs) {
    double mean = accumulate(aCovs.begin(), aCovs.end(), 0.0) / aCovs.size();
    vector<double> diff(aCovs.size());
    transform(aCovs.begin(), aCovs.end(), diff.begin(),
              [mean](double v) {return v - mean + 1;});
    vector<double> diff_s(diff.size());
    transform(diff.begin(), diff.end(), diff_s.begin(),
              [](double v) {return pow(v, 2);});
    double stdev = sqrt(accumulate(diff_s.begin(), diff_s.end(), 0.0) / diff_s.size());
    vector<double> scaled_cov(aCovs.size());
    transform(diff.begin(), diff.end(), scaled_cov.begin(),
              [stdev](double v) {return abs(v / stdev);});
    return scaled_cov;
}

int LocalGenomicMap::balancerILP(const char * lpFn) {
    OsiClpSolverInterface * si = new OsiClpSolverInterface();
    
    // variable structure:
    // {segments (nSeg), junctions (nJunc), segment epsilon, junction epsilon}
    vector<Segment *> * segs = mGraph->getSegments();
    vector<Junction *> * juncs = mGraph->getJunctions();

    int numSegsJuncs = segs->size() + juncs->size();
    int numEpsilons = numSegsJuncs + 2;    // e(t(seg), c(seg)) 
                                           // e(t(junc), c(junc)) 
                                           // e(source, 2) 
                                           // e(sink, 2)
    // int numEpsilons = segs->size() * 3 + juncs->size() * 3 + 2;    
    // e(t(seg), c(seg)) 
    // e(t(seg), 0.4c(seg)) 
    // e(t(seg), 0.4c(seg)) 
    // e(t(junc), c(junc)) 
    // e(t(junc), 0.4c(junc)) 
    // e(t(junc), 0.8c(junc)) 
    // e(source, 2) 
    // e(sink, 2)
    int numVariables = numSegsJuncs + juncs->size() + numEpsilons;
    int numConstrains = segs->size() * 4 + 4 * juncs->size() + 2;  // "2" for source & sink      ------  t(seg)+c(seg), t(seg)-c(seg), in=t(seg), out=t(seg, t(junc)+c(junc), t(junc)-c(junc), e, source, sink
    // int numConstrains = segs->size() * 8 + juncs->size() * 6 + numEpsilons + 2;
    // t(seg)+c(seg), t(seg)-c(seg)
    // t(seg)+0.4c(seg), t(seg)-0.4c(seg)
    // t(seg)+0.8c(seg), t(seg)-0.8c(seg)
    // in=t(seg), out=t(seg)
    // t(junc)+c(junc), t(junc)-c(junc)
    // t(junc)+0.4c(junc), t(junc)-0.4c(junc)
    // t(junc)+0.8c(junc), t(junc)-0.8c(junc)
    // e
    // source, sink

    double * objective = new double[numVariables];
    double * variableLowerBound = new double[numVariables];
    double * variableUpperBound = new double[numVariables];
    double * constrainLowerBound = new double[numConstrains];
    double * constrainUpperBound = new double[numConstrains];
    cout << "Declare done" << endl;

    CoinPackedMatrix * matrix = new CoinPackedMatrix(false, 0, 0);

    double hap_cov = mGraph->getHaploidDepth();
    int maxCopy = 999999;
    double max_cov = -1;
    double max_cp = -1;
    double min_cov = 999999;
    vector<double> covs;
    // constrains for segments
    for (int i = 0; i < segs->size(); i++) {
        covs.push_back((*segs)[i]->getWeight()->getCoverage());
        if ((*segs)[i]->getWeight()->getCoverage() > max_cov) {
            // max_cov = (*segs)[i]->getWeight()->getCopyNum();
            max_cov = (*segs)[i]->getWeight()->getCoverage();
        }
        if ((*segs)[i]->getWeight()->getCopyNum() > max_cp) {
            // max_cov = (*segs)[i]->getWeight()->getCopyNum();
            max_cp = (*segs)[i]->getWeight()->getCopyNum();
        }
        if ((*segs)[i]->getWeight()->getCoverage() < min_cov) {
            // max_cov = (*segs)[i]->getWeight()->getCopyNum();
            min_cov = (*segs)[i]->getWeight()->getCoverage();
        }
        // t_i + e_i_1 >= c_i
        CoinPackedVector constrain1;
        // constrain1.insert(i, 1);  // t_i
        // constrain1.insert(numSegsJuncs + juncs->size() + i, 1);  // e_i_1
        // constrainLowerBound[4 * i] = (*segs)[i]->getWeight()->getCopyNum();
        // constrainUpperBound[4 * i] = si->getInfinity();
        constrain1.insert(i, hap_cov);  // t_i
        constrain1.insert(numSegsJuncs + juncs->size() + i, 1);  // e_i_1
        constrainLowerBound[4 * i] = (*segs)[i]->getWeight()->getCoverage();
        constrainUpperBound[4 * i] = si->getInfinity();
        matrix->appendRow(constrain1);

        // t_i - e_i_1 <= c_i
        CoinPackedVector constrain2;
        // constrain2.insert(i, 1);
        // constrain2.insert(numSegsJuncs + juncs->size() + i, -1);
        // constrainLowerBound[4 * i + 1] = -1 * si->getInfinity();
        // constrainUpperBound[4 * i + 1] = (*segs)[i]->getWeight()->getCopyNum();
        constrain2.insert(i, hap_cov);
        constrain2.insert(numSegsJuncs + juncs->size() + i, -1);
        constrainLowerBound[4 * i + 1] = -1 * si->getInfinity();
        constrainUpperBound[4 * i + 1] = (*segs)[i]->getWeight()->getCoverage();
        matrix->appendRow(constrain2);

        // // t_i + 0.4*c_i*e_i_2 >= c_i
        // CoinPackedVector constrain3;
        // constrain3.insert(i, 1);  // t_i
        // constrain3.insert(numSegsJuncs + 3 * i + 1, 0.4 * (*segs)[i]->getWeight()->getCopyNum());  // e_i_2
        // constrainLowerBound[8 * i + 2] = (*segs)[i]->getWeight()->getCopyNum();
        // constrainUpperBound[8 * i + 2] = si->getInfinity();
        // matrix->appendRow(constrain3);

        // // t_i - 0.4*c_i*e_i_2 <= c_i
        // CoinPackedVector constrain4;
        // constrain4.insert(i, 1);
        // constrain4.insert(numSegsJuncs + 3 * i + 1, -0.4 * (*segs)[i]->getWeight()->getCopyNum());
        // constrainLowerBound[8 * i + 3] = -1 * si->getInfinity();
        // constrainUpperBound[8 * i + 3] = (*segs)[i]->getWeight()->getCopyNum();
        // matrix->appendRow(constrain4);

        // // t_i + 0.8*c_i*e_i_3 >= c_i
        // CoinPackedVector constrain5;
        // constrain5.insert(i, 1);  // t_i
        // constrain5.insert(numSegsJuncs + 3 * i + 2, 0.8 * (*segs)[i]->getWeight()->getCopyNum());  // e_i_3
        // constrainLowerBound[8 * i + 4] = (*segs)[i]->getWeight()->getCopyNum();
        // constrainUpperBound[8 * i + 4] = si->getInfinity();
        // matrix->appendRow(constrain5);

        // // t_i - 0.8*c_i*e_i_2 <= c_i
        // CoinPackedVector constrain6;
        // constrain6.insert(i, 1);
        // constrain6.insert(numSegsJuncs + 3 * i + 2, -0.8 * (*segs)[i]->getWeight()->getCopyNum());
        // constrainLowerBound[8 * i + 5] = -1 * si->getInfinity();
        // constrainUpperBound[8 * i + 5] = (*segs)[i]->getWeight()->getCopyNum();
        // matrix->appendRow(constrain6);

        // constrains for in=out=copynum
        vector<int> visitedIdx;
        vector<int> visitedCount;
        vector<int>::iterator found;
        CoinPackedVector constrain7;
        // constrain7.insert(i, 1);
        // for (Edge * e : *((*segs)[i]->getPositiveVertex()->getEdgesAsTarget())) {
        //     int juncIndex = mGraph->getJunctionIndexByEdge(e);
        //     if (find(visitedIdx.begin(), visitedIdx.end(), juncIndex) 
        //         == visitedIdx.end()) {
        //         constrain7.insert(segs->size() + juncIndex, -1);
        //         visitedIdx.push_back(juncIndex);
        //     }
        // }
        constrain7.insert(i, 1);
        for (Edge * e : *((*segs)[i]->getPositiveVertex()->getEdgesAsTarget())) {
            // cout << e->getInfo() << endl;
            int juncIndex = mGraph->getJunctionIndexByEdge(e);
            found = find(visitedIdx.begin(), visitedIdx.end(), juncIndex);
            if (found == visitedIdx.end()) {
                // constrain7.insert(segs->size() + juncIndex, -1);
                visitedIdx.push_back(juncIndex);
                visitedCount.push_back(1);
            } else {
                if (e->getTarget() != e->getSource()) {
                    visitedCount[found - visitedIdx.begin()]++;
                }
            }
        }
        for (int j = 0; j < visitedIdx.size(); j++) {
            cout << visitedIdx[j] << ":" << visitedCount[j] << " ";
            constrain7.insert(segs->size() + visitedIdx[j], -visitedCount[j]);
        }
        cout << endl;
        constrainLowerBound[4 * i + 2] = 0;
        constrainUpperBound[4 * i + 2] = 0;
        matrix->appendRow(constrain7);
        visitedIdx.clear();
        visitedCount.clear();

        CoinPackedVector constrain8;
        // constrain8.insert(i, 1);
        // for (Edge * e : *((*segs)[i]->getPositiveVertex()->getEdgesAsSource())) {
        //     int juncIndex = mGraph->getJunctionIndexByEdge(e);
        //     if (find(visitedIdx.begin(), visitedIdx.end(), juncIndex) 
        //         == visitedIdx.end()) {
        //         constrain8.insert(segs->size() + juncIndex, -1);
        //         visitedIdx.push_back(juncIndex);
        //     }
        // }
        constrain8.insert(i, 1);
        for (Edge * e : *((*segs)[i]->getPositiveVertex()->getEdgesAsSource())) {
            int juncIndex = mGraph->getJunctionIndexByEdge(e);
            found = find(visitedIdx.begin(), visitedIdx.end(), juncIndex);
            if (found == visitedIdx.end()) {
                // constrain8.insert(segs->size() + juncIndex, -1);
                visitedIdx.push_back(juncIndex);
                visitedCount.push_back(1);
            } else {
                if (e->getTarget() != e->getSource()) {
                    visitedCount[found - visitedIdx.begin()]++;
                }
            }
        }
        for (int j = 0; j < visitedIdx.size(); j++) {
            cout << visitedIdx[j] << ":" << visitedCount[j] << " ";
            constrain8.insert(segs->size() + visitedIdx[j], -visitedCount[j]);
        }
        cout << endl;
        constrainLowerBound[4 * i + 3] = 0;
        constrainUpperBound[4 * i + 3] = 0;
        matrix->appendRow(constrain8);
        visitedIdx.clear();
        visitedCount.clear();
    }
    cout << "Segment done" << endl;

    // constrains for junctions
    for (int i = 0; i < juncs->size(); i++) {
        covs.push_back((*juncs)[i]->getWeight()->getCoverage());
        if ((*juncs)[i]->getWeight()->getCoverage() > max_cov) {
            // max_cov = (*juncs)[i]->getWeight()->getCopyNum();
            max_cov = (*juncs)[i]->getWeight()->getCoverage();
        }
        if ((*juncs)[i]->getWeight()->getCopyNum() > max_cp) {
            // max_cov = (*juncs)[i]->getWeight()->getCopyNum();
            max_cp = (*juncs)[i]->getWeight()->getCopyNum();
        }
        if ((*juncs)[i]->getWeight()->getCoverage() < min_cov) {
            // max_cov = (*juncs)[i]->getWeight()->getCopyNum();
            min_cov = (*juncs)[i]->getWeight()->getCoverage();
        }
        // t_i - c_ix_i + e_i_1 >= 0
        // double copy = (*juncs)[i]->getWeight()->getCopyNum() + 0.01;
        double cov = (*juncs)[i]->getWeight()->getCoverage() + 0.05;
        CoinPackedVector constrain1;
        // constrain1.insert(segs->size() + i, 1);
        // constrain1.insert(numSegsJuncs + i, -copy);
        // constrain1.insert(numSegsJuncs + juncs->size() + segs->size() + i, 1);
        // constrainLowerBound[4 * segs->size() + 3 * i] = 0;
        // constrainUpperBound[4 * segs->size() + 3 * i] = si->getInfinity();
        constrain1.insert(segs->size() + i, hap_cov);
        constrain1.insert(numSegsJuncs + i, -cov);
        constrain1.insert(numSegsJuncs + juncs->size() + segs->size() + i, 1);
        constrainLowerBound[4 * segs->size() + 4 * i] = 0;
        constrainUpperBound[4 * segs->size() + 4 * i] = si->getInfinity();
        matrix->appendRow(constrain1);

        // t_i - c_ix_i - e_i_1 <= 0
        CoinPackedVector constrain2;
        // constrain2.insert(segs->size() + i, 1);
        // constrain2.insert(numSegsJuncs + i, -copy);
        // constrain2.insert(numSegsJuncs + juncs->size() + segs->size() + i, -1);
        // constrainLowerBound[4 * segs->size() + 3 * i + 1] = -1 * si->getInfinity();
        // constrainUpperBound[4 * segs->size() + 3 * i + 1] = 0;
        constrain2.insert(segs->size() + i, hap_cov);
        constrain2.insert(numSegsJuncs + i, -cov);
        constrain2.insert(numSegsJuncs + juncs->size() + segs->size() + i, -1);
        constrainLowerBound[4 * segs->size() + 4 * i + 1] = -1 * si->getInfinity();
        constrainUpperBound[4 * segs->size() + 4 * i + 1] = 0;
        matrix->appendRow(constrain2);

        // t_i - x_iM <= 0
        CoinPackedVector constrain3;
        constrain3.insert(segs->size() + i, 1);
        constrain3.insert(numSegsJuncs + i, -maxCopy);
        constrainLowerBound[4 * segs->size() + 4 * i + 2] = -1 * si->getInfinity();
        constrainUpperBound[4 * segs->size() + 4 * i + 2] = 0;
        matrix->appendRow(constrain3);

        // t_i - x_i >= 0
        CoinPackedVector constrain4;
        constrain4.insert(segs->size() + i, 1);
        constrain4.insert(numSegsJuncs + i, -1);
        constrainLowerBound[4 * segs->size() + 4 * i + 3] = 0;
        constrainUpperBound[4 * segs->size() + 4 * i + 3] = si->getInfinity();
        matrix->appendRow(constrain4);

        // // t_i + 0.4*c_i*e_i_2 >= c_i
        // CoinPackedVector constrain3;
        // constrain3.insert(segs->size() + i, 1);
        // constrain3.insert(numSegsJuncs + segs->size() * 3 + 3 * i + 1, 0.4 * copy);
        // constrainLowerBound[8 * segs->size() + 6 * i + 2] = copy;
        // constrainUpperBound[8 * segs->size() + 6 * i + 2] = si->getInfinity();
        // matrix->appendRow(constrain3);

        // // t_i - 0.4*c_i*e_i_2 <= c_i
        // CoinPackedVector constrain4;
        // constrain4.insert(segs->size() + i, 1);
        // constrain4.insert(numSegsJuncs + segs->size() * 3 + 3 * i + 1, -0.4 * copy);
        // constrainLowerBound[8 * segs->size() + 6 * i + 3] = -1 * si->getInfinity();
        // constrainUpperBound[8 * segs->size() + 6 * i + 3] = copy;
        // matrix->appendRow(constrain4);

        // // t_i + 0.8*c_i*e_i_3 >= c_i
        // CoinPackedVector constrain5;
        // constrain5.insert(segs->size() + i, 1);
        // constrain5.insert(numSegsJuncs + segs->size() * 3 + 3 * i + 2, 0.8 * copy);
        // constrainLowerBound[8 * segs->size() + 6 * i + 4] = copy;
        // constrainUpperBound[8 * segs->size() + 6 * i + 4] = si->getInfinity();
        // matrix->appendRow(constrain5);

        // // t_i - 0.8*c_i*e_i_3 <= c_i
        // CoinPackedVector constrain6;
        // constrain6.insert(segs->size() + i, 1);
        // constrain6.insert(numSegsJuncs + segs->size() * 3 + 3 * i + 2, -0.8 * copy);
        // constrainLowerBound[8 * segs->size() + 6 * i + 5] = -1 * si->getInfinity();
        // constrainUpperBound[8 * segs->size() + 6 * i + 5] = copy;
        // matrix->appendRow(constrain6);
    }
    cout << "Junc done" << endl;

    // constrains for epsilons
    // for (int i = 0; i < numEpsilons; i++) {
    //     CoinPackedVector constrain;
    //     constrain.insert(numSegsJuncs + i, 1);
    //     constrainLowerBound[4 * segs->size() + 3 * juncs->size() + i] = 0;
    //     constrainUpperBound[4 * segs->size() + 3 * juncs->size() + i] = si->getInfinity();
    //     matrix->appendRow(constrain);
    // }

    // constrains for source and sink
    CoinPackedVector constrainSource;
    constrainSource.insert(mGraph->getSource()->getId() - 1, 1);
    constrainSource.insert(numVariables - 2, -1);
    // cout << mGraph->getExpectedPloidy() << endl;
    constrainLowerBound[numConstrains - 2] = mGraph->getExpectedPloidy();
    constrainUpperBound[numConstrains - 2] = mGraph->getExpectedPloidy();
    matrix->appendRow(constrainSource);

    CoinPackedVector constrainSink;
    constrainSink.insert(mGraph->getSink()->getId() - 1, 1);
    constrainSink.insert(numVariables - 1, -1);
    constrainLowerBound[numConstrains - 1] = mGraph->getExpectedPloidy();
    constrainUpperBound[numConstrains - 1] = mGraph->getExpectedPloidy();
    matrix->appendRow(constrainSink);
    cout << "Source sink done" << endl;

    // objective function: coefficient is 0 except for epsilon variables
    // TODO
    vector<double> coefs = scaleILPCoef(covs);
    double max_coef = *max_element(coefs.begin(), coefs.end());
    double min_coef = *min_element(coefs.begin(), coefs.end());
    min_coef = (min_coef > 0) ? min_coef : 0.1;
    max_cov += 1000;
    for (int i = 0; i < numVariables; i++) {
        if (i >= numSegsJuncs) {
            if (i < numVariables - 2) {
                // double cred;
                if (i < numSegsJuncs + juncs->size()) {
                    if ((*juncs)[i - numSegsJuncs]->isInferred()) {
                        objective[i] = max_coef;
                        // objective[i] = max_cov;
                        // objective[i] = max_cp + 100;
                        // objective[i] = exp10((max_cov - min_cov) / (max_cov - min_cov));
                        // objective[i] = log10(max_cov + 1);
                    } else {
                        objective[i] = 0;
                    }
                } else if (i < numSegsJuncs + juncs->size() + segs->size()) {
                    // cred = (*segs)[i - numSegsJuncs]->getCredibility();
                    // objective[i] = (*segs)[i - numSegsJuncs - juncs->size()]->getWeight()->getCoverage();
                    // objective[i] = (*segs)[i - numSegsJuncs - juncs->size()]->getWeight()->getCopyNum();
                    // objective[i] = exp10(((*segs)[i - numSegsJuncs - juncs->size()]->getWeight()->getCoverage() - min_cov) / (max_cov - min_cov));
                    // objective[i] = log10((*segs)[i - numSegsJuncs - juncs->size()]->getWeight()->getCoverage() + 1);
                    // objective[i] = 1;
                    objective[i] = coefs[i - numSegsJuncs - juncs->size()];
                } else {
                    // cred = (*juncs)[i - numSegsJuncs - segs->size()]->getCredibility();
                    // objective[i] = exp10(((*juncs)[i - numSegsJuncs - juncs->size() - segs->size()]->getWeight()->getCoverage() - min_cov) / (max_cov - min_cov));
                    // objective[i] = log10((*juncs)[i - numSegsJuncs - juncs->size() - segs->size()]->getWeight()->getCoverage() + 1);
                    // objective[i] = (*juncs)[i - numSegsJuncs - juncs->size() - segs->size()]->getWeight()->getCoverage();
                    // objective[i] = (*juncs)[i - numSegsJuncs - juncs->size() - segs->size()]->getWeight()->getCopyNum();
                    // ]->getWeight()->getCoverage();
                    objective[i] = min_coef;
                    // objective[i] = coefs[i - numSegsJuncs - juncs->size() - segs->size()];
                }
                // objective[i] = 1 * cred;
                // objective[i] = 1;
            } else {
                objective[i] = max_coef;
                // objective[i] = max_cov;
                // objective[i] = max_cp + 100;
                // objective[i] = exp10((max_cov - min_cov) / (max_cov - min_cov));
                // objective[i] = log10(max_cov + 1);
                // objective[i] = 1;
            }
        } else {
            objective[i] = 0;
        }
    }
    cout << "Objective done" << endl;
    // for (int i = 0; i < numSegsJuncs; i++) {
    //     objective[i] = 0;
    // }
    // for (int i = 0; i < segs->size(); i++) {
    //     objective[numSegsJuncs + 3 * i] = 1;
    //     objective[numSegsJuncs + 3 * i + 1] = 100;
    //     objective[numSegsJuncs + 3 * i + 2] = 10000;
    // }
    // for (int i = 0; i < juncs->size(); i++) {
    //     objective[numSegsJuncs + 3 * segs->size() + 3 * i] = 1;
    //     objective[numSegsJuncs + 3 * segs->size() + 3 * i + 1] = 100;
    //     objective[numSegsJuncs + 3 * segs->size() + 3 * i + 2] = 10000;
    // }
    
    // lower & upper bounds for segments and corresponding epsilons variables
    for (int i = 0; i < segs->size(); i++) {
        if ((*segs)[i]->hasLowerBoundLimit()) {  // inferred junctions have lower bound 1
            variableLowerBound[i] = 1;
        } else {
            variableLowerBound[i] = 0;
        }
        variableUpperBound[i] = si->getInfinity();

        // for epsilons
        variableLowerBound[numSegsJuncs + juncs->size() + i] = 0;
        variableUpperBound[numSegsJuncs + juncs->size() + i] = si->getInfinity();
        // variableLowerBound[numSegsJuncs + 3 * i] = 0;
        // variableUpperBound[numSegsJuncs + 3 * i] = si->getInfinity();
        // variableLowerBound[numSegsJuncs + 3 * i + 1] = 0;
        // variableUpperBound[numSegsJuncs + 3 * i + 1] = si->getInfinity();
        // variableLowerBound[numSegsJuncs + 3 * i + 2] = 0;
        // variableUpperBound[numSegsJuncs + 3 * i + 2] = si->getInfinity();
    }
    cout << "LU seg done" << endl;

    // lower & upper bounds for junctions and corresponding epsilons variables
    for (int i = 0; i < juncs->size(); i++) {
        // if ((*juncs)[i]->hasLowerBoundLimit()) {
        //     variableLowerBound[segs->size() + i] = 1;
        // } else {
        //     variableLowerBound[segs->size() + i] = 0;
        // }
        variableLowerBound[segs->size() + i] = 0;
        variableUpperBound[segs->size() + i] = si->getInfinity();

        // for x
        // if (!(*juncs)[i]->hasLowerBoundLimit()) {
        if ((*juncs)[i]->isInferred()) {
            variableLowerBound[numSegsJuncs + i] = 0;
            cout << "infer info - 0: " << i << " " << segs->size() + i << " " << numSegsJuncs + i << " " << (*juncs)[i]->getInfo()[0] << " " << endl;
        } else {
            variableLowerBound[numSegsJuncs + i] = 1;
            cout << "infer info - 1: " << i << " " << segs->size() + i << " " << numSegsJuncs + i << " " << (*juncs)[i]->getInfo()[0] << " " << endl;
        }
        variableUpperBound[numSegsJuncs + i] = 1;

        // for epsilons
        variableLowerBound[numSegsJuncs + juncs->size() + segs->size() + i] = 0;
        variableUpperBound[numSegsJuncs + juncs->size() + segs->size() + i] = si->getInfinity();
        // variableLowerBound[numSegsJuncs + 3 * segs->size() + 3 * i] = 0;
        // variableUpperBound[numSegsJuncs + 3 * segs->size() + 3 * i] = si->getInfinity();
        // variableLowerBound[numSegsJuncs + 3 * segs->size() + 3 * i + 1] = 0;
        // variableUpperBound[numSegsJuncs + 3 * segs->size() + 3 * i + 1] = si->getInfinity();
        // variableLowerBound[numSegsJuncs + 3 * segs->size() + 3 * i + 2] = 0;
        // variableUpperBound[numSegsJuncs + 3 * segs->size() + 3 * i + 2] = si->getInfinity();
    }
    cout << "LU junc done" << endl;
    variableLowerBound[numVariables - 2] = 0;
    variableUpperBound[numVariables - 2] = si->getInfinity();
   //  variableUpperBound[numVariables - 2] = 0;
    variableLowerBound[numVariables - 1] = 0;
    variableUpperBound[numVariables - 1] = si->getInfinity();
    // variableUpperBound[numVariables - 1] = 0;

    si->loadProblem(*matrix, variableLowerBound, variableUpperBound, objective, constrainLowerBound, constrainUpperBound);
    // set integer for segment and junction variables
    for (int i = 0; i < numSegsJuncs + juncs->size(); i++) {
        si->setInteger(i);
    }
    // if (access(lpFn, F_OK) == -1)    {
    //     cout << "Cannot open file " << lpFn << ": no such file or directory" << endl;
    //     exit(1);
    // }
    si->writeMps(lpFn);
    si->writeLp(lpFn);
    // si->initialSolve();

    CbcModel model(*si);

    // const double * lower = model.getColLower();
    // const double * upper = model.getColUpper();

    // model.setLogLevel(1);
    // model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    // if (model.isProvenInfeasible()) {
    //     return -1;
    //     // throw ILPBalancerInfeasibleEuxception();
    // } else {
    //     cout << "Feasible" << endl;
    // }
    // model.branchAndBound();
    // const double * solution = model.solver()->getColSolution();
    // for (int i = 0; i < numVariables; i++) {
    //     if (i < segs->size()) {
    //         if (model.solver()->isInteger(i)) {
    //             cout << "S_" << (*segs)[i]->getId() << "(I) = " << solution[i] << endl;
    //         } else {
    //             cout << "S_" << (*segs)[i]->getId() << "(C) = " << solution[i] << endl;
    //         }
    //     }
    //     if (i >= segs->size() && i < numSegsJuncs) {
    //         vector<string> juncInfo = (*juncs)[i - segs->size()]->getInfo();
    //         if (model.solver()->isInteger(i)) {
    //             cout << "J_" << i - segs->size() << "(I) = " << solution[i] << ", " << juncInfo[0] << " " << juncInfo[1] << endl;
    //         } else {
    //             cout << "J_" << i - segs->size() << "(C) = " << solution[i] << ", " << juncInfo[0] << " " << juncInfo[1] << endl;
    //         }
    //     }
    //     if (i >= numSegsJuncs) {
    //         if (model.solver()->isInteger(i)) {
    //             cout << "e_" << i - numSegsJuncs << "(I) = " << solution[i] << endl;
    //         } else {
    //             cout << "e_" << i - numSegsJuncs << "(C) = " << solution[i] << endl;
    //         }
    //     }
    // }
    // cout << endl;

    // // write result to graph
    // for (int i = 0; i < segs->size(); i++) {
    //     (*segs)[i]->getWeight()->setCopyNum(solution[i]);
    // }
    // for (int i = 0; i < juncs->size(); i++) {
    //     (*juncs)[i]->getWeight()->setCopyNum(solution[segs->size() + i]);
    // }
    return 0;
}

void LocalGenomicMap::addAllJuncsFromDB(JunctionDB *aJuncDB) {
    // bool added = false;
    // while (!this->isReachable()) {
    for (Record *rec: *(aJuncDB->getRecords())) {
        // added = false;
        for (entry_t *ent: *(rec->getForwardEntries())) {
            // if (ent->support < 2) continue;
            Vertex *currentVertex;
            cout << "Rec forward: " << rec->getChrom() << " " << rec->getPos()
                 << rec->getStrand() << endl;
            try {
                if (rec->getStrand() == '+') {
                    currentVertex = mGraph->getSegmentByChromEnd(
                        rec->getChrom(), rec->getPos()
                    )->getPositiveVertex();
                } else {
                    currentVertex = mGraph->getSegmentByChromStart(
                        rec->getChrom(), rec->getPos()
                    )->getNegativeVertex();
                }
            } catch (SegmentDoesNotExistException &e) {
                cout << e.what() << endl;
                // char t;
                // cin >> t;
                continue;
            }
            Vertex *nextVertex;
            cout << "Ent forward: " << ent->chrom << " " << ent->pos << endl;
            try {
                if (ent->strand == '+') {
                    nextVertex = mGraph->getSegmentByChromStart(
                        ent->chrom, ent->pos
                    )->getPositiveVertex();
                } else {
                    nextVertex = mGraph->getSegmentByChromEnd(
                        ent->chrom, ent->pos
                    )->getNegativeVertex();
                }
            } catch (SegmentDoesNotExistException &e) {
                cout << e.what() << endl;
                // char t;
                // cin >> t;
                continue;
            }
            try {
                Junction *inferredJunc = mGraph->addJunction(
                    currentVertex, nextVertex,
                    this->inferCoverage(currentVertex, nextVertex),
                    this->inferCredibility(currentVertex, nextVertex),
                    -1, true, false, false
                );
                if (inferredJunc == NULL) {
                    continue;
                }
                cout << "add junction with sample juncdb: "
                     << inferredJunc->getInfo()[0] << ", "
                     << inferredJunc->getInfo()[1] << endl;
                // added = true;
                // break;
            } catch (DuplicateJunctionException &e) {
                cout << e.what() << endl;
            }
        }
        for (entry_t *ent: *(rec->getBackwardEntries())) {
            // if (ent->support < 2) continue;
            Vertex *currentVertex;
            cout << "Rec backward: " << rec->getChrom() << " " << rec->getPos()
                 << rec->getStrand() << endl;
            try {
                if (rec->getStrand() == '+') {
                    currentVertex = mGraph->getSegmentByChromStart(
                        rec->getChrom(), rec->getPos()
                    )->getPositiveVertex();
                } else {
                    currentVertex = mGraph->getSegmentByChromEnd(
                        rec->getChrom(), rec->getPos()
                    )->getNegativeVertex();
                }
            } catch (SegmentDoesNotExistException &e) {
                cout << e.what() << endl;
                // char t;
                // cin >> t;
                continue;
            }
            Vertex *prevVertex;
            cout << "Ent backward: " << ent->chrom << " " << ent->pos << endl;
            try {
                if (ent->strand == '+') {
                    prevVertex = mGraph->getSegmentByChromEnd(
                        ent->chrom, ent->pos
                    )->getPositiveVertex();
                } else {
                    prevVertex = mGraph->getSegmentByChromStart(
                        ent->chrom, ent->pos
                    )->getNegativeVertex();
                }
            } catch (SegmentDoesNotExistException &e) {
                cout << e.what() << endl;
                // char t;
                // cin >> t;
                continue;
            }
            try {
                Junction *inferredJunc = mGraph->addJunction(
                    prevVertex, currentVertex,
                    this->inferCoverage(prevVertex, currentVertex),
                    this->inferCredibility(prevVertex, currentVertex),
                    -1, true, false, false
                );
                if (inferredJunc == NULL) {
                    continue;
                }
                cout << "add junction with sample juncdb: "
                     << inferredJunc->getInfo()[0] << ", "
                     << inferredJunc->getInfo()[1] << endl;
                // added = true;
                // break;
            } catch (DuplicateJunctionException &e) {
                cout << e.what() << endl;
            }
        }
    }
    // if (!added) {
    //     break;
    // }
    // }
}

int LocalGenomicMap::strandCross(Vertex * aVertex) {
    int flag = 0;
    Vertex * startVertex = aVertex;
    VertexPath vertexStack;
    vertexStack.push_back(startVertex);
    Edge * prevEdge = this->selectPrevEdge(startVertex);
    while (true) {
        if (prevEdge == NULL) {
            vertexStack.pop_back();
            if (vertexStack.size() == 0) {
                break;
            } else {
                startVertex = vertexStack.back();
                prevEdge = this->selectPrevEdge(startVertex);
            }
        } else {
            Vertex * prevVertex = prevEdge->getSource();
            // cout << nextEdge->getInfo() << ", forward" << endl;

            if (prevVertex->getDir() != aVertex->getDir()) {
                flag |= 1 << 1;
                break;
            }

            prevEdge->setVisited();
            // nextEdge->getComplementEdge()->setVisited();
            vertexStack.push_back(prevVertex);
            startVertex = prevVertex;
            prevEdge = this->selectPrevEdge(startVertex);
        }
    }
    mGraph->resetJunctionVisitFlag();

    startVertex = aVertex;
    vertexStack.push_back(startVertex);
    Edge * nextEdge = this->selectNextEdge(startVertex);
    while (true) {
        if (nextEdge == NULL) {
            vertexStack.pop_back();
            if (vertexStack.size() == 0) {
                break;
            } else {
                startVertex = vertexStack.back();
                nextEdge = this->selectNextEdge(startVertex);
            }
        } else {
            Vertex * nextVertex = nextEdge->getTarget();
            // cout << nextEdge->getInfo() << ", forward" << endl;

            if (nextVertex->getDir() != aVertex->getDir()) {
                flag |= 1;
                break;
            }

            nextEdge->setVisited();
            // nextEdge->getComplementEdge()->setVisited();
            vertexStack.push_back(nextVertex);
            startVertex = nextVertex;
            nextEdge = this->selectNextEdge(startVertex);
        }
    }
    mGraph->resetJunctionVisitFlag();

    return flag;
}

bool LocalGenomicMap::doesPathExists(Vertex * aStartVertex, Vertex * aEndVertex) {
    bool isReach = false;
    Vertex * startVertex = aStartVertex;
    VertexPath vertexStack;
    vertexStack.push_back(startVertex);
    Edge * nextEdge = this->selectNextEdge(startVertex);
    while (true) {
        if (nextEdge == NULL) {
            vertexStack.pop_back();
            if (vertexStack.size() == 0) {
                break;
            } else {
                startVertex = vertexStack.back();
                nextEdge = this->selectNextEdge(startVertex);
            }
        } else {
            Vertex * nextVertex = nextEdge->getTarget();
            // cout << nextEdge->getInfo() << ", forward" << endl;

            // if (nextVertex == mGraph->getSink()->getNegativeVertex()) {
            //     if (startVertex != mGraph->getSource()->getNegativeVertex()) {
            //         throw ForwardReachSinkNegativeException(aStartVertex);
            //     }
            // }
            // if (nextVertex == mGraph->getSource()->getPositiveVertex()) {
            //     if (startVertex != mGraph->getSink()->getPositiveVertex()) {
            //         throw ForwardReachSourcePositiveException(aStartVertex);
            //     }
            // }
            if (nextVertex == aEndVertex) {
                isReach = true;
                break;
            }
            nextEdge->setVisited();
            // nextEdge->getComplementEdge()->setVisited();
            vertexStack.push_back(nextVertex);
            startVertex = nextVertex;
            nextEdge = this->selectNextEdge(startVertex);
        }
    }
    mGraph->resetJunctionVisitFlag();
    return isReach;
}

double LocalGenomicMap::inferCoverage(Vertex * aSource, Vertex * aTarget) {
    double inCoverage = aTarget->getWeight()->getCoverage() - aTarget->getInCoverage();
    double outCoverage = aSource->getWeight()->getCoverage() - aSource->getOutCoverage();
    return max(1.0, (inCoverage + outCoverage) / 2.0);
}

double LocalGenomicMap::weightedCredibility(Vertex * aVertex, bool aIsSource) {
    // the coefficient is modified by a constant setting
    if (aIsSource) {
        return 1.0 * aVertex->getCredibility() * max(0.0, aVertex->getWeight()->getCoverage() - aVertex->getOutCoverage()) / mGraph->getAvgCoverage();
    } else {
        return 1.0 * aVertex->getCredibility() * max(0.0, aVertex->getWeight()->getCoverage() - aVertex->getInCoverage()) / mGraph->getAvgCoverage();
    }
}

double LocalGenomicMap::inferCredibility(Vertex * aSource, Vertex * aTarget) {
    return (this->weightedCredibility(aSource, true) + this->weightedCredibility(aTarget, false)) / 2;
}

void LocalGenomicMap::connectSourceSink() {
    mGraph->addJunction(mGraph->getSink()->getPositiveVertex(), mGraph->getSource()->getPositiveVertex(), (mGraph->getSource()->getWeight()->getCoverage() + mGraph->getSink()->getWeight()->getCoverage()) / 2, 1.0, -1, true, false, true);
}

void LocalGenomicMap::countReachability(int * nReachability, VertexPath & notReachableVertices, Vertex * targetVertex) {
    cout << "count reach" << endl;
    // mGraph->print();
    vector< pair<Vertex *, int> > vertex_reach_count;
    for (VertexPath::const_iterator i = notReachableVertices.begin(); i != notReachableVertices.end(); i++) {
        // cout << "Counting Reachability for " << (*i)->getInfo() << endl;
        int count = 0;
        for (VertexPath::const_iterator j = notReachableVertices.begin(); j != notReachableVertices.end(); j++) {
            if (j == i) continue;
            if (targetVertex == mGraph->getSource()->getPositiveVertex() || targetVertex == mGraph->getSink()->getNegativeVertex()) {
                if (mGraph->BFS(*i, *j) >= 0) {
                // if (this->doesPathExists(*i, *j)) {
                    count++;
                    // cout << (*i)->getInfo() << "=>" << (*j)->getInfo() << " connectable" << endl;
                } else {
                    // cout << (*i)->getInfo() << "=>" << (*j)->getInfo() << " unconnectable" << endl;
                }
            } else if (targetVertex == mGraph->getSource()->getNegativeVertex() || targetVertex == mGraph->getSink()->getPositiveVertex()) {
                if (mGraph->BFS(*j, *i) >= 0) {
                // if (this->doesPathExists(*j, *i)) {
                    count++;
                    // cout << (*i)->getInfo() << "=>" << (*j)->getInfo() << " connectable" << endl;
                } else {
                    // cout << (*i)->getInfo() << "=>" << (*j)->getInfo() << " unconnectable" << endl;
                }
            }
        }
        vertex_reach_count.push_back(pair<Vertex *, int>(*i, count));
        nReachability[i - notReachableVertices.begin()] = count;
    }
    sort(vertex_reach_count.begin(), vertex_reach_count.end(), [](pair<Vertex *, int> e1, pair<Vertex *, int> e2) { 
        if (e1.second < e2.second) {
            return true;
        } else if (e1.second == e2.second) {
            if (e1.first->getWeight()->getCoverage() < e2.first->getWeight()->getCoverage()) {
                return true;
            }
            // if (e1.first->getId() > e2.first->getId()) {
            //     return true;
            // }
        }
        return false;
        // return e1.second < e2.second; 
        });
    notReachableVertices.clear();
    for (vector< pair<Vertex *, int> >::reverse_iterator i = vertex_reach_count.rbegin(); i != vertex_reach_count.rend(); i++) {
        notReachableVertices.push_back(i->first);
    }
}

void LocalGenomicMap::findWithReachability(int * nReachability, VertexPath & notReachableVertices, int reachability) {
    VertexPath backup = notReachableVertices;
    notReachableVertices.clear();
    for (int i = 0; i < backup.size(); i++) {
        if (nReachability[i] == reachability) {
            notReachableVertices.push_back(backup[i]);
        }
    }
}

void LocalGenomicMap::reconnectNegativeToPositive(JunctionDB * aJuncDB, bool verbose) {
    Vertex * v;
    Record * rec;
    Junction * inferredJunc;
    VertexPath asSourceVertices, asTargetVertices;
    for (Segment * seg: *(mGraph->getSegments())) {
        if (seg->getId() > mGraph->getSink()->getId()) {
            cout << "Continue: " << seg->getId() << endl;
            continue;
        }
        v = seg->getNegativeVertex();
        bool asTargetFromPositive = false, asSourceToPositive = false;
        for (Edge * e: *(v->getEdgesAsSource())) {
            if (e->getTarget()->getDir() == '+') {
                asSourceToPositive = true;
                break;
            }
        }
        for (Edge * e: *(v->getEdgesAsTarget())) {
            if (e->getSource()->getDir() == '+') {
                asTargetFromPositive = true;
                break;
            }
        }
        if (asTargetFromPositive != asSourceToPositive) {
            if (asTargetFromPositive) {
                asTargetVertices.push_back(v);
            } else {
                asSourceVertices.push_back(v);
            }
        }
    }
    // cout << asTargetVertices.size() << " " << asSourceVertices.size() << endl;
    // int q;
    // cin >> q;
    while (asTargetVertices.size() > 0 || asSourceVertices.size() > 0) {
    for (Vertex * v: asSourceVertices) {
        cout << "Current vertex: " << v->getInfo() << " " << v->getSegment()->getChrom() << " " << v->getStart() << v->getDir() << endl;
        cout << "Finding previous vertex" << endl;
        // v = seg->getNegativeVertex();
        Segment * seg = v->getSegment();
        // for (Edge * e: *(v->getEdgesAsSource())) {
            // V- => V+, add junction s.t. V+ => V-
            Vertex * prevVertex;

            // Add junction to prev and next segment
            // cout << "Seg id: " << seg->getId() << endl;
            // prevVertex = mGraph->getSegmentById(seg->getId() - 1)->getPositiveVertex();
            // try {
            //     if (verbose)
            //         cout << "prev: " << prevVertex->getInfo() << endl;
            //     inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex, v), this->inferCredibility(prevVertex, v), true, false);
            //     if (verbose)
            //         cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //     cout << "Done" << endl;
            // } catch (DuplicateJunctionException &e) {
            //     cout << e.what() << endl;
            // }
            // prevVertex = mGraph->getSegmentById(seg->getId() + 1)->getPositiveVertex();
            // try {
            //     if (verbose)
            //         cout << "prev: " << prevVertex->getInfo() << endl;
            //     inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex, v), this->inferCredibility(prevVertex, v), true, false);
            //     if (verbose)
            //         cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            // } catch (DuplicateJunctionException &e) {
            //     cout << e.what() << endl;
            // }

            bool added = false;
            rec = aJuncDB->findRecord(v->getSegment()->getChrom(), v->getStart(), v->getDir());
            if (rec != NULL && rec->getBackwardEntries()->size() > 0) {
                for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin(); riter != rec->getBackwardEntries()->rend(); riter++) {
                    cout << (*riter)->chrom << " " << (*riter)->pos << (*riter)->strand << endl;
                    if ((*riter)->strand == '+') {
                        try {
                            prevVertex = mGraph->getSegmentByChromEnd((*riter)->chrom, (*riter)->pos)->getPositiveVertex();
                            if (verbose)
                                cout << "prev: " << prevVertex->getInfo() << " " << prevVertex->getSegment()->getChrom() << " " << v->getEnd() << v->getDir() << endl;
                            inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex,v ), this->inferCredibility(prevVertex, v), -1, true, false, false);
                            if (inferredJunc == NULL) continue;
                            if (verbose)
                                cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                            added = true;
                            break;
                        } catch (DuplicateJunctionException& e) {
                            cout << e.what() << endl;
                        } catch (SegmentDoesNotExistException& e) {
                            cout << e.what() << endl;
                        }
                    }
                }
            }
            if (!added) {
                cout << " As source N2P not added" << endl;
                prevVertex = v;
                while (true) {
                    if (prevVertex->getId() >= mGraph->getSink()->getId() ||prevVertex->getId() == 1) {
                        // prevVertex = mGraph->getSink()->getPositiveVertex();
                        // prevVertex = mGraph->getSegmentById(mGraph->getSink()->getId() - 1)->getPositiveVertex();
                        prevVertex = mGraph->getSegmentById(mGraph->getSink()->getId() - 1)->getPositiveVertex();
                        // prevVertex = mGraph->getSink()->getPositiveVertex();
                    } else {
                        prevVertex = mGraph->getSegmentById(prevVertex->getId() - 1)->getPositiveVertex();
                    }
                try {
                    if (verbose)
                        cout << "prev: " << prevVertex->getInfo() << " " << prevVertex->getSegment()->getChrom() << " " << v->getEnd() << v->getDir() << endl;
                    inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex,v ), this->inferCredibility(prevVertex, v), -1, true, false, false);
                    if (inferredJunc == NULL) {
                        continue;
                    };
                    if (verbose)
                        cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                    break;
                } catch (DuplicateJunctionException& e) {
                    cout << e.what() << endl;
                }
                }
            }
            // if (rec == NULL || rec->getBackwardEntries()->size() == 0) {
            //     prevVertex = mGraph->getSegmentById(v->getId() - 1)->getPositiveVertex();
            //     try {
            //         if (verbose)
            //             cout << "prev: " << prevVertex->getInfo() << " " << prevVertex->getSegment()->getChrom() << " " << v->getEnd() << v->getDir() << endl;
            //         inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex,v ), this->inferCredibility(prevVertex, v), true, false);
            //         if (verbose)
            //             cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //         break;
            //     } catch (DuplicateJunctionException& e) {
            //         cout << e.what() << endl;
            //     }
            // } else {
            //     for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin(); riter != rec->getBackwardEntries()->rend(); riter++) {
            //         cout << (*riter)->chrom << " " << (*riter)->pos << (*riter)->strand << endl;
            //         if ((*riter)->strand == '+') {
            //             prevVertex = mGraph->getSegmentByChromEnd((*riter)->chrom, (*riter)->pos)->getPositiveVertex();
            //             try {
            //                 if (verbose)
            //                     cout << "prev: " << prevVertex->getInfo() << " " << prevVertex->getSegment()->getChrom() << " " << v->getEnd() << v->getDir() << endl;
            //                 inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex,v ), this->inferCredibility(prevVertex, v), true, false);
            //                 if (verbose)
            //                     cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //                 break;
            //             } catch (DuplicateJunctionException& e) {
            //                 cout << e.what() << endl;
            //             }
            //         }
            //     }
            // }
        // }
    }
    for (Vertex * v: asTargetVertices) {
        cout << "Current vertex: " << v->getInfo() << " " << v->getSegment()->getChrom() << " " << v->getEnd() << v->getDir() << endl;
        cout << "Finding next vertex" << endl;
        Segment * seg = v->getSegment();
        // for (Edge * e: *(v->getEdgesAsTarget())) {
            // V- => V+, add junction s.t. V+ => V-
            Vertex * nextVertex;

            // Add junction to prev and next segment
            // nextVertex = mGraph->getSegmentById(seg->getId() + 1) -> getPositiveVertex();
            // try {
            //     if (verbose)
            //         cout << "prev: " << nextVertex->getInfo() << endl;
            //     inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex), this->inferCredibility(v, nextVertex), true, false);
            //     if (verbose)
            //         cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            // } catch (DuplicateJunctionException &e) {
            //     cout << e.what() << endl;
            // }
            // nextVertex = mGraph->getSegmentById(seg->getId() - 1) -> getPositiveVertex();
            // try {
            //     if (verbose)
            //         cout << "prev: " << nextVertex->getInfo() << endl;
            //     inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex), this->inferCredibility(v, nextVertex), true, false);
            //     if (verbose)
            //         cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            // } catch (DuplicateJunctionException &e) {
            //     cout << e.what() << endl;
            // }

            bool added = false;
            rec = aJuncDB->findRecord(seg->getChrom(), seg->getEnd(), '-');
            if (rec != NULL && rec->getForwardEntries()->size() > 0) {
                for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin(); riter != rec->getForwardEntries()->rend(); riter++) {
                    if ((*riter)->strand == '+') {
                        try {
                            nextVertex = mGraph->getSegmentByChromStart((*riter)->chrom, (*riter)->pos)->getPositiveVertex();
                            if (verbose)
                                cout << "next: " << nextVertex->getInfo() << " " << nextVertex->getSegment()->getChrom() << " " << v->getStart() << v->getDir() << endl;
                                // cout << "next: " << nextVertex->getInfo() << endl;
                            inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex), this->inferCredibility(v, nextVertex), -1, true, false, false);
                            if (inferredJunc == NULL) continue;
                            if (verbose)
                                cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                            added = true;
                            break;
                        } catch (DuplicateJunctionException& e) {
                            cout << e.what() << endl;
                        } catch (SegmentDoesNotExistException& e) {
                            cout << e.what() << endl;
                        }
                    }
                }
            }
            if (!added) {
                cout << " As target N2P not added" << endl;
                nextVertex = v;
                while (true) {
                    if (nextVertex->getId() >= mGraph->getSink()->getId() || nextVertex->getId() == mGraph->getSink()->getId()) {
                        // nextVertex = mGraph->getSource()->getPositiveVertex();
                        nextVertex = mGraph->getSegmentById(mGraph->getSource()->getId() + 1)->getPositiveVertex();
                    } else {
                        // nextVertex = mGraph->getSegmentById(mGraph->getSource()->getId() + 1)->getPositiveVertex();
                        nextVertex = mGraph->getSegmentById(nextVertex->getId() + 1)->getPositiveVertex();
                    }
                    // if (nextVertex->getId() > mGraph->getSink()->getId()) {
                    //     continue;
                    // } // TODO
                    // if (nextVertex->getId() == mGraph->getSink()->getId()) {
                    //     nextVertex = mGraph->getSource()->getPositiveVertex();
                    // }
                // nextVertex = mGraph->getSegmentById(nextVertex->getId() + 1)->getPositiveVertex();
                try {
                    if (verbose)
                        cout << "next: " << nextVertex->getInfo() << " " << nextVertex->getSegment()->getChrom() << " " << v->getStart() << v->getDir() << endl;
                    inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex), this->inferCredibility(v, nextVertex), -1, true, false, false);
                    if (inferredJunc == NULL) {
                        continue;
                    };
                    if (verbose)
                        cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                    break;
                } catch (DuplicateJunctionException& e) {
                    cout << e.what() << endl;
                }
                }
            }
            // if (rec == NULL || rec->getForwardEntries()->size() == 0) {
            //     nextVertex = mGraph->getSegmentById(v->getId() + 1)->getPositiveVertex();
            //     try {
            //         if (verbose)
            //             cout << "next: " << nextVertex->getInfo() << " " << nextVertex->getSegment()->getChrom() << " " << v->getStart() << v->getDir() << endl;
            //         inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex), this->inferCredibility(v, nextVertex), true, false);
            //         if (verbose)
            //             cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //         break;
            //     } catch (DuplicateJunctionException& e) {
            //         cout << e.what() << endl;
            //     }
            // } else {
            //     for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin(); riter != rec->getForwardEntries()->rend(); riter++) {
            //         if ((*riter)->strand == '+') {
            //             nextVertex = mGraph->getSegmentByChromStart((*riter)->chrom, (*riter)->pos)->getPositiveVertex();
            //             try {
            //                 if (verbose)
            //                     cout << "next: " << nextVertex->getInfo() << " " << nextVertex->getSegment()->getChrom() << " " << v->getStart() << v->getDir() << endl;
            //                     // cout << "next: " << nextVertex->getInfo() << endl;
            //                 inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex), this->inferCredibility(v, nextVertex), true, false);
            //                 if (verbose)
            //                     cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //                 break;
            //             } catch (DuplicateJunctionException& e) {
            //                 cout << e.what() << endl;
            //             }
            //         }
            //     }
            // }
        //}
    }
    asTargetVertices.clear();
    asSourceVertices.clear();
    for (Segment * seg: *(mGraph->getSegments())) {
        v = seg->getNegativeVertex();
        bool asTargetFromPositive = false, asSourceToPositive = false;
        for (Edge * e: *(v->getEdgesAsSource())) {
            if (e->getTarget()->getDir() == '+') {
                asSourceToPositive = true;
                break;
            }
        }
        for (Edge * e: *(v->getEdgesAsTarget())) {
            if (e->getSource()->getDir() == '+') {
                asTargetFromPositive = true;
                break;
            }
        }
        if (asTargetFromPositive != asSourceToPositive) {
            if (asTargetFromPositive) {
                asTargetVertices.push_back(v);
            } else {
                asSourceVertices.push_back(v);
            }
        }
    }
    // cout << asTargetVertices.size() << " " << asSourceVertices.size() << endl;
    // int t;
    // cin >> t;
    }
}

Vertex * LocalGenomicMap::getMostReachable(VertexPath & notReachableVertices, Vertex * targetVertex) {
    Vertex * mostReachable = NULL;
    if (targetVertex == mGraph->getSource()->getPositiveVertex()) {
        // use the least id
        int mostReachableId = mGraph->getSegments()->size();
        for (Vertex * v : notReachableVertices) {
            if (v->getId() < mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    } else if (targetVertex == mGraph->getSink()->getNegativeVertex()) {
        // use the greatest id
        int mostReachableId = -1;
        for (Vertex * v : notReachableVertices) {
            if (v->getId() > mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    } else if (targetVertex == mGraph->getSource()->getNegativeVertex()) {
        // use the least id
        int mostReachableId = mGraph->getSegments()->size();
        for (Vertex * v : notReachableVertices) {
            if (v->getId() < mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    } else if (targetVertex == mGraph->getSink()->getPositiveVertex()) {
        // use the greatest id
        int mostReachableId = -1;
        for (Vertex * v : notReachableVertices) {
            if (v->getId() > mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    }
    return mostReachable;
}

bool LocalGenomicMap::adjustReachability(VertexPath & notReachableVertices, Vertex * targetVertex, JunctionDB * aJuncDB, bool verbose) {
    // sort(notReachableVertices.begin(), notReachableVertices.end(), [](Vertex *v1, Vertex *v2) { return v1->getWeight()->getCoverage() > v2->getWeight()->getCoverage(); });
    // for (Vertex * v: notReachableVertices) {
    //     cout << v->getInfo() << ":" << v->getWeight()->getCoverage() << " ";
    // }
    // cout << endl;
    int * nReachability = new int[notReachableVertices.size()];
    this->countReachability(nReachability, notReachableVertices, targetVertex);
    if (verbose) {
        for (int i = 0; i < notReachableVertices.size(); i++) {
            cout << nReachability[i] << " ";
        }
        cout << endl;
    }
    // int maxReachability = *max_element(nReachability, nReachability + notReachableVertices.size());
    // this->findWithReachability(nReachability, notReachableVertices, maxReachability);
    // if (verbose) {
    //     cout << "max: " << maxReachability << endl;
    //     this->print(notReachableVertices);
    // }
    // delete nReachability;

    // Vertex * mostReachable = this->getMostReachable(notReachableVertices, targetVertex);
    // if (verbose) cout << "Most reachable vertex: " << mostReachable->getInfo() << endl;
    // if (verbose) cout << "Most reachable vertex: " << mostReachable->getInfo() << endl;

    Record *rec;
    Vertex *mostReachable;
    bool hasAdded;
    for (VertexPath::iterator i = notReachableVertices.begin(); i != notReachableVertices.end(); i++) {
        hasAdded = false;
        mostReachable = *i;
        // if  (!mostReachable->getSegment()->hasLowerBoundLimit()) {
        //     cout << "Not bounded seg" << endl;
        //     continue;
        // }

        if (targetVertex == mGraph->getSource()->getPositiveVertex() || targetVertex == mGraph->getSink()->getNegativeVertex()) {
            cout << "Current vertex: " << mostReachable->getId() << " " << mostReachable->getSegment()->getChrom() << " " << mostReachable->getStart() << mostReachable->getDir() << endl;

            rec = aJuncDB->findRecord(mostReachable->getSegment()->getChrom(), mostReachable->getStart(), mostReachable->getDir());
            // cout << rec->getChrom() << " " << rec->getPos() << " " << rec->getStrand() << endl;
            Vertex *prevVertex;
            Junction *inferredJunc;
            // if (rec == NULL || rec->getBackwardEntries()->size() == 0) {
            if (rec != NULL && rec->getBackwardEntries()->size() > 0) {
            //     cout << "Record is NULL" << endl;
            //     if (mostReachable->getDir() == '+') {
            //         prevVertex = mGraph->getSegmentById(mostReachable->getId() - 1)->getPositiveVertex();
            //     } else {
            //         prevVertex = mGraph->getSegmentById(mostReachable->getId() + 1)->getNegativeVertex();
            //     }
            //     try {
            //         if (verbose)
            //             cout << "prev: " << prevVertex->getInfo() << endl;
            //         inferredJunc = mGraph->addJunction(prevVertex, mostReachable, this->inferCoverage(prevVertex, mostReachable), this->inferCredibility(prevVertex, mostReachable), true, false);
            //         if (verbose)
            //             cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //         hasAdded = true;
            //     } catch (DuplicateJunctionException& e) {
            //         cout << e.what() << endl;
            //     }
            // } else {
                for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin(); riter != rec->getBackwardEntries()->rend(); riter++) {
                    cout << "backward " << (*riter)->chrom << " " << (*riter)->pos << " " << (*riter)->strand << endl;
                    try {
                        if ((*riter)->strand == '+') {
                            prevVertex = mGraph->getSegmentByChromEnd((*riter)->chrom, (*riter)->pos)->getPositiveVertex();
                        } else {
                            prevVertex = mGraph->getSegmentByChromStart((*riter)->chrom, (*riter)->pos)->getNegativeVertex();
                        }
                        if (verbose)
                            cout << "prev: " << prevVertex->getInfo() << endl;
                        inferredJunc = mGraph->addJunction(prevVertex, mostReachable, this->inferCoverage(prevVertex, mostReachable), this->inferCredibility(prevVertex, mostReachable), -1, true, false, false);
                        if (inferredJunc == NULL) continue;
                        if (verbose)
                            cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                        hasAdded = true;
                        break;
                    } catch (SegmentDoesNotExistException& e) {
                        cout << e.what() << endl;
                    } catch (DuplicateJunctionException& e) {
                        cout << e.what() << endl;
                    }
                }
            }
            if (!hasAdded) {
                cout << "Record is NULL or no edge can be added from juncdb" << endl;
                // TODO: need refinement
                prevVertex = mostReachable;
                while (true) {
                if (mostReachable->getDir() == '+') {
                    // prevVertex = mGraph->getSegmentById(mostReachable->getId() - 1)->getPositiveVertex();
                    if (prevVertex->getId() == mGraph->getSource()->getId()) {
                        prevVertex = mGraph->getSink()->getPositiveVertex();
                        // break;
                    } else if (prevVertex->getId() > mGraph->getSink()->getId()) {
                        break;
                    } else {
                        prevVertex = mGraph->getSegmentById(prevVertex->getId() - 1)->getPositiveVertex();
                    }
                } else {
                    // prevVertex = mGraph->getSegmentById(mostReachable->getId() + 1)->getNegativeVertex();
                    if (prevVertex->getId() == mGraph->getSink()->getId()) {
                        prevVertex = mGraph->getSource()->getNegativeVertex();
                        // break;
                    } else if (prevVertex->getId() > mGraph->getSink()->getId()) {
                        break;
                    } else {
                        prevVertex = mGraph->getSegmentById(prevVertex->getId() + 1)->getNegativeVertex();
                    }
                }
                try {
                    if (verbose)
                        cout << "prev: " << prevVertex->getInfo() << endl;
                    inferredJunc = mGraph->addJunction(prevVertex, mostReachable, this->inferCoverage(prevVertex, mostReachable), this->inferCredibility(prevVertex, mostReachable), -1, true, false, false);
                    if (inferredJunc == NULL) continue;
                    if (verbose)
                        cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                    hasAdded = true;
                    break;
                } catch (DuplicateJunctionException& e) {
                    cout << e.what() << endl;
                }
                }
            }
        } else {
            rec = aJuncDB->findRecord(mostReachable->getSegment()->getChrom(), mostReachable->getEnd(), mostReachable->getDir());
            Vertex *nextVertex;
            Junction *inferredJunc;
            // if (rec == NULL || rec->getForwardEntries()->size() == 0) {
            if (rec != NULL && rec->getForwardEntries()->size() > 0) {
            //     cout << "Record is NULL" << endl;
            //     if (mostReachable->getDir() == '+') {
            //         nextVertex = mGraph->getSegmentById(mostReachable->getId() + 1)->getPositiveVertex();
            //     } else {
            //         nextVertex = mGraph->getSegmentById(mostReachable->getId() - 1)->getNegativeVertex();
            //     }
            //     try {
            //         if (verbose)
            //             cout << "prev: " << nextVertex->getInfo() << endl;
            //         inferredJunc = mGraph->addJunction(mostReachable, nextVertex, this->inferCoverage(mostReachable, nextVertex), this->inferCredibility(mostReachable, nextVertex), true, false);
            //         if (verbose)
            //             cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
            //         hasAdded = true;
            //     } catch (DuplicateJunctionException& e) {
            //         cout << e.what() << endl;
            //     }
            // } else {
                for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin(); riter != rec->getForwardEntries()->rend(); riter++) {
                    cout << "forward " << (*riter)->chrom << " " << (*riter)->pos << endl;
                    try {
                        if ((*riter)->strand == '+') {
                            nextVertex = mGraph->getSegmentByChromStart((*riter)->chrom, (*riter)->pos)->getPositiveVertex();
                        } else {
                            nextVertex = mGraph->getSegmentByChromEnd((*riter)->chrom, (*riter)->pos)->getNegativeVertex();
                        }
                        if (verbose)
                            cout << "prev: " << nextVertex->getInfo() << endl;
                        inferredJunc = mGraph->addJunction(mostReachable, nextVertex, this->inferCoverage(mostReachable, nextVertex), this->inferCredibility(mostReachable, nextVertex), -1, true, false, false);
                        if (inferredJunc == NULL) continue;
                        if (verbose)
                            cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                        hasAdded = true;
                        break;
                    } catch (SegmentDoesNotExistException& e) {
                        cout << e.what() << endl;
                    } catch (DuplicateJunctionException& e) {
                        cout << e.what() << endl;
                    }
                }
            }
            if (!hasAdded) {
                cout << "Record is NULL" << endl;
                // TODO: need refinement
                nextVertex = mostReachable;
                while (true) {
                if (mostReachable->getDir() == '+') {
                    // nextVertex = mGraph->getSegmentById(mostReachable->getId() + 1)->getPositiveVertex();
                    if (nextVertex->getId() == mGraph->getSink()->getId()) {
                        nextVertex = mGraph->getSource()->getPositiveVertex();
                        // break;
                    } else if (nextVertex->getId() > mGraph->getSink()->getId()) {
                        break;
                    } else {
                        nextVertex = mGraph->getSegmentById(nextVertex->getId() + 1)->getPositiveVertex();
                    }
                } else {
                    // nextVertex = mGraph->getSegmentById(mostReachable->getId() - 1)->getNegativeVertex();
                    if (nextVertex->getId() == mGraph->getSource()->getId()) {
                        nextVertex = mGraph->getSink()->getNegativeVertex();
                        // break;
                    } else if (nextVertex->getId() > mGraph->getSink()->getId()) {
                        break;
                    } else {
                        nextVertex = mGraph->getSegmentById(nextVertex->getId() - 1)->getNegativeVertex();
                    }
                }
                try {
                    if (verbose)
                        cout << "next: " << nextVertex->getInfo() << endl;
                    inferredJunc = mGraph->addJunction(mostReachable, nextVertex, this->inferCoverage(mostReachable, nextVertex), this->inferCredibility(mostReachable, nextVertex), -1, true, false, false);
                    if (inferredJunc == NULL) continue;
                    if (verbose)
                        cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
                    hasAdded = true;
                    break;
                } catch (DuplicateJunctionException& e) {
                    cout << e.what() << endl;
                }
                }
            }
        }
        if (hasAdded) break;
    }
    return hasAdded;

    // Record * rec = aJuncDB->findRecord(, mostReachable->getDir());
    // if (rec == NULL) {
    //     // no related records in database, add juntion to nearby segment
    // } else {
        // if (targetVertex == mGraph->getSource()->getPositiveVertex()) {
        //     // backward to source+
        //     Vertex * prevVertex;
        //     while (true) {
        //         try {
        //             Vertex * selectedVertex = this->selectPrevVertex(mostReachable, mGraph->getSource()->getPositiveVertex(), aJuncDB);
        //             if (selectedVertex != NULL) {
        //                 prevVertex = selectedVertex;
        //             } else {
        //                 cout << "No available record for backward source" << endl;
        //                 if (mostReachable->getDir() == '+') {
        //                     prevVertex = mGraph->getPrevVertexById(mostReachable);
        //                 } else {
        //                     prevVertex = mGraph->getPrevVertexById(mostReachable->getComplementVertex());
        //                 }
        //             }
        //             if (verbose) cout << "prev: " << prevVertex->getInfo() << endl;
        //             Junction * inferredJunc = mGraph->addJunction(prevVertex, mostReachable, this->inferCoverage(prevVertex, mostReachable), this->inferCredibility(prevVertex, mostReachable), true, false);
        //             if (verbose) cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
        //             break;
        //         } catch (DuplicateJunctionException& e) {
        //             cout << e.what() << endl;
        //             mostReachable = prevVertex;
        //             // prevVertex = mGraph->getPrevVertexById(mostReachable);
        //         }
        //     }
        // } else if (targetVertex == mGraph->getSink()->getNegativeVertex()) {
        //     // backward to sink-
        //     Vertex * prevVertex;
        //     while (true) {
        //         try {
        //             Vertex * selectedVertex = this->selectPrevVertex(mostReachable, mGraph->getSink()->getNegativeVertex(), aJuncDB);
        //             if (selectedVertex != NULL) {
        //                 prevVertex = selectedVertex;
        //             } else {
        //                 if (mostReachable->getDir() == '+') {
        //                     prevVertex = mGraph->getPrevVertexById(mostReachable->getComplementVertex());
        //                 } else {
        //                     prevVertex = mGraph->getPrevVertexById(mostReachable);
        //                 }
        //             }
        //             if (verbose) cout << "prev: " << prevVertex->getInfo() << endl;
        //             Junction * inferredJunc = mGraph->addJunction(prevVertex, mostReachable, this->inferCoverage(prevVertex, mostReachable), this->inferCredibility(prevVertex, mostReachable), true, false);
        //             if (verbose) cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
        //             break;
        //         } catch (DuplicateJunctionException& e) {
        //             cout << e.what() << endl;
        //             mostReachable = prevVertex;
        //             // prevVertex = mGraph->getPrevVertexById(mostReachable);
        //         }
        //     }
        // } else if (targetVertex == mGraph->getSink()->getPositiveVertex()) {
        //     // forward to sink+
        //     Vertex * nextVertex;
        //     while (true) {
        //         try {
        //             Vertex * selectedVertex = this->selectNextVertex(mostReachable, mGraph->getSink()->getPositiveVertex(), aJuncDB);
        //             if (selectedVertex != NULL) {
        //                 nextVertex = selectedVertex;
        //             } else {
        //                 if (mostReachable->getDir() == '+') {
        //                     nextVertex = mGraph->getNextVertexById(mostReachable);
        //                 } else {
        //                     nextVertex = mGraph->getNextVertexById(mostReachable->getComplementVertex());
        //                 }
        //             }
        //             if (verbose) cout << "prev: " << nextVertex->getInfo() << endl;
        //             Junction * inferredJunc = mGraph->addJunction(mostReachable, nextVertex, this->inferCoverage(mostReachable, nextVertex), this->inferCredibility(mostReachable, nextVertex), true, false);
        //             if (verbose) cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
        //             break;
        //         } catch (DuplicateJunctionException& e) {
        //             cout << e.what() << endl;
        //             mostReachable = nextVertex;
        //             // prevVertex = mGraph->getPrevVertexById(mostReachable);
        //         }
        //     }
        // } else if (targetVertex == mGraph->getSource()->getNegativeVertex()) {
        //     // forward to source-
        //     Vertex * nextVertex;
        //     while (true) {
        //         try {
        //             Vertex * selectedVertex = this->selectNextVertex(mostReachable, mGraph->getSource()->getNegativeVertex(), aJuncDB);
        //             if (selectedVertex != NULL) {
        //                 nextVertex = selectedVertex;
        //             } else {
        //                 if (mostReachable->getDir() == '+') {
        //                     nextVertex = mGraph->getNextVertexById(mostReachable->getComplementVertex());
        //                 } else {
        //                     nextVertex = mGraph->getNextVertexById(mostReachable);
        //                 }
        //             }
        //             if (verbose) cout << "prev: " << nextVertex->getInfo() << endl;
        //             Junction * inferredJunc = mGraph->addJunction(mostReachable, nextVertex, this->inferCoverage(mostReachable, nextVertex), this->inferCredibility(mostReachable, nextVertex), true, false);
        //             if (verbose) cout << "add junction: " << inferredJunc->getInfo()[0] << ", " << inferredJunc->getInfo()[1] << endl;
        //             break;
        //         } catch (DuplicateJunctionException& e) {
        //             cout << e.what() << endl;
        //             mostReachable = nextVertex;
        //             // prevVertex = mGraph->getPrevVertexById(mostReachable);
        //         }
        //     }
        // }
    // }
}

bool LocalGenomicMap::isReachable() {
    VertexPath backwardSourceNotReachableVertices;
    VertexPath backwardSinkNotReachableVertices;
    VertexPath forwardSourceNotReachableVertices;
    VertexPath forwardSinkNotReachableVertices;

        mGraph->checkOrphan();

        // check reachability for each vertex
        bool isSinkSourceConnectedOriginally;
        try {
            this->connectSourceSink();
            isSinkSourceConnectedOriginally = false;
        } catch (DuplicateJunctionException &e) {
            isSinkSourceConnectedOriginally = true;
        }
        for (Segment * seg: *(mGraph->getSegments())) {
            if (seg == mGraph->getSource() || seg == mGraph->getSink()) {
                continue;
            }

            // heuristic rule to keep a segment
            if (seg->isOrphan()) {
                if (seg->getWeight()->getCoverage() < 0.25 * mGraph->getAvgCoverage()) {
                    continue;   
                }
            }

            
            VertexPath segV;
            segV.push_back(seg->getPositiveVertex());
            segV.push_back(seg->getNegativeVertex());
            for (Vertex * v : segV) {
                // cout << "source+ -> " << v->getInfo() << endl;
                bool isBackwardSourceReachable = this->doesPathExists(mGraph->getSource()->getPositiveVertex(), v);
                // cout << "sink- -> " << v->getInfo() << endl;
                bool isBackwardSinkReachable = this->doesPathExists(mGraph->getSink()->getNegativeVertex(), v);
                // cout << v->getInfo() << " -> source-" << endl;
                bool isForwardSourceReachable = this->doesPathExists(v, mGraph->getSource()->getNegativeVertex());
                // cout << v->getInfo() << " -> sink+" << endl;
                bool isForwardSinkReachable = this->doesPathExists(v, mGraph->getSink()->getPositiveVertex());
                if (!isBackwardSourceReachable && !isForwardSinkReachable
                        && !isBackwardSinkReachable && !isForwardSourceReachable) {
                    if (v->getDir() == '+') {
                        backwardSourceNotReachableVertices.push_back(v);
                    } else {
                        backwardSinkNotReachableVertices.push_back(v);
                    }
                }
                if ((isBackwardSourceReachable ^ isForwardSinkReachable) != 0) {
                    if (!isBackwardSourceReachable) {
                        backwardSourceNotReachableVertices.push_back(v);
                    } else {
                        forwardSinkNotReachableVertices.push_back(v);
                    }
                }
                if ((isBackwardSinkReachable ^ isForwardSourceReachable) != 0) {
                    if (!isBackwardSinkReachable) {
                        backwardSinkNotReachableVertices.push_back(v);
                    } else {
                        forwardSourceNotReachableVertices.push_back(v);
                    }
                }
            }
        }
        if (!isSinkSourceConnectedOriginally) {
            delete mGraph->getJunctions()->back();
            mGraph->getJunctions()->pop_back();
        }
    
    return backwardSourceNotReachableVertices.size() + backwardSinkNotReachableVertices.size() + forwardSourceNotReachableVertices.size() + forwardSinkNotReachableVertices.size() == 0;
}

void LocalGenomicMap::checkReachability(JunctionDB * aJuncDB, bool verbose) {
    cout << "Checking reachability..." << endl;
    VertexPath backwardSourceNotReachableVertices;
    VertexPath backwardSinkNotReachableVertices;
    VertexPath forwardSourceNotReachableVertices;
    VertexPath forwardSinkNotReachableVertices;
    do {
        backwardSourceNotReachableVertices.clear();
        backwardSinkNotReachableVertices.clear();
        forwardSourceNotReachableVertices.clear();
        forwardSinkNotReachableVertices.clear();

        mGraph->checkOrphan();

        // check reachability for each vertex
        bool isSinkSourceConnectedOriginally;
        try {
            this->connectSourceSink();
            isSinkSourceConnectedOriginally = false;
        } catch (DuplicateJunctionException &e) {
            isSinkSourceConnectedOriginally = true;
        }
        for (Segment * seg: *(mGraph->getSegments())) {
            if (seg == mGraph->getSource() || seg == mGraph->getSink()) {
                continue;
            }

            // heuristic rule to keep a segment
            if (seg->isOrphan()) {
                if (seg->getWeight()->getCoverage() < 0.25 * mGraph->getAvgCoverage()) {
                    continue;   
                }
            }

            if (!seg->hasLowerBoundLimit()) continue;
            
            VertexPath segV;
            segV.push_back(seg->getPositiveVertex());
            segV.push_back(seg->getNegativeVertex());
            for (Vertex * v : segV) {
                // cout << "source+ -> " << v->getInfo() << endl;
                bool isBackwardSourceReachable = this->doesPathExists(mGraph->getSource()->getPositiveVertex(), v);
                // cout << "sink- -> " << v->getInfo() << endl;
                bool isBackwardSinkReachable = this->doesPathExists(mGraph->getSink()->getNegativeVertex(), v);
                // cout << v->getInfo() << " -> source-" << endl;
                bool isForwardSourceReachable = this->doesPathExists(v, mGraph->getSource()->getNegativeVertex());
                // cout << v->getInfo() << " -> sink+" << endl;
                bool isForwardSinkReachable = this->doesPathExists(v, mGraph->getSink()->getPositiveVertex());
                if (!isBackwardSourceReachable && !isForwardSinkReachable
                        && !isBackwardSinkReachable && !isForwardSourceReachable) {
                    if (v->getDir() == '+') {
                        backwardSourceNotReachableVertices.push_back(v);
                    } else {
                        backwardSinkNotReachableVertices.push_back(v);
                    }
                }
                if ((isBackwardSourceReachable ^ isForwardSinkReachable) != 0) {
                    if (!isBackwardSourceReachable) {
                        backwardSourceNotReachableVertices.push_back(v);
                    } else {
                        forwardSinkNotReachableVertices.push_back(v);
                    }
                }
                if ((isBackwardSinkReachable ^ isForwardSourceReachable) != 0) {
                    if (!isBackwardSinkReachable) {
                        backwardSinkNotReachableVertices.push_back(v);
                    } else {
                        forwardSourceNotReachableVertices.push_back(v);
                    }
                }
            }
        }
        if (!isSinkSourceConnectedOriginally) {
            delete mGraph->getJunctions()->back();
            mGraph->getJunctions()->pop_back();
        }
        // mGraph->print();
        
        if (verbose) {
            cout << "backwardSourceNotReachableVertices: ";
            this->print(backwardSourceNotReachableVertices);
            cout << "backwardSinkNotReachableVertices: ";
            this->print(backwardSinkNotReachableVertices);
            cout << "forwardSourceNotReachableVertices: ";
            this->print(forwardSourceNotReachableVertices);
            cout << "forwardSinkNotReachableVertices: ";
            this->print(forwardSinkNotReachableVertices);
            cout << endl;
        }

        bool hasAdded = false;
        if (!hasAdded && backwardSourceNotReachableVertices.size() > 0) {
            if (verbose) cout << "for backward source: " << endl;
            hasAdded = this->adjustReachability(backwardSourceNotReachableVertices, mGraph->getSource()->getPositiveVertex(), aJuncDB, verbose);
            // continue;
        }
        if (!hasAdded && backwardSinkNotReachableVertices.size() > 0) {
            if (verbose) cout << "for backward sink: " << endl;
            hasAdded = this->adjustReachability(backwardSinkNotReachableVertices, mGraph->getSink()->getNegativeVertex(), aJuncDB, verbose);
            // continue;
        }
        if (!hasAdded && forwardSourceNotReachableVertices.size() > 0) {
            if (verbose) cout << "for forward source: " << endl;
            hasAdded = this->adjustReachability(forwardSourceNotReachableVertices, mGraph->getSource()->getNegativeVertex(), aJuncDB, verbose);
            // continue;
        }
        if (!hasAdded && forwardSinkNotReachableVertices.size() > 0) {
            if (verbose) cout << "for forward sink: " << endl;
            hasAdded = this->adjustReachability(forwardSinkNotReachableVertices, mGraph->getSink()->getPositiveVertex(), aJuncDB, verbose);
        }
        // break;
    } while (backwardSourceNotReachableVertices.size() + backwardSinkNotReachableVertices.size() + forwardSourceNotReachableVertices.size() + forwardSinkNotReachableVertices.size() > 0);

    cout << "Graph reachable." << endl;
}

void LocalGenomicMap::addNormalJunctions() {
    for (int i = 0; i < mGraph->getSegments()->size() - 1; i++) {
        Vertex * sourceVertex = (*mGraph->getSegments())[i]->getPositiveVertex();
        Vertex * targetVertex = mGraph->getNextVertexById(sourceVertex);
        try {
            mGraph->addJunction(sourceVertex, targetVertex, this->inferCoverage(sourceVertex, targetVertex), inferCredibility(sourceVertex, targetVertex), -1, true, false, false);
        } catch(DuplicateJunctionException& e) {
            continue;
        }
    }
}

void LocalGenomicMap::checkInferredJunctionCredibility() {
    double maxCred = 1.0;
    cout << "Read" << endl;
    for (Junction * junc : *(mGraph->getJunctions())) {
        if (junc->isInferred()) {
            maxCred = max(pow(junc->getCredibility(), 0.5), maxCred);
        }
    }

    cout << "Assign" << endl;
    for (Junction * junc : *(mGraph->getJunctions())) {
        if (junc->isInferred()) {
            cout << junc->getInfo()[0] << endl;
            junc->setCredibility(pow(junc->getCredibility(), 0.5) / maxCred);
        }
    }
    cout << "Assign done" << endl;
}

void LocalGenomicMap::clearSegmentJunctionCredibility(Segment * aSegment) {
    cout << "Clear credibility of segment " << aSegment->getId() << " and its related junctions" << endl;
    aSegment->setCredibility(0);
    aSegment->resetHasLowerBoundLimit();
    for (Edge * e : *(aSegment->getPositiveVertex()->getEdgesAsSource())) {
        Junction * junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }
    
    for (Edge * e : *(aSegment->getPositiveVertex()->getEdgesAsTarget())) {
        Junction * junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }

    for (Edge * e : *(aSegment->getNegativeVertex()->getEdgesAsSource())) {
        Junction * junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }
    
    for (Edge * e : *(aSegment->getNegativeVertex()->getEdgesAsTarget())) {
        Junction * junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }
}

Vertex * LocalGenomicMap::selectPrevVertex(Vertex * currentVertex, Vertex * targetVertex, JunctionDB * aJuncDB) {
    Record *rec;
    rec = aJuncDB->findRecord(currentVertex->getSegment()->getChrom(), currentVertex->getStart(), currentVertex->getDir());
    Vertex * selectedVertex = NULL;
    for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin(); riter != rec->getBackwardEntries()->rend(); riter++) {
        Segment * seg = mGraph->getSegmentByChromEnd((*riter)->chrom, (*riter)->pos);
        if (targetVertex == mGraph->getSource()->getPositiveVertex()) {
            if ((*riter)->strand == '+') {
                // prefer the positive vertex that directly connects to + strand
                selectedVertex = seg->getPositiveVertex();
                break;
            } else {
                // if no proper positive, - vertex that can cross back to + strand
                int crossFlag = this->strandCross(seg->getNegativeVertex());
                if ((crossFlag & 0b10) != 0) {
                    // can cross to - strand
                    selectedVertex = seg->getNegativeVertex();
                    break;
                }
            }
        } else if (targetVertex == mGraph->getSink()->getNegativeVertex()) {
            if ((*riter)->strand == '-') {
                selectedVertex = seg->getNegativeVertex();
                break;
            } else {
                int crossFlag = this->strandCross(seg->getPositiveVertex());
                if ((crossFlag & 0b10) != 0) {
                    selectedVertex = seg->getPositiveVertex();
                    break;
                }
            }
        }
    }
    return selectedVertex;
}

Vertex * LocalGenomicMap::selectNextVertex(Vertex * currentVertex, Vertex * targetVertex, JunctionDB * aJuncDB) {
    Record *rec;
    rec = aJuncDB->findRecord(currentVertex->getSegment()->getChrom(), currentVertex->getEnd(), currentVertex->getDir());
    Vertex * selectedVertex = NULL;
    for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin(); riter != rec->getForwardEntries()->rend(); riter++) {
        Segment * seg = mGraph->getSegmentByChromStart((*riter)->chrom, (*riter)->pos);
        if (targetVertex == mGraph->getSink()->getPositiveVertex()) {
            if ((*riter)->strand == '+') {
                // prefer positive vertex
                selectedVertex = seg->getPositiveVertex();
                break;
            } else {
                int crossFlag = this->strandCross(seg->getNegativeVertex());
                if ((crossFlag & 0b01) != 0) {
                    // can cross to + strand
                    selectedVertex = seg->getNegativeVertex();
                    break;
                }
            }
        } else if (targetVertex == mGraph->getSource()->getNegativeVertex()) {
            if ((*riter)->strand == '-') {
                selectedVertex = seg->getNegativeVertex();
                break;
            } else {
                int crossFlag = this->strandCross(seg->getPositiveVertex());
                if ((crossFlag & 0b01) != 0) {
                    selectedVertex = seg->getPositiveVertex();
                    break;
                }
            }
        }
    }
    return selectedVertex;
}

Edge * LocalGenomicMap::selectPrevEdge(Vertex * aTargetVertex, bool isTraversing) {
    for (Edge * e : *(aTargetVertex->getEdgesAsTarget())) {
        if (!e->isVisited()) {
            if (isTraversing) {
                if (e->hasCopy()) {
                    return e;
                } else {
                    continue;
                }
            } else {
                return e;
            }
        }
        // if (skipVisited) {
        //     if (!e->isVisited()) {
        //         return e;
        //     }
        // } else {
        //     if (e->hasCopy()) {
        //         return e;
        //     }
        // }
    }
    return NULL;
}

Edge * LocalGenomicMap::selectNextEdge(Vertex * aSourceVertex, bool isTraversing) {
    for (Edge * e : *(aSourceVertex->getEdgesAsSource())) {
        if (!e->isVisited()) {
            if (isTraversing) {
                if (e->hasCopy()) {
                    return e;
                } else {
                    continue;
                }
            } else {
                return e;
            }
        }
        // if (skipVisited) {
        //     if (!e->isVisited()) {
        //         return e;
        //     }
        // } else {
        //     if (e->hasCopy()) {
        //         return e;
        //     }
        // }
    }
    return NULL;
}

int LocalGenomicMap::findCircuit(Vertex * aVertex, VertexPath & pathVertices, EdgePath & pathEdges) {
    cout << "Search start: " << aVertex->getInfo() << endl;
    Vertex * currentVertex = aVertex;
    Edge * currentEdge = this->selectNextEdge(currentVertex, true);

    pathVertices.push_back(currentVertex);
    currentVertex->traverse();
    while (true) {
        if (currentEdge == NULL) {
            // Not found available edge to next vertex
            // return to last vertex and search again
            currentVertex->recover();  // recover the copy
            pathVertices.pop_back();
            
            if (pathVertices.size() == 0) {
                // no vertex left that can continue being used
                return -1;
            } else {
                pathEdges.back()->recover();
                pathEdges.pop_back();

                currentVertex = pathVertices.back();
                currentEdge = this->selectNextEdge(currentVertex, true);
            }
        } else {
            // cout << "Attempt: " << currentEdge->getInfo() << endl;
            Vertex * nextVertex = currentEdge->getTarget();
            
            if (nextVertex == aVertex) {
                // found
                // DO NOT traverse the vertex again
                pathVertices.push_back(nextVertex);

                currentEdge->traverse();
                pathEdges.push_back(currentEdge);
                return 0;
            } else {
                currentVertex = nextVertex;

                currentVertex->traverse();
                currentEdge->setVisited();
                currentEdge->traverse();

                pathVertices.push_back(currentVertex);
                pathEdges.push_back(currentEdge);

                currentEdge = this->selectNextEdge(currentVertex, true);
            }
        }
    }

    // Vertex * startVertex = aStartVertex;
    // VertexPath vertexStack;
    // vertexStack.push_back(startVertex);
    // Edge * nextEdge = this->selectNextEdge(startVertex, true);
    // Edge * lastEdge = NULL;
    // if (nextEdge != NULL) {
    //     cout << "Traverse " << startVertex->getInfo() << endl;
    //     startVertex->traverse();
    // }
    // while (true) {
    //     if (nextEdge == NULL) {
    //         vertexStack.pop_back();
    //         if (vertexStack.size() == 0) {
    //             return -1;  // not found
    //         } else {
    //             pathVertices.pop_back();
    //             pathEdges.pop_back();

    //             cout << "Recover " << startVertex->getInfo() << endl;
    //             startVertex->recover();
    //             if (lastEdge != NULL) {
    //                 cout << "Recover " << lastEdge->getInfo() << endl;
    //                 lastEdge->recover();
    //                 lastEdge = NULL;
    //             }

    //             startVertex = vertexStack.back();
    //             nextEdge = this->selectNextEdge(startVertex, true);
    //         }
    //     } else {
    //         Vertex * nextVertex = nextEdge->getTarget();
    //         cout << "Attempt: " << nextEdge->getInfo() << endl;
    //         pathVertices.push_back(startVertex);
    //         pathEdges.push_back(nextEdge);
    //         cout << "Traverse " << nextEdge->getInfo() << endl;
    //         nextEdge->traverse();
    //         
    //         if (nextVertex == aEndVertex) {
    //             // DO NOT traverse the vertex
    //             pathVertices.push_back(nextVertex);
    //             return 0;
    //         } else {
    //             cout << "Traverse " << nextVertex->getInfo() << endl;
    //             nextVertex->traverse();

    //             nextEdge->setVisited();
    //             vertexStack.push_back(nextVertex);
    //             startVertex = nextVertex;
    //             lastEdge = nextEdge;
    //             nextEdge = this->selectNextEdge(startVertex, true);
    //         }
    //     }
    // }
}

// void LocalGenomicMap::traverse(Vertex * startVertex) {
//     VertexPath * pathV = new VertexPath();
//     Vertex * currentVertex = startVertex;
//     Edge * nextEdge = this->selectNextEdge(currentVertex, true);
//     while (nextEdge != NULL) {
//         pathV->push_back(currentVertex);
//         // cout << currentVertex->getInfo() << " ";
//         currentVertex->traverse();
//         currentVertex = nextEdge->getTarget();
//         nextEdge->traverse();
//         nextEdge = this->selectNextEdge(currentVertex, true);
//     }
//     pathV->push_back(currentVertex);
//     // cout << currentVertex->getInfo() << endl;
//     this->print(*pathV);
//     mCircuits->push_back(pathV);
// }

Edge * LocalGenomicMap::traverseNextEdge(Vertex * aStartVertex, JunctionDB * aJuncDB)  {
    Edge * selectedEdge = NULL;
    Record * rec = aJuncDB->findRecord(aStartVertex->getSegment()->getChrom(), aStartVertex->getEnd(), aStartVertex->getDir());
    if (rec == NULL) {
        for (Edge * e: *(aStartVertex->getEdgesAsSource())) {
            if (e->hasCopy()) {
                selectedEdge = e;
                break;
            }
        }
    } else {
        int support = 0;
        entry_t * entry;
        for (Edge * e: *(aStartVertex->getEdgesAsSource())) {
            if (!e->hasCopy()) continue;
            entry = rec->findForwardEntry(e->getTarget()->getSegment()->getChrom(), e->getTarget()->getStart(), e->getTarget()->getDir());
            if (entry != NULL) {
                if (entry->support > support) {
                    support = entry->support;
                    selectedEdge = e;
                }
            } else {
                if (support == 0) {
                    selectedEdge = e;
                }
            }
        }
    }

    return selectedEdge;
}

void LocalGenomicMap::traverse(Vertex * aStartVertex, JunctionDB * aJuncDB) {
    VertexPath *vp = new VertexPath();
    EdgePath *ep = new EdgePath();
    Vertex * currentVertex;
    if (mLongFrags == nullptr) {
        currentVertex = aStartVertex;
    } else {
        currentVertex = traverseLongPath(aStartVertex, vp);
    }
    vp->push_back(currentVertex);
    Edge * nextEdge = this->traverseNextEdge(currentVertex, aJuncDB);
    while (nextEdge != NULL) {
        if (currentVertex->getId() == 24) {
            cout << mGraph->getSegmentById(24)->getWeight()->getCopyNum() << endl;
        }
        currentVertex->traverse();
        // if (nextEdge->getSource()->getSegment() == nextEdge->getTarget()->getSegment()) {
        //     nextEdge->traverse();
        //     nextEdge->traverse();
        // } else {
            nextEdge->traverse();
        // }
        ep->push_back(nextEdge);
        vp->push_back(nextEdge->getTarget());
        // if (currentVertex->getId() == 24) {
        //     int vcount = 0;
        //     for (VertexPath::iterator it = vp->begin(); it < vp->end() - 1; it++) {
        //         if ((*it)->getId() == 24) {
        //             vcount++;
        //         }
        //     }
        //     this->print(*vp);
        //     cout << "\033[1;31m" << vcount << "\033[0m " << mGraph->getSegmentById(24)->getWeight()->getCopyNum() << endl;
        // }

        currentVertex = nextEdge->getTarget();
        nextEdge = this->traverseNextEdge(currentVertex, aJuncDB);
    }

    mCircuits->push_back(vp);

    cout << "Traversed path: ";
    for (Vertex * v: *vp) {
        cout << v->getInfo() << " ";
    }
    cout << endl;
}

void LocalGenomicMap::traverseGraph(JunctionDB * aJuncDB) {
    Vertex * currentVertex;
    while (!mGraph->isCopyExhaustive()) {
        for (Segment * seg : *(mGraph->getSegments())) {
            if (seg->hasCopy()) {
                currentVertex = seg->getPositiveVertex();
                break;
            }
        }
        this->traverse(currentVertex, aJuncDB);
        // mGraph->print();
    }
}
Vertex * LocalGenomicMap::traverseLongPath(Vertex *aStartVertex, VertexPath* vPath) {
    int pathN = -1;
    int maxL = 0;
    int i = 0;
    for(auto path: *this->mLongFrags) {
        if (path[0][0] == aStartVertex) {
            int l = longPathLenInGraph(path);
            if (l > maxL) {
                maxL = l;
                pathN = i;
            }
        }
        i++;
    }
    i = 0;
    if (maxL <= 0) {
        return aStartVertex;
    }
//    auto m = this->mLongFrags[pathN][0];
    for (i = 0; i < maxL - 1; i++) {
        auto v = (*((*this->mLongFrags)[pathN]))[i];
        vPath->push_back(v);
        v->traverse();
        auto vNext = (*((*this->mLongFrags)[pathN]))[i+1];
        for(Edge* e: *(v->getEdgesAsSource())) {
            if (e->getTarget() == vNext) {
                e->traverse();
                break;
            }
        }
    }
    return (*(this->mLongFrags)[0][pathN])[maxL-1];
}

int LocalGenomicMap::longPathLenInGraph(VertexPath *longPath) {
    int n = 1;
    for (auto v = longPath->begin(); v != longPath->end()-1;v++) {
        bool flag = false;
//        auto t = *(v+1);
        for (Edge * e : *(*v)->getEdgesAsSource()) {
            if (e->getTarget() == *(v+1)) {
                flag = true;
                break;
            }
        }
        if(flag) n++;
    }
    return n;
}

// void LocalGenomicMap::traverseGraph() {
//     while (!mGraph->isCopyExhaustive()) {
//         for (Segment * seg : *(mGraph->getSegments())) {
//             if (!seg->hasCopy()) continue;
//             // this->traverse(seg->getPositiveVertex());
//             VertexPath * pathV = new VertexPath();
//             EdgePath pathE;
//             if (seg == mGraph->getSource()) {
//                 this->BFS(mGraph->getSource()->getPositiveVertex());
//                 int flag = this->findShortestPath(mGraph->getSource()->getPositiveVertex(), mGraph->getSink()->getPositiveVertex(), *pathV, pathE);
//                 mGraph->resetVertexVisitFlag();
//                 if (flag < 0) {
//                     for (Vertex * v : *pathV) v->recover();
//                     for (Edge * e : pathE) e->recover();
// 
//                     pathV->clear();
//                     pathE.clear();
//                     flag = this->findCircuit(seg->getPositiveVertex(), *pathV, pathE);
//                     mGraph->resetJunctionVisitFlag();
//                     if (flag < 0) {
//                         // for (Vertex * v : *pathV) v->recover();
//                         // for (Edge * e : pathE) e->recover();
//                         continue;
//                     }
//                 }
//             } else {
//                 int flag = this->findCircuit(seg->getPositiveVertex(), *pathV, pathE);
//                 mGraph->resetJunctionVisitFlag();
//                 if (flag < 0) {
//                     // for (Vertex * v : *pathV) v->recover();
//                     // for (Edge * e : pathE) e->recover();
//                     continue;
//                 }
//             }
//             if (pathV->size() > 0) {
//                 mCircuits->push_back(pathV);
//             }
//             if (seg == mGraph->getSource()) {
//                 mGraph->print();
//             }
//             break;
//         }
//     }
//     // for (VertexPath * pathV : *mCircuits) {
//     //     this->print(*pathV);
//     // }
//     // mGraph->print();
// }

void LocalGenomicMap::isCircuitSimple(VertexPath * circuit, pair<int, int> & notSimpleIdx) {
    notSimpleIdx = make_pair(-1, -1);
    VertexPath::iterator it;
    for (int i = 0; i < circuit->size(); i++) {
        it = find(circuit->begin() + i + 1, circuit->end(), (*circuit)[i]);
        if (i == 0 && it == circuit->end() - 1) continue;
        if (it != circuit->end()) {
            notSimpleIdx.first = i;
            notSimpleIdx.second = it - circuit->begin();
            break;
        }
    }
}

void LocalGenomicMap::allCircuitsSimple(vector< tuple<int, int, int> > & notSimpleIdx) {
    for (int i = 0; i < mCircuits->size(); i++) {
        pair<int, int> notSimpleVertexIdx;
        this->isCircuitSimple((*mCircuits)[i], notSimpleVertexIdx);
        if (notSimpleVertexIdx.first >= 0) {
            notSimpleIdx.push_back(make_tuple(i, notSimpleVertexIdx.first, notSimpleVertexIdx.second));
        }
    }
}

void LocalGenomicMap::extractCircuits() {
    vector< tuple<int, int, int> > notSimpleCircuitIdx;
    this->allCircuitsSimple(notSimpleCircuitIdx);
    while (notSimpleCircuitIdx.size() > 0) {
        for (tuple<int, int, int> idxPair : notSimpleCircuitIdx) {
            int circuitIdx, begin, end;
            tie(circuitIdx, begin, end) = idxPair;
            VertexPath * subCircuit = new VertexPath((*mCircuits)[circuitIdx]->begin() + begin, (*mCircuits)[circuitIdx]->begin() + end + 1);
            mCircuits->push_back(subCircuit);
            (*mCircuits)[circuitIdx]->erase((*mCircuits)[circuitIdx]->begin() + begin + 1, (*mCircuits)[circuitIdx]->begin() + end + 1);
        }
        notSimpleCircuitIdx.clear();
        this->allCircuitsSimple(notSimpleCircuitIdx);
    }
}

void LocalGenomicMap::sortCircuits() {
    sort(mCircuits->begin(), mCircuits->end(), [](VertexPath * c1, VertexPath * c2) { return (*c1)[0]->getId() < (*c2)[0]->getId(); });
}

void LocalGenomicMap::writeCircuits(const char * outFn) {
    cout << "Write circuits" << endl;
    ofstream fout(outFn);
    if (!fout) {
        cout << "Cannot open file " << outFn << ": no such file or directory" << endl;
        exit(1);
    }
    for (VertexPath * pathV : *mCircuits) {
        for (Vertex * v : *pathV) {
            fout << v->getInfo() << " ";
        }
        fout << endl;
    }
    fout.close();
}

void LocalGenomicMap::generateHaploids() {
    cout << "sort circuits" << endl;
    this->sortCircuits();
    srand(time(NULL));

    vector<bool> is_inserted(mCircuits->size(), false);
    is_inserted[0] = true;
    VertexPath *mainPath = (*mCircuits)[0];
    int i = 1;
    VertexPath * current_circuit = new VertexPath();
    VertexPath * comp_circuit;
    bool is_comp = false;
    bool all_inserted = false;
    while (!all_inserted) {
        while (i < mCircuits->size()) {
            if (is_inserted[i]) {
                i++;
                continue;
            }
        // for (int i = 1; i < mCircuits->size(); i++) {
            if (!is_comp) {
                current_circuit->assign(mCircuits->at(i)->begin(), mCircuits->at(i)->end());
            } else {
                comp_circuit = this->get_complement(mCircuits->at(i));
                current_circuit->assign(comp_circuit->begin(), comp_circuit->end());
                // comp_circuit->clear();
            }
            deque<Vertex *> vq = deque<Vertex *>(current_circuit->begin(), current_circuit->end() - 1);
            bool isInserted = false;
            VertexPath::iterator cItr;
            VertexPath::iterator foundItr;
            for (int j = 0; j <= vq.size(); j++) {
                Vertex * startV = vq.front();
                // while (!isInserted) {
                cItr = mainPath->begin();
                // cout << "i=" << i << ": " << startV->getInfo() << " " << (*cItr)->getInfo() << endl;
                while (cItr != mainPath->end() && cItr - 1 != mainPath->end()) {
                    foundItr = find(cItr, mainPath->end(), startV);
                    if (foundItr != mainPath->end()) {
                    // if (rand() * 1.0 / RAND_MAX < 0.5 && foundItr != mainPath->end()) {
                        cout << "Insert at " << (*foundItr)->getInfo() << endl;
                        isInserted = true;
                        break;
                    }
                    cItr = foundItr + 1;
                }
                // }
                if (isInserted) {
                    // mainPath->insert(foundItr, vq.begin(), vq.end());
                    // int count = 0;
                    // for (Vertex * v: *mainPath) {
                    //     if (v->getId() == 24) {
                    //         count++;
                    //     }
                    // }
                    // this->print(*mainPath);
                    // cout << count << endl;
                    break;
                    // mainPath->insert(foundItr + 1, (*mCircuits)[i]->begin() + 1, (*mCircuits)[i]->end());
                // } else {
                //     mHaploids->push_back((*mCircuits)[i]);
                }
                vq.pop_front();
                vq.push_back(startV);
            }
            if (isInserted) {
                mainPath->insert(foundItr, vq.begin(), vq.end());
                is_inserted[i] = true;
                i++;
                is_comp = false;
            } else {
                if (is_comp) {
                    i++;
                    is_comp = false;
                } else {
                    is_comp = true;
                }
                // cout << "NOT INSERTED: " << endl;
                // this->print(*(*mCircuits)[i]);
            }
        }
        all_inserted = true;
        for (i = 0; i < is_inserted.size(); i++) {
            if (!is_inserted[i]) {
                this->print(*mainPath);
                this->print(*mCircuits->at(i));
                all_inserted = false;
                break;
            }
        }
    }
    mHaploids->push_back(mainPath);
    cout << "Mainpath done" << endl;
    cout << mHaploids->size() << endl;

    for (Vertex * v: *mainPath) {
        cout << v->getInfo() << " ";
    }
    cout << endl;

    // VertexPath::iterator itr = find(mainPath->begin(), mainPath->end(), mGraph->getSink()->getPositiveVertex());
    // mHaploids->push_back(new VertexPath(mainPath->begin(), itr + 1));
    // mHaploids->push_back(new VertexPath(itr + 1, mainPath->end() - 1));
}

// void LocalGenomicMap::generateHaploids() {
//     this->sortCircuits();
// 
//     mHaploids->push_back(new VertexPath());    // by default, take first two circuits as main ones
//     mHaploids->push_back(new VertexPath());
//     *((*mHaploids)[0]) = *((*mCircuits)[0]);
//     *((*mHaploids)[1]) = *((*mCircuits)[1]);
// 
//     VertexPath * firstMain = (*mHaploids)[0];
//     VertexPath * secondMain = (*mHaploids)[1];
//     for (int i = 2; i < mCircuits->size(); i++) {
//         VertexPath::iterator foundOnFirst = find(firstMain->begin(), firstMain->end(), (*mCircuits)[i]->front());
//         VertexPath::iterator foundOnSecond = find(secondMain->begin(), secondMain->end(), (*mCircuits)[i]->front());
//         if (foundOnFirst != firstMain->end() && foundOnSecond != secondMain->end()) {
//             // seen on both main circuits
//             if (firstMain->size() <= secondMain->size()) {
//                 // insert to first
//                 firstMain->insert(foundOnFirst, (*mCircuits)[i]->begin(), (*mCircuits)[i]->end() - 1);
//             } else {
//                 secondMain->insert(foundOnSecond, (*mCircuits)[i]->begin(), (*mCircuits)[i]->end() - 1);
//             }
//         } else {
//             if (foundOnFirst != firstMain->end()) {
//                 firstMain->insert(foundOnFirst, (*mCircuits)[i]->begin(), (*mCircuits)[i]->end() - 1);
//             }
//             if (foundOnSecond != secondMain->end()) {
//                 secondMain->insert(foundOnSecond, (*mCircuits)[i]->begin(), (*mCircuits)[i]->end() - 1);
//             }
//         }
//     }
// }

void LocalGenomicMap::writeHaploids(const char * outFn) {
    ofstream fout(outFn);
    if (!fout) {
        cout << "Cannot open file " << outFn << " : no such file or directory" << endl;
        exit(1);
    }
    for (VertexPath * pathV : *mHaploids) {
        for (Vertex * v : *pathV) {
            fout << v->getInfo() << " ";
        }
        fout << endl;
    }
    fout.close();
}

void LocalGenomicMap::print(VertexPath& path) {
    for (Vertex * v : path) {
        cout << v->getInfo() << " ";
    }
    cout << endl;
}

void LocalGenomicMap::print(EdgePath& path) {
    for (Edge * e : path) {
        cout << e->getInfo() << " ";
    }
    cout << endl;
}

void LocalGenomicMap::printCircuits() {
    int count = 1;
    for (VertexPath * pathV : *mCircuits) {
        cout << count << " ---- ";
        count++;
        this->print(*pathV);
    }
}

void LocalGenomicMap::printHaploids() {
    for (VertexPath * pathV : *mHaploids) {
        this->print(*pathV);
    }
}
