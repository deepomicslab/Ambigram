#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <unistd.h>
#include <deque>
#include <set>
#include <random>
#include <stack>
#include <map>

#include "LocalGenomicMap.hpp"
#include "Exceptions.hpp"
#include <coin/CbcModel.hpp>
#include <coin/OsiClpSolverInterface.hpp>

using namespace std;

LocalGenomicMap::LocalGenomicMap(Graph *aGraph) {
    mGraph = aGraph;

    mCircuits = new vector<VertexPath *>();
    mHaploids = new vector<VertexPath *>();
    dividedCircuits = new unordered_map<int, vector<VertexPath *> *>();
    dividedHaploids = new unordered_map<int, vector<VertexPath *> *>();
    traversedCircuits = new unordered_map<int, vector<VertexPath *> *>();
    usingLong = false;
    usingHic = false;
}

LocalGenomicMap::~LocalGenomicMap() { ; }

Graph *LocalGenomicMap::getGraph() { return mGraph; }

void LocalGenomicMap::setGraph(Graph *aGraph) { mGraph = aGraph; }

bool LocalGenomicMap::isUsingLong() { return usingLong; }

bool LocalGenomicMap::isUsingHic() { return usingHic; }

void LocalGenomicMap::setUsingLong(bool v) { usingLong = v; }

void LocalGenomicMap::setUsingHic(bool v) { usingHic = v; }

unordered_map<int, vector<VertexPath *> *> *LocalGenomicMap::get_long_frags() {
    return mLongFrags;
}

void LocalGenomicMap::read_long_frags(const char *fn) {
//    init long path
    mLongFrags = new unordered_map<int, vector<VertexPath *> *>();
    for(auto i: *mGraph->getMSources()) {
        (*mLongFrags)[i->getId()] = new vector<VertexPath *>();
    }
//    (*mLongFrags)[0] = new vector<VertexPath *>();
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
//        check parition
        auto partitionPair = findLongPathPartition(p);
        if(partitionPair.first == -1) {
//            cout<<  <<"\n";
            continue;
        }
//        (*mLongFrags)
        if (p->front()->getDir() == '-') {
            if (p->back()->getDir() == '-') {
                (*mLongFrags)[partitionPair.first]->push_back(this->get_complement(p));
            } else {
                if (p->back()->getId() < p->front()->getId()) {
                    (*mLongFrags)[partitionPair.first]->push_back(this->get_complement(p));
                }
            }
            p->clear();
        } else {
            (*mLongFrags)[partitionPair.first]->push_back(p);
        }
    }

    in.close();


    for(auto &it : *mLongFrags) {
        sort(it.second->begin(), it.second->end(), [](VertexPath *p1, VertexPath *p2) {
            int size = min(p1->size(), p2->size());
            for (int i = 0; i < size; i++) {
                if (p1->at(i)->getId() != p2->at(i)->getId()) {
                    return p1->at(i)->getId() < p2->at(i)->getId();
                }
            }
        });
        vector<VertexPath *> *long_frags_new = new vector<VertexPath *>();
        vector<VertexPath *> *temp;
        while (!this->equal_frags(it.second, long_frags_new)) {
            temp = it.second;
            it.second = long_frags_new;
            long_frags_new = temp;
            it.second->clear();
            for (VertexPath *frag : *long_frags_new) {
                this->merge_long_frags(it.second,frag);
            }
            sort(it.second->begin(), it.second->end(), [](VertexPath *p1, VertexPath *p2) {
                int size = min(p1->size(), p2->size());
                for (int i = 0; i < size; i++) {
                    if (p1->at(i)->getId() != p2->at(i)->getId()) {
                        return p1->at(i)->getId() < p2->at(i)->getId();
                    }
                }
            });
        }
        sort(it.second->begin(), it.second->end(), [](VertexPath *p1, VertexPath *p2) {
            return p1->size() > p2->size();
        });
    }
//    while (!this->equal_frags(mLongFrags, long_frags_new)) {
//        temp = mLongFrags;
//        mLongFrags = long_frags_new;
//        long_frags_new = temp;
//        mLongFrags->clear();
//        for (VertexPath *frag : *long_frags_new) {
//            this->merge_long_frags(frag);
//        }
//        sort(mLongFrags->begin(), mLongFrags->end(), [](VertexPath *p1, VertexPath *p2) {
//            int size = min(p1->size(), p2->size());
//            for (int i = 0; i < size; i++) {
//                if (p1->at(i)->getId() != p2->at(i)->getId()) {
//                    return p1->at(i)->getId() < p2->at(i)->getId();
//                }
//            }
//        });
//    }
//    sort(mLongFrags->begin(), mLongFrags->end(), [](VertexPath *p1, VertexPath *p2) {
//        return p1->size() > p2->size();
//    });
}

void LocalGenomicMap::read_hic_matrix(const char *fn) {
    int len = this->mGraph->getSegments()->size();
    this->hicMatrix = new double *[len + 1];
    this->decreaseMatrix = new double *[len + 1];
    for (int i = 0; i < len + 1; i++) {
        this->hicMatrix[i] = new double[len + 1];
    }
    for (int i = 0; i < len + 1; i++) {
        this->decreaseMatrix[i] = new double[len + 1];
    }
    ifstream in(fn);
    string line, value;
    stringstream ss;
//    The first line is the copy number information
    int *copys = new int[len + 1];
    getline(in, line);
    ss = stringstream(line);
    while (ss >> value) {
        int copy = stoi(value);
        *(copys++) = copy;
    }
//    interactions
    int i = 0;
    int j = 0;
    while (getline(in, line)) {
        ss = stringstream(line);
        while (ss >> value) {
            double matrixV = stod(value);
            this->hicMatrix[i][j++] = matrixV;
            this->decreaseMatrix[i][j++] = matrixV / (copys[i]);
//            *(++(*(this->hicMatrix))) = matrixV;
//            cout<<**(this->hicMatrix)<<endl;
//            *(++(*(this->decreaseMatrix))) = matrixV/(*(++copys));
        }
        i++;
        j = 0;
    }
}

void LocalGenomicMap::merge_long_frags(vector<VertexPath *>* partitionLong,VertexPath *aFrag) {
    VertexPath::iterator found_iter, iter;
    bool found = false;
    for (iter = aFrag->end(); iter != aFrag->begin(); iter--) {
        for (int i = 0; i < partitionLong->size(); i++) {
            found_iter = find_end(partitionLong->at(i)->begin(), partitionLong->at(i)->end(), aFrag->begin(), iter);
            if (found_iter != partitionLong->at(i)->end()) {
                found = true;
                VertexPath *new_f = new VertexPath(partitionLong->at(i)->begin(), found_iter);
                new_f->insert(new_f->end(), aFrag->begin(), aFrag->end());
                found_iter = find_end(partitionLong->at(i)->begin(), partitionLong->at(i)->end(), new_f->begin(),
                                      new_f->end());
                if (found_iter == partitionLong->at(i)->end()) {
                    // new_f in
                    //     break;
                    // } else {
                    // new_f not in
                    found_iter = find_end(new_f->begin(), new_f->end(), partitionLong->at(i)->begin(),
                                          partitionLong->at(i)->end());
                    if (found_iter != new_f->end()) {
                        // new_f include
                        partitionLong->at(i)->clear();
                        partitionLong->at(i)->insert(partitionLong->at(i)->end(), new_f->begin(), new_f->end());
                    } else {
                        // new_f not include
                        partitionLong->push_back(new_f);
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
            for (int i = 0; i < partitionLong->size(); i++) {
                found_iter = find_end(partitionLong->at(i)->begin(), partitionLong->at(i)->end(), comp->begin(), iter);
                if (found_iter != partitionLong->at(i)->end()) {
                    found = true;
                    VertexPath *new_f = new VertexPath(partitionLong->at(i)->begin(), found_iter);
                    new_f->insert(new_f->end(), comp->begin(), comp->end());
                    found_iter = find_end(partitionLong->at(i)->begin(), partitionLong->at(i)->end(), new_f->begin(),
                                          new_f->end());
                    if (found_iter == partitionLong->at(i)->end()) {
                        // new_f in
                        //     break;
                        // } else {
                        // new_f not in
                        found_iter = find_end(new_f->begin(), new_f->end(), partitionLong->at(i)->begin(),
                                              partitionLong->at(i)->end());
                        if (found_iter != new_f->end()) {
                            // new_f include
                            partitionLong->at(i)->clear();
                            partitionLong->at(i)->insert(partitionLong->at(i)->end(), new_f->begin(), new_f->end());
                        } else {
                            // new_f not include
                            partitionLong->push_back(new_f);
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
        partitionLong->push_back(aFrag);
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
              [mean](double v) { return v - mean + 1; });
    vector<double> diff_s(diff.size());
    transform(diff.begin(), diff.end(), diff_s.begin(),
              [](double v) { return pow(v, 2); });
    double stdev = sqrt(accumulate(diff_s.begin(), diff_s.end(), 0.0) / diff_s.size());
    vector<double> scaled_cov(aCovs.size());
    transform(diff.begin(), diff.end(), scaled_cov.begin(),
              [stdev](double v) { return abs(v / stdev); });
    return scaled_cov;
}

int LocalGenomicMap::balancerILP(const char *lpFn) {
    OsiClpSolverInterface *si = new OsiClpSolverInterface();
    int sourceSize = this->getGraph()->getMSinks()->size() - 1;
//    int sourceSize = 2;

    // variable structure:
    // {segments (nSeg), junctions (nJunc), segment epsilon, junction epsilon}
    vector<Segment *> *segs = mGraph->getSegments();
    vector<Junction *> *juncs = mGraph->getJunctions();
    int sourceSinkNum = 2 * (sourceSize);


    int numSegsJuncs = segs->size() + juncs->size();
    int numEpsilons = numSegsJuncs;    // e(t(seg), c(seg))

    int numVariables = numSegsJuncs + juncs->size() + numEpsilons;
    int numConstrains = segs->size() * 4 + 4 * juncs->size();  // "2" for source & sink      ------  t(seg)+c(seg), t(seg)-c(seg), in=t(seg), out=t(seg, t(junc)+c(junc), t(junc)-c(junc), e, source, sink
    // int numConstrains = segs->size() * 8 + juncs->size() * 6 + numEpsilons + 2;

    double *objective = new double[numVariables];
    double *variableLowerBound = new double[numVariables];
    double *variableUpperBound = new double[numVariables];
    double *constrainLowerBound = new double[numConstrains];
    double *constrainUpperBound = new double[numConstrains];
    cout << "Declare done" << endl;

    CoinPackedMatrix *matrix = new CoinPackedMatrix(false, 0, 0);

    double hap_cov = mGraph->getHaploidDepth();
//    double v_hap_cov = mGraph->getAvgRawCoverage()/2;
    double v_hap_cov = mGraph->getHaploidDepth();
    int v_seg_start = mGraph->getVirusSegStart();
    int maxCopy = 999999;
    double max_cov = -1;
    double max_cp = -1;
    double min_cov = 999999;
    vector<double> covs;
    // constrains for segments, virus and host difference
    for (int i = 0; i < segs->size(); i++) {
        covs.push_back((*segs)[i]->getWeight()->getCorrectedCoverage());
        if ((*segs)[i]->getWeight()->getCorrectedCoverage() > max_cov) {
            // max_cov = (*segs)[i]->getWeight()->getCopyNum();
            max_cov = (*segs)[i]->getWeight()->getCorrectedCoverage();
        }
        if ((*segs)[i]->getWeight()->getCopyNum() > max_cp) {
            // max_cov = (*segs)[i]->getWeight()->getCopyNum();
            max_cp = (*segs)[i]->getWeight()->getCopyNum();
        }
        if ((*segs)[i]->getWeight()->getCorrectedCoverage() < min_cov) {
            // max_cov = (*segs)[i]->getWeight()->getCopyNum();
            min_cov = (*segs)[i]->getWeight()->getCorrectedCoverage();
        }
        // t_i + e_i_1 >= c_i
        CoinPackedVector constrain1;
        // constrain1.insert(i, 1);  // t_i
        // constrain1.insert(numSegsJuncs + juncs->size() + i, 1);  // e_i_1
        // constrainLowerBound[4 * i] = (*segs)[i]->getWeight()->getCopyNum();
        // constrainUpperBound[4 * i] = si->getInfinity();
        if (i < v_seg_start)
            constrain1.insert(i, hap_cov);  // t_i
        else
            constrain1.insert(i, v_hap_cov);
        constrain1.insert(numSegsJuncs + juncs->size() + i, 1);  // e_i_1
        constrainLowerBound[4 * i] = (*segs)[i]->getWeight()->getCorrectedCoverage();
        constrainUpperBound[4 * i] = si->getInfinity();
        matrix->appendRow(constrain1);

        // t_i - e_i_1 <= c_i
        CoinPackedVector constrain2;
        // constrain2.insert(i, 1);
        // constrain2.insert(numSegsJuncs + juncs->size() + i, -1);
        // constrainLowerBound[4 * i + 1] = -1 * si->getInfinity();
        // constrainUpperBound[4 * i + 1] = (*segs)[i]->getWeight()->getCopyNum();
        if (i < v_seg_start)
            constrain2.insert(i, hap_cov);  // t_i
        else
            constrain2.insert(i, v_hap_cov);
        constrain2.insert(numSegsJuncs + juncs->size() + i, -1);
        constrainLowerBound[4 * i + 1] = -1 * si->getInfinity();
        constrainUpperBound[4 * i + 1] = (*segs)[i]->getWeight()->getCorrectedCoverage();
        matrix->appendRow(constrain2);

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
        for (Edge *e : *((*segs)[i]->getPositiveVertex()->getEdgesAsTarget())) {
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
        for (Edge *e : *((*segs)[i]->getPositiveVertex()->getEdgesAsSource())) {
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
        covs.push_back((*juncs)[i]->getWeight()->getCorrectedCoverage());
        if ((*juncs)[i]->getWeight()->getCorrectedCoverage() > max_cov) {
            // max_cov = (*juncs)[i]->getWeight()->getCopyNum();
            max_cov = (*juncs)[i]->getWeight()->getCorrectedCoverage();
        }
        if ((*juncs)[i]->getWeight()->getCopyNum() > max_cp) {
            // max_cov = (*juncs)[i]->getWeight()->getCopyNum();
            max_cp = (*juncs)[i]->getWeight()->getCopyNum();
        }
        if ((*juncs)[i]->getWeight()->getCorrectedCoverage() < min_cov) {
            // max_cov = (*juncs)[i]->getWeight()->getCopyNum();
            min_cov = (*juncs)[i]->getWeight()->getCorrectedCoverage();
        }
        // t_i - c_ix_i + e_i_1 >= 0
        // double copy = (*juncs)[i]->getWeight()->getCopyNum() + 0.01;
        double cov = (*juncs)[i]->getWeight()->getCorrectedCoverage() + 0.05;
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
    auto sources = mGraph->getMSources();
    auto sinks = mGraph->getMSinks();

    // for (int i = 0; i < sourceSize; i++) {
    //     CoinPackedVector constrainSource;
    //     CoinPackedVector constrainSink;
    //     constrainSource.insert((*sources)[i]->getId() - 1, 1);
    //     constrainSource.insert(numVariables - 2 * (i + 1), -1);
    //     // cout << mGraph->getExpectedPloidy() << endl;
    //     constrainLowerBound[numConstrains - 2 * (i + 1)] = mGraph->getExpectedPloidy();
    //     constrainUpperBound[numConstrains - 2 * (i + 1)] = mGraph->getExpectedPloidy();
    //     matrix->appendRow(constrainSource);

    //     constrainSink.insert((*sinks)[i]->getId() - 1, 1);
    //     constrainSink.insert(numVariables - 2*i - 1, -1);
    //     constrainLowerBound[numConstrains - 2*i - 1] = mGraph->getExpectedPloidy();
    //     constrainUpperBound[numConstrains - 2*i - 1] = mGraph->getExpectedPloidy();
    //     matrix->appendRow(constrainSink);
    // }
//    CoinPackedVector constrainSource;
//    constrainSource.insert(mGraph->getFirstSource()->getId() - 1, 1);
//    constrainSource.insert(numVariables - 2, -1);
//    // cout << mGraph->getExpectedPloidy() << endl;
//    constrainLowerBound[numConstrains - 2] = mGraph->getExpectedPloidy();
//    constrainUpperBound[numConstrains - 2] = mGraph->getExpectedPloidy();
//    matrix->appendRow(constrainSource);
//
//    CoinPackedVector constrainSink;
//    constrainSink.insert(mGraph->getFirstSink()->getId() - 1, 1);
//    constrainSink.insert(numVariables - 1, -1);
//    constrainLowerBound[numConstrains - 1] = mGraph->getExpectedPloidy();
//    constrainUpperBound[numConstrains - 1] = mGraph->getExpectedPloidy();
//    matrix->appendRow(constrainSink);
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
            if (i < numVariables) {
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
            cout << "infer info - 0: " << i << " " << segs->size() + i << " " << numSegsJuncs + i << " "
                 << (*juncs)[i]->getInfo()[0] << " " << endl;
        } else {
            variableLowerBound[numSegsJuncs + i] = 1;
            cout << "infer info - 1: " << i << " " << segs->size() + i << " " << numSegsJuncs + i << " "
                 << (*juncs)[i]->getInfo()[0] << " " << endl;
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
    // for (int i = 0; i < sourceSize; i++) {
    //     variableLowerBound[numVariables - 2 * (i + 1)] = 0;
    //     variableUpperBound[numVariables - 2 * (i + 1)] = si->getInfinity();
    //     //  variableUpperBound[numVariables - 2] = 0;
    //     variableLowerBound[numVariables - 2*i - 1] = 0;
    //     variableUpperBound[numVariables - 2*i - 1] = si->getInfinity();
    // }
//    variableLowerBound[numVariables - 2] = 0;
//    variableUpperBound[numVariables - 2] = si->getInfinity();
//    //  variableUpperBound[numVariables - 2] = 0;
//    variableLowerBound[numVariables - 1] = 0;
//    variableUpperBound[numVariables - 1] = si->getInfinity();
    // variableUpperBound[numVariables - 1] = 0;

    si->loadProblem(*matrix, variableLowerBound, variableUpperBound, objective, constrainLowerBound,
                    constrainUpperBound);
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

int LocalGenomicMap::strandCross(Vertex *aVertex) {
    int flag = 0;
    Vertex *startVertex = aVertex;
    VertexPath vertexStack;
    vertexStack.push_back(startVertex);
    Edge *prevEdge = this->selectPrevEdge(startVertex);
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
            Vertex *prevVertex = prevEdge->getSource();
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
    Edge *nextEdge = this->selectNextEdge(startVertex);
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
            Vertex *nextVertex = nextEdge->getTarget();
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

bool LocalGenomicMap::doesPathExists(Vertex *aStartVertex, Vertex *aEndVertex) {
    bool isReach = false;
    Vertex *startVertex = aStartVertex;
    VertexPath vertexStack;
    vertexStack.push_back(startVertex);
    auto partitionPair = findPartition(aStartVertex->getId());
    Edge *nextEdge = this->selectNextEdgeByPartition(partitionPair.first, partitionPair.second, startVertex);
    while (true) {
        if (nextEdge == NULL) {
            vertexStack.pop_back();
            if (vertexStack.size() == 0) {
                break;
            } else {
                startVertex = vertexStack.back();
                nextEdge = this->selectNextEdgeByPartition(partitionPair.first, partitionPair.second, startVertex);
            }
        } else {
            Vertex *nextVertex = nextEdge->getTarget();
            // cout << nextEdge->getInfo() << ", forward" << endl;

            // if (nextVertex == mGraph->getFirstSink()->getNegativeVertex()) {
            //     if (startVertex != mGraph->getFirstSource()->getNegativeVertex()) {
            //         throw ForwardReachSinkNegativeException(aStartVertex);
            //     }
            // }
            // if (nextVertex == mGraph->getFirstSource()->getPositiveVertex()) {
            //     if (startVertex != mGraph->getFirstSink()->getPositiveVertex()) {
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
            nextEdge = this->selectNextEdgeByPartition(partitionPair.first, partitionPair.second, startVertex);
        }
    }
    mGraph->resetJunctionVisitFlag();
    return isReach;
}

double LocalGenomicMap::inferCoverage(Vertex *aSource, Vertex *aTarget) {
    double inCoverage = aTarget->getWeight()->getCoverage() - aTarget->getInCoverage();
    double outCoverage = aSource->getWeight()->getCoverage() - aSource->getOutCoverage();
    return max(1.0, (inCoverage + outCoverage) / 2.0);
}

double LocalGenomicMap::weightedCredibility(Vertex *aVertex, bool aIsSource) {
    // the coefficient is modified by a constant setting
    if (aIsSource) {
        return 1.0 * aVertex->getCredibility() *
               max(0.0, aVertex->getWeight()->getCoverage() - aVertex->getOutCoverage()) / mGraph->getAvgCoverage();
    } else {
        return 1.0 * aVertex->getCredibility() *
               max(0.0, aVertex->getWeight()->getCoverage() - aVertex->getInCoverage()) / mGraph->getAvgCoverage();
    }
}

double LocalGenomicMap::inferCredibility(Vertex *aSource, Vertex *aTarget) {
    return (this->weightedCredibility(aSource, true) + this->weightedCredibility(aTarget, false)) / 2;
}

void LocalGenomicMap::connectSourceSink() {
    int i = 0;
    auto sources = mGraph->getMSources();
    auto sinks = mGraph->getMSinks();
    for (; i < sources->size(); i++) {
        mGraph->addJunction((*sinks)[i]->getPositiveVertex(), (*sources)[i]->getPositiveVertex(),
                            ((*sources)[i]->getWeight()->getCoverage() +
                             (*sinks)[i]->getWeight()->getCoverage()) / 2, 1.0, -1, true, false, true);
    }
}

void LocalGenomicMap::countReachability(int *nReachability, VertexPath &notReachableVertices, Vertex *targetVertex) {
    cout << "count reach" << endl;
    // mGraph->print();
    vector<pair<Vertex *, int> > vertex_reach_count;
    for (VertexPath::const_iterator i = notReachableVertices.begin(); i != notReachableVertices.end(); i++) {
        // cout << "Counting Reachability for " << (*i)->getInfo() << endl;
        int count = 0;
        for (VertexPath::const_iterator j = notReachableVertices.begin(); j != notReachableVertices.end(); j++) {
            if (j == i) continue;
            if (targetVertex == mGraph->getFirstSource()->getPositiveVertex() ||
                targetVertex == mGraph->getFirstSink()->getNegativeVertex()) {
                if (mGraph->BFS(*i, *j) >= 0) {
                    // if (this->doesPathExists(*i, *j)) {
                    count++;
                    // cout << (*i)->getInfo() << "=>" << (*j)->getInfo() << " connectable" << endl;
                } else {
                    // cout << (*i)->getInfo() << "=>" << (*j)->getInfo() << " unconnectable" << endl;
                }
            } else if (targetVertex == mGraph->getFirstSource()->getNegativeVertex() ||
                       targetVertex == mGraph->getFirstSink()->getPositiveVertex()) {
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
    for (vector<pair<Vertex *, int> >::reverse_iterator i = vertex_reach_count.rbegin();
         i != vertex_reach_count.rend(); i++) {
        notReachableVertices.push_back(i->first);
    }
}

void LocalGenomicMap::findWithReachability(int *nReachability, VertexPath &notReachableVertices, int reachability) {
    VertexPath backup = notReachableVertices;
    notReachableVertices.clear();
    for (int i = 0; i < backup.size(); i++) {
        if (nReachability[i] == reachability) {
            notReachableVertices.push_back(backup[i]);
        }
    }
}

void LocalGenomicMap::reconnectNegativeToPositive(JunctionDB *aJuncDB, bool verbose) {
    Vertex *v;
    Record *rec;
    Junction *inferredJunc;
    VertexPath asSourceVertices, asTargetVertices;
    for (Segment *seg: *(mGraph->getSegments())) {
        if (seg->getId() > mGraph->getFirstSink()->getId()) {
            cout << "Continue: " << seg->getId() << endl;
            continue;
        }
        v = seg->getNegativeVertex();
        bool asTargetFromPositive = false, asSourceToPositive = false;
        for (Edge *e: *(v->getEdgesAsSource())) {
            if (e->getTarget()->getDir() == '+') {
                asSourceToPositive = true;
                break;
            }
        }
        for (Edge *e: *(v->getEdgesAsTarget())) {
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
        for (Vertex *v: asSourceVertices) {
            cout << "Current vertex: " << v->getInfo() << " " << v->getSegment()->getChrom() << " " << v->getStart()
                 << v->getDir() << endl;
            cout << "Finding previous vertex" << endl;
            // v = seg->getNegativeVertex();
            Segment *seg = v->getSegment();
            // for (Edge * e: *(v->getEdgesAsSource())) {
            // V- => V+, add junction s.t. V+ => V-
            Vertex *prevVertex;

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
                for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin();
                     riter != rec->getBackwardEntries()->rend(); riter++) {
                    cout << (*riter)->chrom << " " << (*riter)->pos << (*riter)->strand << endl;
                    if ((*riter)->strand == '+') {
                        try {
                            prevVertex = mGraph->getSegmentByChromEnd((*riter)->chrom,
                                                                      (*riter)->pos)->getPositiveVertex();
                            if (verbose)
                                cout << "prev: " << prevVertex->getInfo() << " " << prevVertex->getSegment()->getChrom()
                                     << " " << v->getEnd() << v->getDir() << endl;
                            inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex, v),
                                                               this->inferCredibility(prevVertex, v), -1, true, false,
                                                               false);
                            if (inferredJunc == NULL) continue;
                            if (verbose)
                                cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", "
                                     << inferredJunc->getInfo()[1] << endl;
                            added = true;
                            break;
                        } catch (DuplicateJunctionException &e) {
                            cout << e.what() << endl;
                        } catch (SegmentDoesNotExistException &e) {
                            cout << e.what() << endl;
                        }
                    }
                }
            }
            if (!added) {
                cout << " As source N2P not added" << endl;
                prevVertex = v;
                while (true) {
                    if (prevVertex->getId() >= mGraph->getFirstSink()->getId() || prevVertex->getId() == 1) {
                        // prevVertex = mGraph->getFirstSink()->getPositiveVertex();
                        // prevVertex = mGraph->getSegmentById(mGraph->getFirstSink()->getId() - 1)->getPositiveVertex();
                        prevVertex = mGraph->getSegmentById(mGraph->getFirstSink()->getId() - 1)->getPositiveVertex();
                        // prevVertex = mGraph->getFirstSink()->getPositiveVertex();
                    } else {
                        prevVertex = mGraph->getSegmentById(prevVertex->getId() - 1)->getPositiveVertex();
                    }
                    try {
                        if (verbose)
                            cout << "prev: " << prevVertex->getInfo() << " " << prevVertex->getSegment()->getChrom()
                                 << " " << v->getEnd() << v->getDir() << endl;
                        inferredJunc = mGraph->addJunction(prevVertex, v, this->inferCoverage(prevVertex, v),
                                                           this->inferCredibility(prevVertex, v), -1, true, false,
                                                           false);
                        if (inferredJunc == NULL) {
                            continue;
                        };
                        if (verbose)
                            cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", "
                                 << inferredJunc->getInfo()[1] << endl;
                        break;
                    } catch (DuplicateJunctionException &e) {
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
        for (Vertex *v: asTargetVertices) {
            cout << "Current vertex: " << v->getInfo() << " " << v->getSegment()->getChrom() << " " << v->getEnd()
                 << v->getDir() << endl;
            cout << "Finding next vertex" << endl;
            Segment *seg = v->getSegment();
            // for (Edge * e: *(v->getEdgesAsTarget())) {
            // V- => V+, add junction s.t. V+ => V-
            Vertex *nextVertex;

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
                for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin();
                     riter != rec->getForwardEntries()->rend(); riter++) {
                    if ((*riter)->strand == '+') {
                        try {
                            nextVertex = mGraph->getSegmentByChromStart((*riter)->chrom,
                                                                        (*riter)->pos)->getPositiveVertex();
                            if (verbose)
                                cout << "next: " << nextVertex->getInfo() << " " << nextVertex->getSegment()->getChrom()
                                     << " " << v->getStart() << v->getDir() << endl;
                            // cout << "next: " << nextVertex->getInfo() << endl;
                            inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex),
                                                               this->inferCredibility(v, nextVertex), -1, true, false,
                                                               false);
                            if (inferredJunc == NULL) continue;
                            if (verbose)
                                cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", "
                                     << inferredJunc->getInfo()[1] << endl;
                            added = true;
                            break;
                        } catch (DuplicateJunctionException &e) {
                            cout << e.what() << endl;
                        } catch (SegmentDoesNotExistException &e) {
                            cout << e.what() << endl;
                        }
                    }
                }
            }
            if (!added) {
                cout << " As target N2P not added" << endl;
                nextVertex = v;
                while (true) {
                    if (nextVertex->getId() >= mGraph->getFirstSink()->getId() ||
                        nextVertex->getId() == mGraph->getFirstSink()->getId()) {
                        // nextVertex = mGraph->getFirstSource()->getPositiveVertex();
                        nextVertex = mGraph->getSegmentById(mGraph->getFirstSource()->getId() + 1)->getPositiveVertex();
                    } else {
                        // nextVertex = mGraph->getSegmentById(mGraph->getFirstSource()->getId() + 1)->getPositiveVertex();
                        nextVertex = mGraph->getSegmentById(nextVertex->getId() + 1)->getPositiveVertex();
                    }
                    // if (nextVertex->getId() > mGraph->getFirstSink()->getId()) {
                    //     continue;
                    // } // TODO
                    // if (nextVertex->getId() == mGraph->getFirstSink()->getId()) {
                    //     nextVertex = mGraph->getFirstSource()->getPositiveVertex();
                    // }
                    // nextVertex = mGraph->getSegmentById(nextVertex->getId() + 1)->getPositiveVertex();
                    try {
                        if (verbose)
                            cout << "next: " << nextVertex->getInfo() << " " << nextVertex->getSegment()->getChrom()
                                 << " " << v->getStart() << v->getDir() << endl;
                        inferredJunc = mGraph->addJunction(v, nextVertex, this->inferCoverage(v, nextVertex),
                                                           this->inferCredibility(v, nextVertex), -1, true, false,
                                                           false);
                        if (inferredJunc == NULL) {
                            continue;
                        };
                        if (verbose)
                            cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", "
                                 << inferredJunc->getInfo()[1] << endl;
                        break;
                    } catch (DuplicateJunctionException &e) {
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
        for (Segment *seg: *(mGraph->getSegments())) {
            v = seg->getNegativeVertex();
            bool asTargetFromPositive = false, asSourceToPositive = false;
            for (Edge *e: *(v->getEdgesAsSource())) {
                if (e->getTarget()->getDir() == '+') {
                    asSourceToPositive = true;
                    break;
                }
            }
            for (Edge *e: *(v->getEdgesAsTarget())) {
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

Vertex *LocalGenomicMap::getMostReachable(VertexPath &notReachableVertices, Vertex *targetVertex) {
    Vertex *mostReachable = NULL;
    if (targetVertex == mGraph->getFirstSource()->getPositiveVertex()) {
        // use the least id
        int mostReachableId = mGraph->getSegments()->size();
        for (Vertex *v : notReachableVertices) {
            if (v->getId() < mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    } else if (targetVertex == mGraph->getFirstSink()->getNegativeVertex()) {
        // use the greatest id
        int mostReachableId = -1;
        for (Vertex *v : notReachableVertices) {
            if (v->getId() > mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    } else if (targetVertex == mGraph->getFirstSource()->getNegativeVertex()) {
        // use the least id
        int mostReachableId = mGraph->getSegments()->size();
        for (Vertex *v : notReachableVertices) {
            if (v->getId() < mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    } else if (targetVertex == mGraph->getFirstSink()->getPositiveVertex()) {
        // use the greatest id
        int mostReachableId = -1;
        for (Vertex *v : notReachableVertices) {
            if (v->getId() > mostReachableId) {
                mostReachable = v;
                mostReachableId = v->getId();
            }
        }
    }
    return mostReachable;
}

bool LocalGenomicMap::adjustReachability(VertexPath &notReachableVertices, Vertex *targetVertex, JunctionDB *aJuncDB,
                                         bool verbose) {
    // sort(notReachableVertices.begin(), notReachableVertices.end(), [](Vertex *v1, Vertex *v2) { return v1->getWeight()->getCoverage() > v2->getWeight()->getCoverage(); });
    // for (Vertex * v: notReachableVertices) {
    //     cout << v->getInfo() << ":" << v->getWeight()->getCoverage() << " ";
    // }
    // cout << endl;
    int *nReachability = new int[notReachableVertices.size()];
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

        if (targetVertex == mGraph->getFirstSource()->getPositiveVertex() ||
            targetVertex == mGraph->getFirstSink()->getNegativeVertex()) {
            cout << "Current vertex: " << mostReachable->getId() << " " << mostReachable->getSegment()->getChrom()
                 << " " << mostReachable->getStart() << mostReachable->getDir() << endl;

            rec = aJuncDB->findRecord(mostReachable->getSegment()->getChrom(), mostReachable->getStart(),
                                      mostReachable->getDir());
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
                for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin();
                     riter != rec->getBackwardEntries()->rend(); riter++) {
                    cout << "backward " << (*riter)->chrom << " " << (*riter)->pos << " " << (*riter)->strand << endl;
                    try {
                        if ((*riter)->strand == '+') {
                            prevVertex = mGraph->getSegmentByChromEnd((*riter)->chrom,
                                                                      (*riter)->pos)->getPositiveVertex();
                        } else {
                            prevVertex = mGraph->getSegmentByChromStart((*riter)->chrom,
                                                                        (*riter)->pos)->getNegativeVertex();
                        }
                        if (verbose)
                            cout << "prev: " << prevVertex->getInfo() << endl;
                        inferredJunc = mGraph->addJunction(prevVertex, mostReachable,
                                                           this->inferCoverage(prevVertex, mostReachable),
                                                           this->inferCredibility(prevVertex, mostReachable), -1, true,
                                                           false, false);
                        if (inferredJunc == NULL) continue;
                        if (verbose)
                            cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", "
                                 << inferredJunc->getInfo()[1] << endl;
                        hasAdded = true;
                        break;
                    } catch (SegmentDoesNotExistException &e) {
                        cout << e.what() << endl;
                    } catch (DuplicateJunctionException &e) {
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
                        if (prevVertex->getId() == mGraph->getFirstSource()->getId()) {
                            prevVertex = mGraph->getFirstSink()->getPositiveVertex();
                            // break;
                        } else if (prevVertex->getId() > mGraph->getFirstSink()->getId()) {
                            break;
                        } else {
                            prevVertex = mGraph->getSegmentById(prevVertex->getId() - 1)->getPositiveVertex();
                        }
                    } else {
                        // prevVertex = mGraph->getSegmentById(mostReachable->getId() + 1)->getNegativeVertex();
                        if (prevVertex->getId() == mGraph->getFirstSink()->getId()) {
                            prevVertex = mGraph->getFirstSource()->getNegativeVertex();
                            // break;
                        } else if (prevVertex->getId() > mGraph->getFirstSink()->getId()) {
                            break;
                        } else {
                            prevVertex = mGraph->getSegmentById(prevVertex->getId() + 1)->getNegativeVertex();
                        }
                    }
                    try {
                        if (verbose)
                            cout << "prev: " << prevVertex->getInfo() << endl;
                        inferredJunc = mGraph->addJunction(prevVertex, mostReachable,
                                                           this->inferCoverage(prevVertex, mostReachable),
                                                           this->inferCredibility(prevVertex, mostReachable), -1, true,
                                                           false, false);
                        if (inferredJunc == NULL) continue;
                        if (verbose)
                            cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", "
                                 << inferredJunc->getInfo()[1] << endl;
                        hasAdded = true;
                        break;
                    } catch (DuplicateJunctionException &e) {
                        cout << e.what() << endl;
                    }
                }
            }
        } else {
            rec = aJuncDB->findRecord(mostReachable->getSegment()->getChrom(), mostReachable->getEnd(),
                                      mostReachable->getDir());
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
                for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin();
                     riter != rec->getForwardEntries()->rend(); riter++) {
                    cout << "forward " << (*riter)->chrom << " " << (*riter)->pos << endl;
                    try {
                        if ((*riter)->strand == '+') {
                            nextVertex = mGraph->getSegmentByChromStart((*riter)->chrom,
                                                                        (*riter)->pos)->getPositiveVertex();
                        } else {
                            nextVertex = mGraph->getSegmentByChromEnd((*riter)->chrom,
                                                                      (*riter)->pos)->getNegativeVertex();
                        }
                        if (verbose)
                            cout << "prev: " << nextVertex->getInfo() << endl;
                        inferredJunc = mGraph->addJunction(mostReachable, nextVertex,
                                                           this->inferCoverage(mostReachable, nextVertex),
                                                           this->inferCredibility(mostReachable, nextVertex), -1, true,
                                                           false, false);
                        if (inferredJunc == NULL) continue;
                        if (verbose)
                            cout << "add junction with juncdb: " << inferredJunc->getInfo()[0] << ", "
                                 << inferredJunc->getInfo()[1] << endl;
                        hasAdded = true;
                        break;
                    } catch (SegmentDoesNotExistException &e) {
                        cout << e.what() << endl;
                    } catch (DuplicateJunctionException &e) {
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
                        if (nextVertex->getId() == mGraph->getFirstSink()->getId()) {
                            nextVertex = mGraph->getFirstSource()->getPositiveVertex();
                            // break;
                        } else if (nextVertex->getId() > mGraph->getFirstSink()->getId()) {
                            break;
                        } else {
                            nextVertex = mGraph->getSegmentById(nextVertex->getId() + 1)->getPositiveVertex();
                        }
                    } else {
                        // nextVertex = mGraph->getSegmentById(mostReachable->getId() - 1)->getNegativeVertex();
                        if (nextVertex->getId() == mGraph->getFirstSource()->getId()) {
                            nextVertex = mGraph->getFirstSink()->getNegativeVertex();
                            // break;
                        } else if (nextVertex->getId() > mGraph->getFirstSink()->getId()) {
                            break;
                        } else {
                            nextVertex = mGraph->getSegmentById(nextVertex->getId() - 1)->getNegativeVertex();
                        }
                    }
                    try {
                        if (verbose)
                            cout << "next: " << nextVertex->getInfo() << endl;
                        inferredJunc = mGraph->addJunction(mostReachable, nextVertex,
                                                           this->inferCoverage(mostReachable, nextVertex),
                                                           this->inferCredibility(mostReachable, nextVertex), -1, true,
                                                           false, false);
                        if (inferredJunc == NULL) continue;
                        if (verbose)
                            cout << "add junction without juncdb: " << inferredJunc->getInfo()[0] << ", "
                                 << inferredJunc->getInfo()[1] << endl;
                        hasAdded = true;
                        break;
                    } catch (DuplicateJunctionException &e) {
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
    // if (targetVertex == mGraph->getFirstSource()->getPositiveVertex()) {
    //     // backward to source+
    //     Vertex * prevVertex;
    //     while (true) {
    //         try {
    //             Vertex * selectedVertex = this->selectPrevVertex(mostReachable, mGraph->getFirstSource()->getPositiveVertex(), aJuncDB);
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
    // } else if (targetVertex == mGraph->getFirstSink()->getNegativeVertex()) {
    //     // backward to sink-
    //     Vertex * prevVertex;
    //     while (true) {
    //         try {
    //             Vertex * selectedVertex = this->selectPrevVertex(mostReachable, mGraph->getFirstSink()->getNegativeVertex(), aJuncDB);
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
    // } else if (targetVertex == mGraph->getFirstSink()->getPositiveVertex()) {
    //     // forward to sink+
    //     Vertex * nextVertex;
    //     while (true) {
    //         try {
    //             Vertex * selectedVertex = this->selectNextVertex(mostReachable, mGraph->getFirstSink()->getPositiveVertex(), aJuncDB);
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
    // } else if (targetVertex == mGraph->getFirstSource()->getNegativeVertex()) {
    //     // forward to source-
    //     Vertex * nextVertex;
    //     while (true) {
    //         try {
    //             Vertex * selectedVertex = this->selectNextVertex(mostReachable, mGraph->getFirstSource()->getNegativeVertex(), aJuncDB);
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
    for (Segment *seg: *(mGraph->getSegments())) {
        if (seg == mGraph->getFirstSource() || seg == mGraph->getFirstSink()) {
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
        for (Vertex *v : segV) {
            // cout << "source+ -> " << v->getInfo() << endl;
            bool isBackwardSourceReachable = this->doesPathExists(mGraph->getFirstSource()->getPositiveVertex(), v);
            // cout << "sink- -> " << v->getInfo() << endl;
            bool isBackwardSinkReachable = this->doesPathExists(mGraph->getFirstSink()->getNegativeVertex(), v);
            // cout << v->getInfo() << " -> source-" << endl;
            bool isForwardSourceReachable = this->doesPathExists(v, mGraph->getFirstSource()->getNegativeVertex());
            // cout << v->getInfo() << " -> sink+" << endl;
            bool isForwardSinkReachable = this->doesPathExists(v, mGraph->getFirstSink()->getPositiveVertex());
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

    return backwardSourceNotReachableVertices.size() + backwardSinkNotReachableVertices.size() +
           forwardSourceNotReachableVertices.size() + forwardSinkNotReachableVertices.size() == 0;
}

void LocalGenomicMap::checkReachability(JunctionDB *aJuncDB, bool verbose) {
    cout << "Checking reachability..." << endl;
    VertexPath backwardSourceNotReachableVertices;
    VertexPath backwardSinkNotReachableVertices;
    VertexPath forwardSourceNotReachableVertices;
    VertexPath forwardSinkNotReachableVertices;
    int i = 0;
    do {
        backwardSourceNotReachableVertices.clear();
        backwardSinkNotReachableVertices.clear();
        forwardSourceNotReachableVertices.clear();
        forwardSinkNotReachableVertices.clear();

        mGraph->checkOrphan();

        // check reachability for each vertex
        bool isSinkSourceConnectedOriginally;
        try {
            this->connectSourceSink();  // connect all the sinks and the corresponding sources
//            isSinkSourceConnectedOriginally = false;
        } catch (DuplicateJunctionException &e) {
            isSinkSourceConnectedOriginally = true;
        }
        for (Segment *seg: *(mGraph->getSegments())) {
            if (std::count(mGraph->getMSources()->begin(), mGraph->getMSources()->end(), seg)
                || std::count(mGraph->getMSinks()->begin(), mGraph->getMSinks()->end(), seg)) {
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
            auto segPartitionPair = findPartition(seg->getId());
//            ??????seg?????????seg??????????????????????????????reachable?????????, ????????????????????????partition???reachable
            if (false) {
            //if (segPartitionPair.first == this->mGraph->getMSources()->back()->getId()) {
                auto sources = this->mGraph->getMSources();
                auto sinks = this->mGraph->getMSinks();
                for (Vertex *v : segV) {
                    bool finalReachable = false;
                    for (int i = 0; i < sources->size(); i++) {
                        auto tSource = (*sources)[i];
                        auto tSinks = (*sinks)[i];
                        // cout << "source+ -> " << v->getInfo() << endl;
                        bool isBackwardSourceReachable = this->doesPathExists(tSource->getPositiveVertex(), v);
                        // cout << "sink- -> " << v->getInfo() << endl;
                        bool isBackwardSinkReachable = this->doesPathExists(tSinks->getNegativeVertex(), v);
                        // cout << v->getInfo() << " -> source-" << endl;
                        bool isForwardSourceReachable = this->doesPathExists(v, tSource->getNegativeVertex());
                        // cout << v->getInfo() << " -> sink+" << endl;
                        bool isForwardSinkReachable = this->doesPathExists(v, tSinks->getPositiveVertex());
                        if (vReachable(isBackwardSourceReachable, isForwardSinkReachable, isBackwardSinkReachable,
                                       isForwardSourceReachable)) {
                            finalReachable = true;
                            break;
                        }
                    }
                    if(!finalReachable) {
                        cout << v->getInfo()<<"-----\n";
                    }
                }
            } else {
                auto segPartitionSource = this->mGraph->getSegmentById(segPartitionPair.first);
                auto segPartitionSink = this->mGraph->getSegmentById(segPartitionPair.second);
//                auto segPartitionSource = this->mGraph->getFirstSource();
//                auto segPartitionSink = this->mGraph->getFirstSink();
                for (Vertex *v : segV) {
                    // cout << "source+ -> " << v->getInfo() << endl;
                    bool isBackwardSourceReachable = this->doesPathExists(segPartitionSource->getPositiveVertex(), v);
                    // cout << "sink- -> " << v->getInfo() << endl;
                    bool isBackwardSinkReachable = this->doesPathExists(segPartitionSink->getNegativeVertex(), v);
                    // cout << v->getInfo() << " -> source-" << endl;
                    bool isForwardSourceReachable = this->doesPathExists(v, segPartitionSource->getNegativeVertex());
                    // cout << v->getInfo() << " -> sink+" << endl;
                    bool isForwardSinkReachable = this->doesPathExists(v, segPartitionSink->getPositiveVertex());
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
//        if (!hasAdded && backwardSourceNotReachableVertices.size() > 0) {
//            if (verbose) cout << "for backward source: " << endl;
//            hasAdded = this->adjustReachability(backwardSourceNotReachableVertices,
//                                                mGraph->getFirstSource()->getPositiveVertex(), aJuncDB, verbose);
//            // continue;
//        }
//        if (!hasAdded && backwardSinkNotReachableVertices.size() > 0) {
//            if (verbose) cout << "for backward sink: " << endl;
//            hasAdded = this->adjustReachability(backwardSinkNotReachableVertices,
//                                                mGraph->getFirstSink()->getNegativeVertex(), aJuncDB, verbose);
//            // continue;
//        }
//        if (!hasAdded && forwardSourceNotReachableVertices.size() > 0) {
//            if (verbose) cout << "for forward source: " << endl;
//            hasAdded = this->adjustReachability(forwardSourceNotReachableVertices,
//                                                mGraph->getFirstSource()->getNegativeVertex(), aJuncDB, verbose);
//            // continue;
//        }
//        if (!hasAdded && forwardSinkNotReachableVertices.size() > 0) {
//            if (verbose) cout << "for forward sink: " << endl;
//            hasAdded = this->adjustReachability(forwardSinkNotReachableVertices,
//                                                mGraph->getFirstSink()->getPositiveVertex(),
//                                                aJuncDB, verbose);
//        }
	i = i+1;
        // break;
    } while (backwardSourceNotReachableVertices.size() + backwardSinkNotReachableVertices.size() +
             forwardSourceNotReachableVertices.size() + forwardSinkNotReachableVertices.size() > 0 && i <=10);

    cout << "Graph reachable." << endl;
}

vector<Segment *> *LocalGenomicMap::extractReachableGraph(bool verbose) {
    vector<Segment *> *res = new vector<Segment *>();
    cout << "Checking reachability..." << endl;
    VertexPath backwardSourceNotReachableVertices;
    VertexPath backwardSinkNotReachableVertices;
    VertexPath forwardSourceNotReachableVertices;
    VertexPath forwardSinkNotReachableVertices;
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
    for (Segment *seg: *(mGraph->getSegments())) {
        if (seg == mGraph->getFirstSource() || seg == mGraph->getFirstSink()) {
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
        bool flag = true;
        for (Vertex *v : segV) {
            // cout << "source+ -> " << v->getInfo() << endl;
            bool isBackwardSourceReachable = this->doesPathExists(mGraph->getFirstSource()->getPositiveVertex(), v);
            // cout << "sink- -> " << v->getInfo() << endl;
            bool isBackwardSinkReachable = this->doesPathExists(mGraph->getFirstSink()->getNegativeVertex(), v);
            // cout << v->getInfo() << " -> source-" << endl;
            bool isForwardSourceReachable = this->doesPathExists(v, mGraph->getFirstSource()->getNegativeVertex());
            // cout << v->getInfo() << " -> sink+" << endl;
            bool isForwardSinkReachable = this->doesPathExists(v, mGraph->getFirstSink()->getPositiveVertex());
            if (isBackwardSinkReachable && isBackwardSourceReachable && isForwardSinkReachable &&
                isForwardSourceReachable) {
            } else {
                flag = false;
            }
        }
        if (flag) {
            res->push_back(seg);
        }
    }
    return res;
}

void LocalGenomicMap::addNormalJunctions() {
    for (int i = 0; i < mGraph->getSegments()->size() - 1; i++) {
        Vertex *sourceVertex = (*mGraph->getSegments())[i]->getPositiveVertex();
        Vertex *targetVertex = mGraph->getNextVertexById(sourceVertex);
        if (sourceVertex->getSegment()->getChrom() != targetVertex->getSegment()->getChrom()) {
            continue;
        }
        try {
            mGraph->addJunction(sourceVertex, targetVertex, this->inferCoverage(sourceVertex, targetVertex),
                                inferCredibility(sourceVertex, targetVertex), -1, true, false, false);
        } catch (DuplicateJunctionException &e) {
            continue;
        }
    }
}

void LocalGenomicMap::checkInferredJunctionCredibility() {
    double maxCred = 1.0;
    cout << "Read" << endl;
    for (Junction *junc : *(mGraph->getJunctions())) {
        if (junc->isInferred()) {
            maxCred = max(pow(junc->getCredibility(), 0.5), maxCred);
        }
    }

    cout << "Assign" << endl;
    for (Junction *junc : *(mGraph->getJunctions())) {
        if (junc->isInferred()) {
            cout << junc->getInfo()[0] << endl;
            junc->setCredibility(pow(junc->getCredibility(), 0.5) / maxCred);
        }
    }
    cout << "Assign done" << endl;
}

void LocalGenomicMap::clearSegmentJunctionCredibility(Segment *aSegment) {
    cout << "Clear credibility of segment " << aSegment->getId() << " and its related junctions" << endl;
    aSegment->setCredibility(0);
    aSegment->resetHasLowerBoundLimit();
    for (Edge *e : *(aSegment->getPositiveVertex()->getEdgesAsSource())) {
        Junction *junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }

    for (Edge *e : *(aSegment->getPositiveVertex()->getEdgesAsTarget())) {
        Junction *junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }

    for (Edge *e : *(aSegment->getNegativeVertex()->getEdgesAsSource())) {
        Junction *junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }

    for (Edge *e : *(aSegment->getNegativeVertex()->getEdgesAsTarget())) {
        Junction *junc = e->getJunction();
        junc->setCredibility(0);
        junc->resetHasLowerBoundLimit();
        e->setCredibility(0);
    }
}

Vertex *LocalGenomicMap::selectPrevVertex(Vertex *currentVertex, Vertex *targetVertex, JunctionDB *aJuncDB) {
    Record *rec;
    rec = aJuncDB->findRecord(currentVertex->getSegment()->getChrom(), currentVertex->getStart(),
                              currentVertex->getDir());
    Vertex *selectedVertex = NULL;
    for (vector<entry_t *>::reverse_iterator riter = rec->getBackwardEntries()->rbegin();
         riter != rec->getBackwardEntries()->rend(); riter++) {
        Segment *seg = mGraph->getSegmentByChromEnd((*riter)->chrom, (*riter)->pos);
        if (targetVertex == mGraph->getFirstSource()->getPositiveVertex()) {
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
        } else if (targetVertex == mGraph->getFirstSink()->getNegativeVertex()) {
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

Vertex *LocalGenomicMap::selectNextVertex(Vertex *currentVertex, Vertex *targetVertex, JunctionDB *aJuncDB) {
    Record *rec;
    rec = aJuncDB->findRecord(currentVertex->getSegment()->getChrom(), currentVertex->getEnd(),
                              currentVertex->getDir());
    Vertex *selectedVertex = NULL;
    for (vector<entry_t *>::reverse_iterator riter = rec->getForwardEntries()->rbegin();
         riter != rec->getForwardEntries()->rend(); riter++) {
        Segment *seg = mGraph->getSegmentByChromStart((*riter)->chrom, (*riter)->pos);
        if (targetVertex == mGraph->getFirstSink()->getPositiveVertex()) {
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
        } else if (targetVertex == mGraph->getFirstSource()->getNegativeVertex()) {
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

Edge *LocalGenomicMap::selectPrevEdge(Vertex *aTargetVertex, bool isTraversing) {
    for (Edge *e : *(aTargetVertex->getEdgesAsTarget())) {
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

Edge *LocalGenomicMap::selectNextEdgeByPartition(int partitionStart, int partitionEnd, Vertex *aSourceVertex,
                                                 bool isTraversing) {
    for (Edge *e : *(aSourceVertex->getEdgesAsSource())) {
        int eTargetId = e->getTarget()->getId();
        auto eTargetPartitionPair = findPartition(eTargetId);
        int lastPartitionId = this->getGraph()->getMSinks()->back()->getId();
        if (eTargetPartitionPair.first != 0 &&
            (eTargetPartitionPair.first != partitionStart && eTargetPartitionPair.second != partitionEnd))
            continue;
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
    }
    return NULL;
}

Edge *LocalGenomicMap::selectNextEdge(Vertex *aSourceVertex, bool isTraversing) {
    for (Edge *e : *(aSourceVertex->getEdgesAsSource())) {
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

int LocalGenomicMap::findCircuit(Vertex *aVertex, VertexPath &pathVertices, EdgePath &pathEdges) {
    cout << "Search start: " << aVertex->getInfo() << endl;
    Vertex *currentVertex = aVertex;
    Edge *currentEdge = this->selectNextEdge(currentVertex, true);

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
            Vertex *nextVertex = currentEdge->getTarget();

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

Edge *LocalGenomicMap::traverseNextEdge(Vertex *aStartVertex, VertexPath *vp, JunctionDB *aJuncDB) {
    Edge *selectedEdge = nullptr;
    if (this->isUsingHic()) {
        selectedEdge = traverseWithHic(vp);
        if (selectedEdge != nullptr) return selectedEdge;
    }
    Record *rec = aJuncDB->findRecord(aStartVertex->getSegment()->getChrom(), aStartVertex->getEnd(),
                                      aStartVertex->getDir());
    if (rec == NULL) {
        for (Edge *e: *(aStartVertex->getEdgesAsSource())) {
            if (e->hasCopy()) {
                selectedEdge = e;
                break;
            }
        }
    } else {
        int support = 0;
        entry_t *entry;
        for (Edge *e: *(aStartVertex->getEdgesAsSource())) {
            if (!e->hasCopy()) continue;
            entry = rec->findForwardEntry(e->getTarget()->getSegment()->getChrom(), e->getTarget()->getStart(),
                                          e->getTarget()->getDir());
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

Edge *LocalGenomicMap::traverseNextEdgeByPartition(Vertex *aStartVertex, VertexPath *vp, JunctionDB *aJuncDB,
                                                   int *partitionStart, int *partitionEnd) {
    Edge *selectedEdge = nullptr;
    if (this->isUsingHic()) {
        selectedEdge = traverseWithHic(vp);
        if (selectedEdge != nullptr) return selectedEdge;
    }
//    select jun in juncdb first
    auto *recs = aJuncDB->findRecords(aStartVertex->getSegment()->getChrom(), aStartVertex->getEnd(),
                                      aStartVertex->getDir());
    if (recs != nullptr) {
        int support = 0;
        entry_t *entry;
        for (auto rec : *recs) {
//             for (Edge *e: *(aStartVertex->getEdgesAsSource())) {
//                 if (e->hasCopy()) {
//                     int eTargetId = e->getTarget()->getId();
// //                    TODO ??????check????????????
//                     auto isInPartition = this->checkCommon(eTargetId, partitionStart, partitionEnd);
//                     if (isInPartition) {
//                         entry = rec->findForwardEntry(e->getTarget()->getSegment()->getChrom(),
//                                                       e->getTarget()->getStart(),
//                                                       e->getTarget()->getDir());
//                         if (entry != nullptr) {
//                             if (entry->support > support) {
//                                 support = entry->support;
//                                 selectedEdge = e;
//                                 return selectedEdge;
//                             }
//                         } else {
//                             if (support == 0) {
//                                 selectedEdge = e;
//                                 return selectedEdge;
//                             }
//                         }
//                     }
//                 }
//             }
            for (Edge *e: *(aStartVertex->getEdgesAsSource())) {
                if (e->hasCopy()) {
                    int eTargetId = e->getTarget()->getId();
//                    TODO ??????check????????????
                    auto isInPartition = this->checkPartition(eTargetId, partitionStart, partitionEnd);
                    if (isInPartition) {
                        entry = rec->findForwardEntry(e->getTarget()->getSegment()->getChrom(),
                                                      e->getTarget()->getStart(),
                                                      e->getTarget()->getDir());
                        if (entry != nullptr) {
                            if (entry->support > support) {
                                support = entry->support;
                                selectedEdge = e;
                                return selectedEdge;
                            }
                        } else {
                            if (support == 0) {
                                selectedEdge = e;
                                return selectedEdge;
                            }
                        }
                    }
                }
            }
        }
    }
// if not return, random select e in partition
    // for (Edge *e: *(aStartVertex->getEdgesAsSource())) {
    //     int eTargetId = e->getTarget()->getId();
    //     if (e->hasCopy()) {
    //         auto isInPartition = this->checkCommon(eTargetId, partitionStart, partitionEnd);
    //         if (isInPartition) {
    //             selectedEdge = e;
    //             return selectedEdge;
    //         }
    //     }
    // }

    for (Edge *e: *(aStartVertex->getEdgesAsSource())) {
        int eTargetId = e->getTarget()->getId();
        if (e->hasCopy()) {
            auto isInPartition = this->checkPartition(eTargetId, partitionStart, partitionEnd);
            if (isInPartition) {
                selectedEdge = e;
                return selectedEdge;
            }
        }
    }
    return selectedEdge;
}

//TODO HIC partition
Edge *LocalGenomicMap::traverseWithHic(VertexPath *vp) {
    auto currentVertex = vp->back();
    Edge *maxEdge;
    double maxV = 0;
    for (auto *e : *(currentVertex->getEdgesAsSource())) {
        if (!e->hasCopy()) continue;
        double v = calculateHicInteraction(vp, e->getTarget());
        if (v > maxV) {
            maxEdge = e;
            maxV = v;
        }
    }
    if (maxV == 0) return nullptr;
    this->decreaseHicMatrix(vp, maxEdge);
    return maxEdge;
}

pair<int, int> LocalGenomicMap::findPartition(int id) {
//    if(id == 7) {
//        int m = 0;
//    }
    auto sources = mGraph->getMSources();
    auto sinks = mGraph->getMSinks();
    for (int i = 0; i < sources->size(); i++) {
        auto sId = (*sources)[i]->getId();
        auto eId = (*sinks)[i]->getId();
        if (id >= sId && id <= eId) {
//            if (sId == sources->back()->getId()) {
//                return make_pair(0,0);
//            } else
                return make_pair(sId, eId);
        }
    }
    return make_pair(0, 0);
}

pair<int, int> LocalGenomicMap::findLongPathPartition(VertexPath* vp) {
//    pair<int, int> parPair = findPartition(vp->front()->getId());
//    for(auto v : *vp) {
//        auto tParPair = findPartition(v->getId());
//        if(tParPair.first == 0 || tParPair.first = parPair.first)
//    }
    int startId = mGraph->getMSources()->back()->getId();
    int endId = mGraph->getMSinks()->back()->getId();
    auto v = vp->begin();
    while(v != vp->end()){
        auto tParPair = findPartition((*v)->getId());
        v++;
        if(tParPair.first == mGraph->getMSources()->back()->getId()) continue;
        else if(startId == mGraph->getMSources()->back()->getId()){
            startId = tParPair.first;
            endId = tParPair.second;
        } else if (tParPair.first != startId) {
            return make_pair(-1,-1);
        }
    }
    return make_pair(startId, endId);
}

double LocalGenomicMap::calculateHicInteraction(VertexPath *vp, Vertex *currentVertex) {
    double res = 0;
    for (auto *v : *vp) {
        int id1 = v->getId();
        int id2 = currentVertex->getId();
        cout << this->hicMatrix[1][1] << endl;
        double hicV = this->hicMatrix[id1 - 1][id2 - 1];
        res += hicV;
    }
    return res;
}

void LocalGenomicMap::traverse(Vertex *aStartVertex, JunctionDB *aJuncDB) {
//    current traverse path partition
    int *partitionStart = new int();
    int *partitionEnd = new int();
    *partitionStart = mGraph->getMSources()->back()->getId();
    *partitionEnd = mGraph->getMSinks()->back()->getId();
    auto *vp = new VertexPath();
    auto *ep = new EdgePath();
    Vertex *currentVertex;
//    if (!usingLong) {
//        currentVertex = aStartVertex;
//        vp->push_back(currentVertex);
//    } else {
//        currentVertex = traverseLongPath(aStartVertex, vp, false);
//    }
    bool start = true;
    currentVertex = aStartVertex;
    if (isUsingLong()) {
        checkPartition(currentVertex->getId(), partitionStart, partitionEnd);
        while (true) {
            currentVertex = traverseLongPath(currentVertex, vp, partitionStart, partitionEnd);
            Edge *nextEdge = this->traverseNextEdgeByPartition(currentVertex, vp, aJuncDB, partitionStart, partitionEnd);
            if (nextEdge == nullptr) break;
            nextEdge->traverse();
            currentVertex = nextEdge->getTarget();
            ep->push_back(nextEdge);
        }
    } else {
        vp->push_back(currentVertex);
        checkPartition(currentVertex->getId(), partitionStart, partitionEnd);
        while (true) {
            auto nextEdge = this->traverseNextEdgeByPartition(currentVertex, vp, aJuncDB, partitionStart, partitionEnd);
            currentVertex->traverse();
            if (nextEdge == nullptr) break;
            nextEdge->traverse();
            ep->push_back(nextEdge);
            vp->push_back(nextEdge->getTarget());
            currentVertex = nextEdge->getTarget();
        }
    }
    if (traversedCircuits->count(*partitionStart) != 0)
        (*traversedCircuits)[*partitionStart]->push_back(vp);
    else {
        (*traversedCircuits)[*partitionStart] = new vector<VertexPath *>();
        (*traversedCircuits)[*partitionStart]->push_back(vp);
    }
    mCircuits->push_back(vp);
    cout << "Traversed path: ";
    for (Vertex *v: *vp) {
        cout << v->getInfo() << " ";
    }
    cout << endl;
}

void LocalGenomicMap::decreaseHicMatrix(VertexPath *vp, Edge *e) {
    this->decreaseHicInteraction(e->getSource(), e->getTarget());
    for (auto *v : *vp) {
        this->decreaseHicInteraction(v, e->getTarget());
    }
}

void LocalGenomicMap::decreaseHicInteraction(Vertex *v1, Vertex *v2) {
    int id1 = v1->getId();
    int id2 = v2->getId();
    double v = this->hicMatrix[id1][id2];
    v = v - this->decreaseMatrix[id1][id2];
    this->hicMatrix[id1][id2] = v;
    this->hicMatrix[id2][id1] = v;
}

void LocalGenomicMap::traverseGraphByPartition(JunctionDB *aJuncDB) {

}

void LocalGenomicMap::traverseGraph(JunctionDB *aJuncDB) {
    Vertex *currentVertex;
    // traverse starts from sources
    auto sources = mGraph->getMSources();
    vector<Segment*> segments;
    for (Segment *seg : *(mGraph->getSegments())) {
        //find all the segments except the sources
        if(find(sources->begin(), sources->end(), seg) == sources->end()){
            segments.push_back(seg);
        }
    }
    while (!mGraph->isCopyExhaustive()) {
        for (Segment *src: *sources) {
            if (src->hasCopy()) {
                currentVertex = src->getPositiveVertex();
                this->traverse(currentVertex, aJuncDB);
            }
        }
        for (Segment *seg : segments) {
            if (seg->hasCopy()) {
                currentVertex = seg->getPositiveVertex();
                this->traverse(currentVertex, aJuncDB);
            }
        }
        // mGraph->print();
    }
}

// No partition check;
Vertex *LocalGenomicMap::traverseLongPath(Vertex *aStartVertex, VertexPath *vPath, int *partitionStart, int *partitionEnd) {
    int pathN = -1;
    int maxL = 0;
    int i = 0;
    auto selectedLongPartition = (*this->mLongFrags)[*partitionStart];
//    find longest path that current graph can cover, choose this graph to traverse.
    for (auto path: *selectedLongPartition) {
        if ((*path)[0] == aStartVertex) {
            int l = longPathLenInGraph(path);
            if (l > maxL) {
                maxL = l;
                pathN = i;
            }
        }
        i++;
    }
    i = 0;
    if (maxL <= 1) {
            vPath->push_back(aStartVertex);
            aStartVertex->traverse();
        return aStartVertex;
    }
//    auto m = this->mLongFrags[pathN][0];
    i = 0;
    for (; i < maxL; i++) {
        auto v = (*((*selectedLongPartition)[pathN]))[i];
        vPath->push_back(v);
        v->traverse();
        auto vNext = (*((*selectedLongPartition)[pathN]))[i + 1];
        for (Edge *e: *(v->getEdgesAsSource())) {
            if (e->getTarget() == vNext) {
                e->traverse();
                break;
            }
        }
    }
    return (*((*selectedLongPartition)[pathN]))[maxL - 1];
}

int LocalGenomicMap::longPathLenInGraph(VertexPath *longPath) {
    int n = 1;
    for (auto v = longPath->begin(); v != longPath->end() - 1; v++) {
        bool flag = false;
//        auto t = *(v+1);
        for (Edge *e : *(*v)->getEdgesAsSource()) {
            if (e->getTarget() == *(v + 1) and e->hasCopy()) {
                flag = true;
                break;
            }
        }
        if (flag) n++;
        else {
            n = 1;
            break;
        }
    }
    return n;
}

bool LocalGenomicMap::checkPartition(int eTargetId, int *partitionStart, int *partitionEnd) {
//    return true;
    int lastPartitionId = mGraph->getMSources()->back()->getId();
    if (eTargetId >= lastPartitionId || (eTargetId >= *partitionStart && eTargetId <= *partitionEnd)) {
        return true;
    } else {
        auto partitionPair = findPartition(eTargetId);
        if (*partitionStart == mGraph->getMSources()->back()->getId()) {
            *partitionStart = partitionPair.first;
            *partitionEnd = partitionPair.second;
            return true;
        } else if (*partitionStart == partitionPair.first && *partitionEnd == partitionPair.second) {
            return true;
        } else return false;
    }
}

bool LocalGenomicMap::checkCommon(int eTargetId, int *partitionStart, int *partitionEnd) {
    int lastPartitionId = mGraph->getMSinks()->back()->getId();
    if (eTargetId >= *partitionStart && eTargetId <= *partitionEnd) {
        return true;
    } else {
        auto partitionPair = findPartition(eTargetId);
        if (*partitionStart == mGraph->getMSources()->back()->getId()) {
            *partitionStart = partitionPair.first;
            *partitionEnd = partitionPair.second;
            return true;
        } else if (*partitionStart == partitionPair.first && *partitionEnd == partitionPair.second) {
            return true;
        } else return false;
    }
}

bool
LocalGenomicMap::vReachable(bool isBackwardSourceReachable, bool isForwardSinkReachable, bool isBackwardSinkReachable,
                            bool isForwardSourceReachable) {
    if (!isBackwardSourceReachable && !isForwardSinkReachable
        && !isBackwardSinkReachable && !isForwardSourceReachable) {
        return false;
    }
    if ((isBackwardSourceReachable ^ isForwardSinkReachable) != 0) {
        return false;
    }
    if ((isBackwardSinkReachable ^ isForwardSourceReachable) != 0) {
        return false;
    }
    return true;
}
// void LocalGenomicMap::traverseGraph() {
//     while (!mGraph->isCopyExhaustive()) {
//         for (Segment * seg : *(mGraph->getSegments())) {
//             if (!seg->hasCopy()) continue;
//             // this->traverse(seg->getPositiveVertex());
//             VertexPath * pathV = new VertexPath();
//             EdgePath pathE;
//             if (seg == mGraph->getFirstSource()) {
//                 this->BFS(mGraph->getFirstSource()->getPositiveVertex());
//                 int flag = this->findShortestPath(mGraph->getFirstSource()->getPositiveVertex(), mGraph->getFirstSink()->getPositiveVertex(), *pathV, pathE);
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
//             if (seg == mGraph->getFirstSource()) {
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

void LocalGenomicMap::isCircuitSimple(VertexPath *circuit, pair<int, int> &notSimpleIdx) {
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

void LocalGenomicMap::allCircuitsSimple(vector<tuple<int, int, int> > &notSimpleIdx) {
    for (int i = 0; i < mCircuits->size(); i++) {
        pair<int, int> notSimpleVertexIdx;
        this->isCircuitSimple((*mCircuits)[i], notSimpleVertexIdx);
        if (notSimpleVertexIdx.first >= 0) {
            notSimpleIdx.push_back(make_tuple(i, notSimpleVertexIdx.first, notSimpleVertexIdx.second));
        }
    }
}

void LocalGenomicMap::extractCircuits() {
    vector<tuple<int, int, int> > notSimpleCircuitIdx;
    this->allCircuitsSimple(notSimpleCircuitIdx);
    while (notSimpleCircuitIdx.size() > 0) {
        for (tuple<int, int, int> idxPair : notSimpleCircuitIdx) {
            int circuitIdx, begin, end;
            tie(circuitIdx, begin, end) = idxPair;
            VertexPath *subCircuit = new VertexPath((*mCircuits)[circuitIdx]->begin() + begin,
                                                    (*mCircuits)[circuitIdx]->begin() + end + 1);
            mCircuits->push_back(subCircuit);
            (*mCircuits)[circuitIdx]->erase((*mCircuits)[circuitIdx]->begin() + begin + 1,
                                            (*mCircuits)[circuitIdx]->begin() + end + 1);
        }
        notSimpleCircuitIdx.clear();
        this->allCircuitsSimple(notSimpleCircuitIdx);
    }
}

void LocalGenomicMap::sortCircuits() {
    sort(mCircuits->begin(), mCircuits->end(),
         [](VertexPath *c1, VertexPath *c2) { return (*c1)[0]->getId() < (*c2)[0]->getId(); });
}

void LocalGenomicMap::divideCircuits() {
//    divide circuits to each partition, for common peace, divide averagely
    auto sources = mGraph->getMSources();
    int index = 0;
    int size = sources->size();
    for (int i = 0; i < size ; i ++) {
        auto seg = (*sources)[i];
        (*dividedCircuits)[seg->getId()] = new vector<VertexPath *>();
    }
    for (int i = 0; i < mCircuits->size(); i++) {
        auto circuit = (*mCircuits)[i];
        auto startId = (*circuit)[0]->getId();
        auto pair = this->findPartition(startId);
        if (pair.first != 0) {
            (*dividedCircuits)[pair.first]->push_back(circuit);
        } else {
            int partitionId = (*sources)[index % size]->getId();
            (*dividedCircuits)[partitionId]->push_back(circuit);
            index++;
        }
    }
}

void LocalGenomicMap::writeCircuits(const char *outFn) {
    cout << "Write circuits" << endl;
    ofstream fout(outFn);
    if (!fout) {
        cout << "Cannot open file " << outFn << ": no such file or directory" << endl;
        exit(1);
    }
    for (auto const &partitions : *dividedCircuits) {
        fout << "partition: " << partitions.first << "\n";
        for (auto circuits : *(partitions.second)) {
            for (Vertex *v : *circuits) {
                fout << v->getInfo() << " ";
            }
            fout << endl;
        }
    }
    fout.close();
}
//TODO if merge the virus path
void LocalGenomicMap::writeTraversedPath(const char *outFn) {
    cout << "Write Traversed Path" << endl;
    ofstream fout(outFn);
    if (!fout) {
        cout << "Cannot open file " << outFn << ": no such file or directory" << endl;
        exit(1);
    }
    for (auto &it : *traversedCircuits) {
        auto circuits = it.second;
        auto fSeg = mGraph->getSegmentById(it.first);
        fout <<fSeg->getChrom()<<"_"<<fSeg->getStart()<<":";
        for (auto circuit : *(circuits)) {
            fout << circuit->front()->getInfo()<<" ";
            int startId = circuit->front()->getId();
            for (auto v = circuit->begin() + 1; v != circuit->end() -1;) {
                if ((*v)->getId() == startId)
                    fout<<":";
                fout << (*v)->getInfo() << " ";
                v++;
            }
        }
        fout << "\n";
    }
}

void LocalGenomicMap::generateHaploids() {
    cout << "sort circuits" << endl;
    this->sortCircuits();
    srand(time(NULL));
    for (auto partitions : *dividedCircuits) {
        int partition = partitions.first;
        auto circuits = partitions.second;
        (*dividedHaploids)[partition] = new vector<VertexPath *>();
        vector<bool> is_inserted(circuits->size(), false);
        is_inserted[0] = true;
        VertexPath *mainPath = (*circuits)[0];
        int i = 1;
        VertexPath *current_circuit = new VertexPath();
        VertexPath *comp_circuit;
        bool is_comp = false;
        bool all_inserted = false;
        while (!all_inserted) {
            while (i < circuits->size()) {
                if (is_inserted[i]) {
                    i++;
                    continue;
                }
                // for (int i = 1; i < mCircuits->size(); i++) {
                if (!is_comp) {
                    current_circuit->assign(circuits->at(i)->begin(), circuits->at(i)->end());
                } else {
                    comp_circuit = this->get_complement(circuits->at(i));
                    current_circuit->assign(comp_circuit->begin(), comp_circuit->end());
                    // comp_circuit->clear();
                }
                deque<Vertex *> vq = deque<Vertex *>(current_circuit->begin(), current_circuit->end() - 1);
                bool isInserted = false;
                VertexPath::iterator cItr;
                VertexPath::iterator foundItr;
                for (int j = 0; j <= vq.size(); j++) {
                    Vertex *startV = vq.front();
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
                    this->print(*circuits->at(i));
                    all_inserted = false;
                    break;
                }
            }
        }
        (*dividedHaploids)[partition]->push_back(mainPath);
        cout << "Mainpath done" << endl;
        cout << (*dividedHaploids)[partition]->size() << endl;

        for (Vertex *v: *mainPath) {
            cout << v->getInfo() << " ";
        }
        cout << endl;
    }

    // VertexPath::iterator itr = find(mainPath->begin(), mainPath->end(), mGraph->getFirstSink()->getPositiveVertex());
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

void LocalGenomicMap::writeHaploids(const char *outFn) {
    ofstream fout(outFn);
    if (!fout) {
        cout << "Cannot open file " << outFn << " : no such file or directory" << endl;
        exit(1);
    }
    for (auto const &partitionHP: *dividedHaploids) {
        fout << "partition: " << partitionHP.first << "\n";
        auto hps = partitionHP.second;
        for (VertexPath *pathV : *hps) {
            for (Vertex *v : *pathV) {
                fout << v->getInfo() << " ";
            }
            fout << endl;
        }
    }
    fout.close();
}

void LocalGenomicMap::print(VertexPath &path) {
    for (Vertex *v : path) {
        cout << v->getInfo() << " ";
    }
    cout << endl;
}

void LocalGenomicMap::print(EdgePath &path) {
    for (Edge *e : path) {
        cout << e->getInfo() << " ";
    }
    cout << endl;
}

void LocalGenomicMap::printCircuits() {
    int count = 1;
    for (VertexPath *pathV : *mCircuits) {
        cout << count << " ---- ";
        count++;
        this->print(*pathV);
    }
}

void LocalGenomicMap::printHaploids() {
    for (VertexPath *pathV : *mHaploids) {
        this->print(*pathV);
    }
}

// Find BFB patterns/loops
void LocalGenomicMap::combinations(int start, int end, int len, vector<vector<int>> &res, vector<int> temp) {
	if (len == 0) {
        res.push_back(temp);
        return;
    }
    for (int i=start;i<=end;i++) {
        temp.push_back(i);
        combinations(i, end, len-1, res, temp);
        temp.pop_back();
    }
}

//comparison function for loops
bool compareLoops(vector<int> a, vector<int> b) {
    int diff1=0, diff2=0;
    if (a.size()>0 && b.size()>0) {
        diff1 = abs(a[0] - a[1]),
        diff2 = abs(b[0] - b[1]);
    }
    return (diff1>diff2); 
}

void LocalGenomicMap::constructDAG(vector<vector<int>> &adj, vector<vector<int>> &node2pat, 
                            vector<vector<int>> &node2loop, map<string, int> &variableIdx, int *elementCN) {
    vector<vector<int>> parents;
    for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) {
        if(elementCN[iter->second] > 0) {     
            cout<<iter->first<<" "<<elementCN[iter->second]<<endl;
            vector<int> temp;
            //push an empty vector<int>
            adj.push_back(temp), parents.push_back(temp);
            //split the string
            string key = iter->first;
            temp.push_back(stoi(key.substr(2, key.find(",")-2)));
            temp.push_back(stoi(key.substr(key.find(",")+1)));
            temp.push_back(elementCN[iter->second]);//copy number
            if (key[0] == 'p') {
                node2pat.push_back(temp);
                temp.clear();//synchronize the index
                node2loop.push_back(temp);
            }       
            else {
                node2loop.push_back(temp);
                temp.clear();
                node2pat.push_back(temp);
            }
        }
    }
    //sort all the loop in the decreasing order of length
    sort(node2loop.begin(), node2loop.end(), compareLoops);
    for (int i=0; i<adj.size(); i++) {
        if (node2pat[i].size()>0)
            cout<<"p "<<node2pat[i][0]<<" "<<node2pat[i][1]<<endl;
        else if(node2loop[i].size()>0)
            cout<<"l "<<node2loop[i][0]<<" "<<node2loop[i][1]<<endl;
    }
    //construct the adjacent lists
    for (int i=0; i<node2pat.size();i++) {
        if (node2pat[i].size()>0) {
            //p1 -> p2
            for (int j=0; j<node2pat.size();j++) {
                if (node2pat[j].size()>0 && (node2pat[i][0] == node2pat[j][0] ||
                    node2pat[i][1] == node2pat[j][1])) {
                    int diff1 = node2pat[i][0]-node2pat[i][1],
                        diff2 = node2pat[j][0]-node2pat[j][1];
                    if (abs(diff1)>abs(diff2)) {
                        adj[i].push_back(j);
                        parents[j].push_back(i);
                    }
                }
            }
            //p -> l
            for (int j=0; j<node2loop.size();j++) {
                if (node2loop[j].size()>0 && (node2pat[i][0] == node2loop[j][0] ||
                    node2pat[i][1] == node2loop[j][1])) {
                    int diff1 = node2pat[i][0]-node2pat[i][1],
                        diff2 = node2loop[j][0]-node2loop[j][1];
                    if (abs(diff1)>abs(diff2)) {
                        adj[i].push_back(j);
                        parents[j].push_back(i);
                    }
                }
            }
        }
    }
    for (int i=0; i<node2loop.size();i++) {
        if (node2loop[i].size()>0) {
            //l -> p
            for (int j=0; j<node2pat.size();j++) {
                if (find(parents[i].begin(), parents[i].end(), j) != parents[i].end()) // the pattern is parent of the loop
                    continue;
                if (node2pat[j].size()>0 && (node2loop[i][0] == node2pat[j][0] ||
                    node2loop[i][1] == node2pat[j][1])) {
                    int diff1 = node2loop[i][0]-node2loop[i][1],
                        diff2 = node2pat[j][0]-node2pat[j][1];
                    if(abs(diff1)>abs(diff2)) {
                        adj[i].push_back(j);
                        parents[j].push_back(i);
                    }
                    else {
                        for (int parent: parents[i]) {//the pattern is the sub-pattern of one of the loop's parents
                            if (find(adj[parent].begin(), adj[parent].end(), j) != adj[parent].end()) {
                                adj[i].push_back(j);
                                parents[j].push_back(i);
                                break;
                            }
                        }
                    }                    
                }
            }
            //l1 -> l2            
            for (int j=0; j<node2loop.size(); j++) {
                if (node2loop[j].size()>0 && (node2loop[i][0] == node2loop[j][0] || 
                    node2loop[i][1] == node2loop[j][1])) {//back
                    int diff1 = node2loop[i][0]-node2loop[i][1],
                        diff2 = node2loop[j][0]-node2loop[j][1];
                    if(abs(diff1)>abs(diff2)) {
                        adj[i].push_back(j);
                        parents[j].push_back(i);
                    }
                }
            }
        }
    }
}

void LocalGenomicMap::allTopologicalOrders(vector<int> &res, bool visited[], int num, int indeg[], vector<vector<int>> &adj, 
    vector<vector<int>> &orders) {
    //All the nodes are visited
    if (res.size() == num) {
        // for (int i: res)
        //     cout<<i+1<<" ";
        // cout<<endl;
        orders.push_back(res);
    }
    for (int i = 0; i < adj.size(); i++) {
        //If indegree is 0 and not yet visited then just choose node i
        if (indeg[i] == 0 && !visited[i]) {
            //Reduce indegree of adjacent vertices
            for (auto j = adj[i].begin(); j != adj[i].end(); j++)
                indeg[*j]--;
 
            //Push node i into the result vector
            res.push_back(i);
            visited[i] = true;
            allTopologicalOrders(res, visited, num, indeg, adj, orders);
 
            //Resetting visited, res and indegree
            visited[i] = false;
            res.erase(res.end() - 1);
            for (auto j = adj[i].begin(); j != adj[i].end(); j++)
                indeg[*j]++;
            //There is at least one unvisited node
        }
    }
}

void LocalGenomicMap::getBFB(vector<vector<int>> &orders, vector<vector<int>> &node2pat, 
                            vector<vector<int>> &node2loop, vector<int> &res) {
    bool forwardDir = true;                                
    for (int n=0; n<orders.size(); n++) {
        vector<int> bfb = orders[n];
        int i;
        vector<int> path;
        for (i=0;i<bfb.size();i++) {
            cout<<bfb[i]+1<<" ";
            if (node2pat[bfb[i]].size()) {
                if (path.size()==0) {
                    if (forwardDir) {
                        path.push_back(node2pat[bfb[i]][0]);
                        path.push_back(node2pat[bfb[i]][1]);
                    }
                    else {
                        path.push_back(node2pat[bfb[i]][1]);
                        path.push_back(node2pat[bfb[i]][0]);
                    }
                }
                else if (path.back() == node2pat[bfb[i]][0]) {
                    path.push_back(node2pat[bfb[i]][0]);
                    path.push_back(node2pat[bfb[i]][1]);
                }
                else if (path.back() == node2pat[bfb[i]][1]) {
                    path.push_back(node2pat[bfb[i]][1]);
                    path.push_back(node2pat[bfb[i]][0]);
                } 
                else
                    break; 
            }
            else if (node2loop[bfb[i]].size()) {
                cout<<"idx: "<<i+1<<" l:"<<node2loop[bfb[i]][0]<<","<<node2loop[bfb[i]][1]<<" CN: "<<node2loop[bfb[i]][2]<<endl;
                int cn = node2loop[bfb[i]][2];
                if (path.size()==0) {
                    while(cn--) {
                        if (forwardDir) {
                            path.push_back(node2loop[bfb[i]][0]);
                            path.push_back(node2loop[bfb[i]][1]);
                            path.push_back(node2loop[bfb[i]][1]);
                            path.push_back(node2loop[bfb[i]][0]);
                        }
                        else {
                            path.push_back(node2loop[bfb[i]][1]);
                            path.push_back(node2loop[bfb[i]][0]);
                            path.push_back(node2loop[bfb[i]][0]);
                            path.push_back(node2loop[bfb[i]][1]);
                        }
                    }
                    continue;
                }
                bool flag = true;
                while(cn--&&flag){
                    int idx = path.size()-1;
                    while(idx>=0) {//find a insertion position
                        if(path[idx]==node2loop[bfb[i]][0]&&path[idx]-path[idx-1]<0) {
                            if(idx == path.size()-1) break;
                            else if(path[idx-1]==path[idx+2]) break;
                        }
                        else if(path[idx]==node2loop[bfb[i]][1]&&path[idx]-path[idx-1]>0) {
                            if(idx == path.size()-1) break;
                            else if(path[idx-1]==path[idx+2]) break;
                        }
                        idx -= 2;
                    }
                    // while (idx>=0 && !(path[idx]== node2loop[bfb[i]][0]&&path[idx]-path[idx-1]<0) && 
                    //             !(path[idx]==node2loop[bfb[i]][1]&&path[idx]-path[idx-1]>0)) idx-=2;
                    if (idx<=0 || (idx<path.size()-1 && path[idx-1]!=path[idx+2])) {
                        flag = false;
                        break;
                    }
                    if (path[idx] == node2loop[bfb[i]][0] && path[idx]-path[idx-1]<0) {
                        path.insert(idx==path.size()-1?path.end():path.begin()+idx+1, {node2loop[bfb[i]][0], node2loop[bfb[i]][1], 
                            node2loop[bfb[i]][1], node2loop[bfb[i]][0]});
                    }
                    else if (path[idx] == node2loop[bfb[i]][1] && path[idx]-path[idx-1]>0) {
                        path.insert(idx==path.size()-1?path.end():path.begin()+idx+1, {node2loop[bfb[i]][1], node2loop[bfb[i]][0], 
                            node2loop[bfb[i]][0], node2loop[bfb[i]][1]});
                    }
                    else {
                        flag = false;
                        break;
                    }
                }
                if(flag == false) break;
            }            
        }
        cout<<"Order: ";
        for(int j=0; j<i; j++)
            cout<<bfb[j]+1<<" ";
        cout<<"\nPossible path: ";
        for (int n: path) {
            cout<<n<<" ";
        }
        cout<<endl;
        if(path.size()>res.size()) res = path;//find the longest result
        if (i == bfb.size()) {
            cout<<"Quit: "<<i<<" "<<bfb.size()<<endl;
            break;
        }
        else if (n == orders.size()-1 && forwardDir) {
            cout<<"Reverse\n";
            n = -1;
            forwardDir = false;
        }
    }
}

void LocalGenomicMap::readBFBProps(string &mainChr, int &insMode, vector<string> &insChr, int &conMode, vector<string> &conChr, 
                                    vector<int> &startSegs, const char *lhRawFn) {
    ifstream lhFile(lhRawFn);
    string line, prop;
    while (getline(lhFile, line)) {
        stringstream ss(line);
        ss>>prop;
        if (prop == "PROP") {
            while (ss>>prop) {
                int pos = 2, lastPos = 2;// ignore the first two chars e.g. "I:chr3:chr5"
                if (prop[0] == 'M')
                    mainChr = prop.substr(2);
                else if (prop[0] == 'I') {
                    if(prop[1]!=':') {// identify which mode (pre-BFB or post-BFB)
                        insMode = prop[1]-'0';
                        lastPos = 3;
                    }
                    else insMode = 2;               
                    while (pos != string::npos) {
                        pos = prop.find(":", lastPos);
                        insChr.push_back(prop.substr(lastPos, pos-lastPos));
                        lastPos = pos+1;
                    }
                }
                else if (prop[0] == 'C') {
                    if(prop[1]!=':') {
                        conMode = prop[1]-'0';
                        lastPos = 3;
                    }
                    else conMode = 2;
                    while (pos != string::npos) {
                        pos = prop.find(":", lastPos);
                        conChr.push_back(prop.substr(lastPos, pos-lastPos));
                        lastPos = pos+1;
                    }
                }
                else if (prop[0] == 'S') {
                    while (pos != string::npos) {
                        pos = prop.find(":", lastPos);
                        startSegs.push_back(stoi(prop.substr(lastPos, pos-lastPos)));
                        lastPos = pos+1;
                    }
                }
            }
        }
    }
    cout<<"Main chr: "<<mainChr<<endl;
    cout<<"Insertion chr: "<<endl;
    for (string ins: insChr)
        cout<<ins<<endl;
    cout<<"Concatenation chr: "<<endl;
    for (string con: conChr)
        cout<<con<<endl;
    cout<<"Starting seg: "<<endl;
    for (int seg: startSegs)
        cout<<seg<<endl;
}

void LocalGenomicMap::getJuncCN(vector<Junction *> &inversions , double** juncCN, Graph &graph,  int startSegID, int endSegID) {
    // find copy number for both normal junctions and fold-back inversions
    inversions.clear();
    for (int i=0; i <= endSegID; i++) {
        juncCN[i] = new double[2];//0: normal junction   1: inversion
        memset(juncCN[i], 0, 2*sizeof(double));
    }
    for (Junction *junc: *graph.getJunctions()) {
        int sourceID = junc->getSource()->getId(), targetID = junc->getTarget()->getId();
        char sourceDir = junc->getSourceDir(), targetDir = junc->getTargetDir();
        if (sourceID<startSegID||sourceID>endSegID||targetID<startSegID||targetID>endSegID) continue;
        double copyNum = junc->getWeight()->getCopyNum();
        if (0.5 < copyNum && copyNum < 1)
            copyNum = 1;// round small CN to 1
        if (sourceDir == targetDir) {// ht or th: not consider intra-chromosomal deletion
            if (sourceID+1 == targetID) {// normal edges
                juncCN[sourceID][0] += copyNum;
            }
            else if (sourceID-1 == targetID) {//normal edges (negative strand)
                juncCN[targetID][0] += copyNum;
            }
        }
        else {//hh or tt (inversion)
            if (abs(sourceID-targetID)<=3) {// fold-back inversion (with error of 3 bp)
                if(sourceID != targetID)
                    inversions.push_back(junc);
                if (sourceDir == '+') {
                    int greaterID = sourceID>targetID? sourceID:targetID;
                    juncCN[greaterID][1] += copyNum;
                }
                else {
                    int smallerID = sourceID<targetID? sourceID:targetID;
                    juncCN[smallerID][1] += copyNum;
                }                        
            }
        }            
    }
}

void LocalGenomicMap::insertBeforeBFB(Graph*& g, vector<string>& insChr, unordered_map<int, int>& originalSegs, vector<Junction *>& unusedSV) {
    unordered_map<int, int> segConversion;
    //segs, juncs, sources, sinks
    vector<Segment *> segs = *g->getSegments(), sources = *g->getMSources(), sinks = *g->getMSinks();
    vector<Junction *> juncs = *g->getJunctions();            
    vector<Segment *> mSegs, mSources, mSinks;
    vector<Junction *> mJuncs;
    //search for segments and junctions involved in insertion
    int sID, eID;
    vector<int> insertionIDs, deletedChrIDs;
    vector<Junction *> visited;
    for(int i=1; i<insChr.size(); i++) {
        for(int j=0; j<juncs.size(); j++) {
            if(find(visited.begin(), visited.end(), juncs[j])!=visited.end()) continue;
            Segment *seg1 = juncs[j]->getSource(), *seg2 = juncs[j]->getTarget();
            string chr1 = seg1->getChrom(), chr2 = seg2->getChrom();
            if((insChr[i-1]==chr1&&insChr[i]==chr2)||
                (insChr[i-1]==chr2&&insChr[i]==chr1)) {
                int id1 = seg1->getId(), id2 = seg2->getId();
                if(insChr[i-1]==chr2&&insChr[i]==chr1) swap(id1, id2);
                if(!insertionIDs.empty() && insertionIDs.back()!=id1) {
                    if(insertionIDs.back()<id1)
                        for(int i=insertionIDs.back(); i<id1; i++) insertionIDs.push_back(i);
                    else
                        for(int i=insertionIDs.back(); i>id1; i--) insertionIDs.push_back(i);
                }
                insertionIDs.push_back(id1), insertionIDs.push_back(id2);
                visited.push_back(juncs[j]);
                break;
            }
        }
    }            
    insertionIDs.erase(unique(insertionIDs.begin(), insertionIDs.end()), insertionIDs.end());//remove duplicates
    if(insertionIDs.front()>insertionIDs.back()) reverse(insertionIDs.begin(), insertionIDs.end());
    cout<<"Insertion ids: ";
    for(int id: insertionIDs) cout<<id<<" ";
    cout<<endl;
    sID = insertionIDs.front(), eID = insertionIDs.back();
    insertionIDs.erase(insertionIDs.begin());
    insertionIDs.pop_back();            
    juncs = *g->getJunctions();//reset juncs vector
    //set deleted chromosome IDs
    for(int id: insertionIDs) deletedChrIDs.push_back(segs[id-1]->getChrId());
    //set mSegs
    for(int i=1; i<=segs.size(); i++) {
        if(i<sID||i>eID) {
            if(find(deletedChrIDs.begin(),deletedChrIDs.end(), segs[i-1]->getChrId())!=deletedChrIDs.end()) continue;
            segConversion.insert(pair<int,int>(i, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[i-1]->getChrId(), segs[i-1]));
        }
        else {
            segConversion.insert(pair<int,int>(sID, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[sID-1]->getChrId(), segs[sID-1]));
            for(int j=sID+1; j<eID; j++) segConversion.insert(pair<int,int>(j, 0));//deleted segments
            for(int id: insertionIDs) {
                segConversion.insert(pair<int,int>(id, mSegs.size()+1));
                mSegs.push_back(new Segment(mSegs.size()+1, segs[sID-1]->getChrId(), segs[id-1]));
            }
            segConversion.insert(pair<int,int>(eID, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[eID-1]->getChrId(), segs[eID-1]));
            i = eID;
        }
    }
    //set mSources & mSinks
    mSources.push_back(mSegs[0]);
    for(int i=1; i<mSegs.size(); i++) {
        if(mSegs[i]->getChrId()!=mSegs[i-1]->getChrId()) {
            mSinks.push_back(mSegs[i-1]);
            mSources.push_back(mSegs[i]);
        }
    }
    mSinks.push_back(mSegs.back());
    //set mJuncs
    for(Junction* junc: juncs) {
        int startSegID = junc->getSource()->getId(), targetSegID = junc->getTarget()->getId();
        int id1 = segConversion[startSegID]-1, id2 = segConversion[targetSegID]-1;
        if(id1==-1 || id2 == -1) {
            unusedSV.push_back(junc);
            continue;
        }
        char dir1 = junc->getSourceDir(), dir2 = junc->getTargetDir();
        if(find(insertionIDs.begin(),insertionIDs.end(),startSegID)!=insertionIDs.end() ||
            find(insertionIDs.begin(),insertionIDs.end(),targetSegID)!=insertionIDs.end()) {
            if(id1>id2) swap(id1, id2);
            dir1 = '+', dir2 = '+';
        }
        cout<<startSegID<<"-"<<targetSegID<<" "<<id1+1<<"-"<<id2+1<<endl;
        mJuncs.push_back(new Junction(mSegs[id1],mSegs[id2],dir1,dir2,junc->getWeight()->getCoverage(),
            junc->getCredibility(),junc->getWeight()->getCopyNum(),junc->isInferred(),junc->hasLowerBoundLimit(),false));
    }
    //construct reversed reference for segments
    cout<<"Seg conversion:\n";
    for(auto iter=segConversion.begin(); iter!=segConversion.end(); iter++) {
        cout<<iter->first<<"-"<<iter->second<<endl;
        originalSegs.insert(pair<int,int>(iter->second, iter->first));
    }
    delete g;
    g = new Graph(mSegs, mJuncs, mSources, mSinks);
    g->writeGraph("./new.lh");
}

void LocalGenomicMap::insertAfterBFB(vector<string>& insChr, string& mainChr, vector<int>& startSegs,
    vector<vector<int>>& bfbPaths) {
    Graph* g = this->getGraph();
    //find other SVs based on graph of segments
    vector<Junction *> insertionSV, concatenationSV;
    vector<Segment *> segments = *g->getSegments();
    int segNum = segments.size()+1;
    vector<vector<int>> connections(segNum, vector<int>(segNum, -1));
    //find a range for normal links of each chromosome
    vector<int> maxSeg, minSeg;
    for (int i=0; i<g->getMSources()->size(); i++) {
        maxSeg.push_back((*g->getMSources())[i]->getId());
        minSeg.push_back((*g->getMSinks())[i]->getId());
    }
    for (Junction *junc: *g->getJunctions()) {
        if (junc->isInferred())
            continue;
        Segment *source = junc->getSource();
        Segment *target = junc->getTarget();
        int chr1 = source->getChrId(), chr2 = target->getChrId();
        if (chr1 != chr2) {
            int s = source->getId(), e = target->getId();
            if (s > maxSeg[chr1]) maxSeg[chr1] = s;
            if (s < minSeg[chr1]) minSeg[chr1] = s;
            if (e > maxSeg[chr2]) maxSeg[chr2] = e;
            if (e < minSeg[chr2]) minSeg[chr2] = e;
        }
    }        
    //get all SVs for insertion
    for (Junction *junc: *g->getJunctions()) {
        if (junc->isInferred())
            continue;
        Segment *source = junc->getSource();
        Segment *target = junc->getTarget();
        int chr1 = source->getChrId(), chr2 = target->getChrId();//index for chromosome
        string chrNum1 = source->getChrom(), chrNum2 = target->getChrom();//chromosome name
        int sourceId = source->getId(), targetId = target->getId();
        //sv for insertion
        if (find(insChr.begin(),insChr.end(),chrNum1)!=insChr.end() &&
            find(insChr.begin(),insChr.end(),chrNum2)!=insChr.end()) {
            if (chr1 != chr2) {
                insertionSV.push_back(junc);
                connections[sourceId][targetId] = insertionSV.size()-1;
                connections[targetId][sourceId] = insertionSV.size()-1;
            }
            else if (chrNum1 != mainChr) {//chr1==chr2: normal junctions for insertion
                if ((minSeg[chr1] <= sourceId && sourceId <= maxSeg[chr1]) &&
                    (minSeg[chr1] <= targetId && targetId <= maxSeg[chr1])) {
                    insertionSV.push_back(junc);
                    connections[sourceId][targetId] = insertionSV.size()-1;
                    connections[targetId][sourceId] = insertionSV.size()-1;
                }
            }
        }
    }
    cout<<"sv for insertion: "<<endl;
    for (Junction *junc: insertionSV) {
        cout<<junc->getSource()->getId()<<" "<<junc->getTarget()->getId()<<endl;
    }
    
    //construct bfb paths with insertion        
    for (int i=0; i<startSegs.size(); i++) {
        if (insertionSV.size() == 0) break;
        //find a path by traversing all the SVs with DFS
        vector<int> segs;//segments in sequence
        bool finished = false;
        bool *visited = new bool[segNum];
        memset(visited, false, segNum*sizeof(bool));
        stack<int> s;
        int startSeg = startSegs[i];
        int startChr = segments[startSeg-1]->getChrId();
        s.push(startSeg);//starting segment
        visited[startSeg] = true;
        while (!s.empty()) {
            int front = s.top();
            // cout<<front<<" "<<segments[front-1]->getChrId()<<endl;
            segs.push_back(front);
            s.pop();
            //adjacent segments
            for (int next=segNum-1; next>=1; next--) {
                if (!visited[next] && connections[front][next] != -1) {
                    s.push(next);
                    visited[next] = true;
                    if (segments[next-1]->getChrId() == startChr) {
                        segs.push_back(next);
                        finished = true;
                        break;
                    }
                }
            }
            if (finished)
                break;
        }
        if (segments[segs.back()-1]->getChrId() != startChr)
            segs.push_back(startSeg);
        cout<<"The sequence of segments"<<endl;
        for (int idx: segs)
            cout<<idx<<" ";
        cout<<endl;

        //construct the bool array for SV directions
        bool *edgeA = new bool[segs.size()];
        memset(edgeA, true, segs.size()*sizeof(bool));
        vector<Junction *> svPath;
        //find valid SVs in sequence
        int cnt = 0;
        for (int i=0; i<segs.size()-1; i+=1) {
            if(connections[segs[i]][segs[i+1]] == -1)
                break;
            Junction *sv = insertionSV[connections[segs[i]][segs[i+1]]];
            Edge *e = sv->getEdgeA();
            int chr1 = sv->getSource()->getChrId(), chr2 = sv->getTarget()->getChrId();
            //cout<<chr1<<" "<<chr2<<" "<<endl;
            if (chr1 == chr2) continue;
            svPath.push_back(sv);
            int sourceId = e->getSource()->getId(), targetId = e->getTarget()->getId();                
            if (sourceId == segs[i+1] && targetId == segs[i])
                edgeA[cnt] = false;
            cnt++;
        }
        //print bfb path with insertions
        if(svPath.empty())
            continue;
        cout<<"bfb path with insertions: "<<i<<endl;
        vector<int> res;
        this->bfbInsertion(svPath, bfbPaths, edgeA, res);
        vector<int> output;
        this->editBFB(bfbPaths, res, output);
        bfbPaths[res[0]] = output;
        //lgm->printBFB(output);
    }
}

void LocalGenomicMap::concatBeforeBFB(Graph*& g, vector<string>& conChr, unordered_map<int, int>& originalSegs, vector<Junction *>& unusedSV) {
    unordered_map<int, int> segConversion;
    //segs, juncs, sources, sinks
    vector<Segment *> segs = *g->getSegments(), sources = *g->getMSources(), sinks = *g->getMSinks();
    vector<Junction *> juncs = *g->getJunctions();            
    vector<Segment *> mSegs, mSources, mSinks;
    vector<Junction *> mJuncs;
    //search for segments and junctions involved in concatenation
    int sID, eID;
    char sDir, eDir;
    for(int i=0; i<juncs.size(); i++) {
        Junction *junc = juncs[i];
        if ((junc->getSource()->getChrom()==conChr[0] && junc->getTarget()->getChrom()==conChr[1]) ||
        (junc->getTarget()->getChrom()==conChr[0] && junc->getSource()->getChrom()==conChr[1])) {                    
            sID = junc->getSource()->getId(), eID = junc->getTarget()->getId();
            sDir = junc->getSourceDir(), eDir = junc->getTargetDir();
            // if(junc->getTarget()->getChrom().substr(0,3)=="chr") {
            //     swap(sID, eID);
            //     if(sDir == eDir) {
            //         if(sDir == '+') sDir = eDir = '-';
            //         else sDir = eDir = '+';
            //     }
            // }
            break;
        }
    }
    cout<<"Concat segs: "<<sID<<sDir<<" "<<eID<<eDir<<endl;
    //set mSegs
    int chrID1 = segs[sID-1]->getChrId();
    if(sDir == '+') {                
        for(int i=sources[chrID1]->getId(); i<=sID; i++) {
            segConversion.insert(pair<int,int>(i, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[sID-1]->getChrId(), segs[i-1]));
        }
        for(int i=sID+1; i<=sinks[chrID1]->getId(); i++) segConversion.insert(pair<int,int>(i, 0));
    }
    else {
        for(int i=sinks[chrID1]->getId(); i>=sID; i--) {
            segConversion.insert(pair<int,int>(i, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[sID-1]->getChrId(), segs[i-1]));
        }
        for(int i=sID-1; i>=sources[chrID1]->getId(); i--) segConversion.insert(pair<int,int>(i, 0));
    }
    int chrID2 = segs[eID-1]->getChrId();
    if(eDir == '+') {
        for(int i=eID; i<=sinks[chrID2]->getId(); i++) {
            segConversion.insert(pair<int,int>(i, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[sID-1]->getChrId(), segs[i-1]));
        }
        for(int i=sources[chrID2]->getId(); i<eID; i++) segConversion.insert(pair<int,int>(i, 0));
    }
    else {
        for(int i=eID; i>=sources[chrID2]->getId(); i--) {
            segConversion.insert(pair<int,int>(i, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[sID-1]->getChrId(), segs[i-1]));
        }
        for(int i=sinks[chrID2]->getId(); i>eID; i--) segConversion.insert(pair<int,int>(i, 0));
    }
    for(int i=1; i<=segs.size(); i++) {
        if(segs[i-1]->getChrId() != chrID1 && segs[i-1]->getChrId() != chrID2) {
            segConversion.insert(pair<int,int>(i, mSegs.size()+1));
            mSegs.push_back(new Segment(mSegs.size()+1, segs[i-1]->getChrId(), segs[i-1]));
        }
    }
    //set mSources & mSinks
    mSources.push_back(mSegs[0]);
    for(int i=1; i<mSegs.size(); i++) {
        if(mSegs[i]->getChrId()!=mSegs[i-1]->getChrId()) {
            mSinks.push_back(mSegs[i-1]);
            mSources.push_back(mSegs[i]);
        }
    }
    mSinks.push_back(mSegs.back());
    //set mJuncs
    for(Junction* junc: juncs) {
        int startSegID = junc->getSource()->getId(), targetSegID = junc->getTarget()->getId();
        int id1 = segConversion[startSegID]-1, id2 = segConversion[targetSegID]-1;
        char dir1 = junc->getSourceDir(), dir2 = junc->getTargetDir();
        cout<<startSegID<<dir1<<" - "<<targetSegID<<dir2<<" "<<id1+1<<"-"<<id2+1<<endl;
        if(id1==-1 || id2 == -1) {
            unusedSV.push_back(junc);
            continue;
        }
        if((startSegID==sID&&targetSegID==eID)||(startSegID==eID&&targetSegID==sID)) {
            if(id1>id2) swap(id1, id2);
            dir1 = '+', dir2 = '+';
        }
        mJuncs.push_back(new Junction(mSegs[id1],mSegs[id2],dir1,dir2,junc->getWeight()->getCoverage(),
            junc->getCredibility(),junc->getWeight()->getCopyNum(),junc->isInferred(),junc->hasLowerBoundLimit(),false));
    }
    //construct reversed reference for segments
    cout<<"Seg conversion:\n";
    for(auto iter=segConversion.begin(); iter!=segConversion.end(); iter++) {
        cout<<iter->first<<"-"<<iter->second<<endl;
        originalSegs.insert(pair<int,int>(iter->second, iter->first));
    }
    delete g;
    g = new Graph(mSegs, mJuncs, mSources, mSinks);
    g->writeGraph("./new.lh");
}

void LocalGenomicMap::concatAfterBFB(vector<string>& conChr, vector<vector<int>>& bfbPaths) {
    //find other SVs based on graph of segments
    vector<Junction *> concatenationSV;
    vector<Segment *> segments = *this->getGraph()->getSegments();
    int segNum = segments.size()+1;    
    for (Junction *junc: *this->getGraph()->getJunctions()) {
        if (junc->isInferred())
            continue;
        Segment *source = junc->getSource();
        Segment *target = junc->getTarget();
        int chr1 = source->getChrId(), chr2 = target->getChrId();//index for chromosome
        string chrNum1 = source->getChrom(), chrNum2 = target->getChrom();//chromosome name
        //sv for concatenation
        if (chr1 != chr2 && find(conChr.begin(),conChr.end(),chrNum1)!=conChr.end() &&
            find(conChr.begin(),conChr.end(),chrNum2)!=conChr.end()) {
            concatenationSV.push_back(junc);             
        }
    }
    cout<<"sv for concatenation: "<<endl;
    for (Junction *junc: concatenationSV) {
        cout<<junc->getSource()->getId()<<" "<<junc->getTarget()->getId()<<endl;
    }
    //construct bfb paths with concatenation        
    for (int i=0; i<concatenationSV.size(); i++) {
        cout<<"bfb paths with concatenation: "<<i<<endl;
        int chrID = concatenationSV[i]->getSource()->getChrId();
        int start = bfbPaths[chrID].size()-4<0? 0 : bfbPaths[chrID].size()-4;
        vector<int> res;
        this->bfbConcate(concatenationSV[i], true, start, 0, bfbPaths, res);//start from position 2 of the main chromosome
        if (!res.empty()) {
            vector<int> output;
            this->editBFB(bfbPaths, res, output);
            //lgm->printBFB(output);
        }         
    }
}

void LocalGenomicMap::editBFB(vector<vector<int>> bfbPaths, vector<int> &posInfo, vector<int> &output) {
    int chr = posInfo[0], pos = posInfo[1], segID = posInfo[2];
    output.insert(output.end(), bfbPaths[chr].begin(), bfbPaths[chr].begin()+pos+1);
    output.push_back(segID);
    for (int i=6; i<=posInfo.size()-4; i+=6) {
        int chr1 = posInfo[i-3], pos1 = posInfo[i-2], sourceID = posInfo[i-1],
            chr2 = posInfo[i], pos2 = posInfo[i+1], targetID = posInfo[i+2];
        //cout<<chr1<<" "<<pos1<<" "<<sourceID<<" "<<chr2<<" "<<pos2<<" "<<targetID<<endl;
        output.push_back(-1), output.push_back(-1);//sign for translocation
        output.push_back(sourceID);
        if (pos1 != pos2)
            output.insert(output.end(), bfbPaths[chr1].begin()+pos1+1, bfbPaths[chr2].begin()+pos2+1);
        output.push_back(targetID);
        //output.push_back(-1), output.push_back(-1);//sign for translocation 
    }
    output.push_back(-1), output.push_back(-1);
    chr = posInfo[posInfo.size()-3], pos = posInfo[posInfo.size()-2], segID = posInfo[posInfo.size()-1];
    output.push_back(segID);
    output.insert(output.end(), bfbPaths[chr].begin()+pos+1, bfbPaths[chr].end());
    //text output for visualization
    string bfbRes = "";
    vector<Junction*> juncs = *this->getGraph()->getJunctions();
    for (int i=1;i<output.size()-1;i+=2) {
        if (output[i] != -1)
            continue; 
        Junction* junc = NULL;
        for(Junction* j: juncs) {
            if((j->getSource()->getId()==output[i-2]&&j->getTarget()->getId()==output[i+1]) ||
            (j->getSource()->getId()==output[i+1]&&j->getTarget()->getId()==output[i-2])) {
                junc = j;
                break;
            }
        }
        if(junc != NULL) {
            char strand_5p = junc->getSourceDir(),
                strand_3p = junc->getTargetDir();
            if(junc->getSource()->getId()==output[i+1]) {
                if(strand_5p=='+'&&strand_3p=='+') {
                    strand_5p = '-';
                    strand_3p = '-';
                }
                else if(strand_5p=='-'&&strand_3p=='-') {
                    strand_5p = '+';
                    strand_3p = '+';
                }
            }
            if(bfbRes != "") bfbRes += "\t";
            bfbRes += to_string(junc->getWeight()->getCopyNum())+":"+to_string(output[i-2])+":"+strand_5p+":"+to_string(output[i+1])+":"+strand_3p;
        }
    }
    bfbRes += "\n";
    ofstream bfbFile;
    bfbFile.open("bfbPaths.txt",std::ios_base::app);
    bfbFile<<bfbRes;
    bfbFile.close();
    printBFB(output);
}

void LocalGenomicMap::editInversions(vector<int> &res, vector<Junction *> &inversions,
                            double** juncCN, int* elementCN, map<string, int> &variableIdx) {
    //text output for visualization
    string bfbRes = "";
    vector<Segment*> segs = *this->getGraph()->getSegments();
    bool isPositive = true;
    if(res[1]<res[0]) isPositive = false;
    //deal with imperfect inversions
    for (int j=1;j<res.size()-1;j+=2) {
        if (res[j] == -1)
            continue;        
        if (abs(res[j]-res[j-1])==abs(res[j+2]-res[j+1])) {
            string key = "l:"+to_string(min(res[j-1],res[j]))+","+to_string(max(res[j-1],res[j]));
            bfbRes += to_string(elementCN[variableIdx[key]])+":";
        }
        else
            bfbRes += "1:";
        for (Junction* junc: inversions) {
            int sourceSegID = junc->getSource()->getId(),
                targetSegID = junc->getTarget()->getId();
            char sourceDir = junc->getSourceDir();                        
            if ((sourceSegID==res[j] || res[j]==targetSegID) ||
                (targetSegID==res[j+1] || res[j+1]==sourceSegID)) {
                //check the copy number of junctions
                int segID = sourceSegID;
                if (sourceDir == '+') 
                    segID = sourceSegID>targetSegID? sourceSegID:targetSegID;
                else
                    segID = sourceSegID<targetSegID? sourceSegID:targetSegID;
                // if (juncCN[segID][1] > 0)
                //     juncCN[segID][1]--;
                // else
                //     continue;
                //edit the fold-back inversion
                if (sourceDir == '+') {
                    if (res[j]>res[j-1]) {
                        res[j] = sourceSegID;
                        res[j+1] = targetSegID;
                    }
                    else {
                        res[j] = targetSegID;
                        res[j+1] = sourceSegID;
                    }
                }
                else {
                    if (res[j]>res[j-1]) {
                        res[j] = targetSegID;
                        res[j+1] = sourceSegID;
                    }
                    else {
                        res[j] = sourceSegID;
                        res[j+1] = targetSegID;
                    }
                }
                break;
            }
        }
        char strand_5p = isPositive? '+':'-',
            strand_3p = isPositive? '-':'+';
        bfbRes += to_string(res[j])+":"+strand_5p+":"+to_string(res[j+1])+":"+strand_3p;
        if(j<res.size()-3) bfbRes += "\t";
        isPositive = !isPositive;
    }
    bfbRes += "\n";    
    ofstream bfbPaths;
    bfbPaths.open("bfbPaths.txt",std::ios_base::app);
    bfbPaths<<bfbRes;
    bfbPaths.close();
    printBFB(res);
}

void LocalGenomicMap::printBFB(vector<int> &res) {
    //text output for visualization
    string bfbRes = "";
    cout<<"find a BFB path"<<endl;
    vector<Segment *> segs = *this->getGraph()->getSegments();
    for (int j=0;j<res.size()-1;j+=2) {
        if (res[j] == -1)
            continue;
        cout<<res[j]<<"-"<<res[j+1]<<"|";
    }
    cout<<endl;
    string bedStr = "";
    for (int j=0;j<res.size()-1;j+=2) {
        if (res[j] == -1) { 
            bfbRes += "->";
            cout<<"->";
            continue;
        }
        if (res[j]<res[j+1]) {
            for (int k=res[j];k<=res[j+1];k++) {
                bfbRes += to_string(k);
                if(k<res[j+1]) bfbRes += ":";
                cout<<k<<":";
                //write down bp sequence
                bedStr += segs[k-1]->getChrom()+" "+to_string(segs[k-1]->getStart())+" "+to_string(segs[k-1]->getEnd())+" forward 1 +\n";
            }
        }
        else {
            for (int k=res[j];k>=res[j+1];k--) {
                bfbRes += to_string(k);
                if(k>res[j+1]) bfbRes += ":";
                cout<<k<<":";
                //write down bp sequence
                bedStr += segs[k-1]->getChrom()+" "+to_string(segs[k-1]->getStart())+" "+to_string(segs[k-1]->getEnd())+" reverse 1 -\n";
            }
        }
        if(j<res.size()-3&&res[j+2]!=-1) {
            bfbRes += "|";
            cout<<"|";
        }
    }
    bfbRes += "\n";
    cout<<endl;
    ofstream bfbPaths;
    bfbPaths.open("bfbPaths.txt",std::ios_base::app);
    bfbPaths<<bfbRes;
    bfbPaths.close();

    //output the bed file
    // ofstream bedFile;
    // bedFile.open("bed.txt");
    // bedFile<<bedStr;
    // bedFile.close();
}

void LocalGenomicMap::printOriginalBFB(vector<int> &res, unordered_map<int, int> &m, vector<Junction *> &unusedSV) {
    string bfbRes = "";
    cout<<"find a BFB path with integration first"<<endl;
    vector<Segment *> segs = *this->getGraph()->getSegments();
    for (int j=0;j<res.size()-1;j+=2) {
        if (res[j] == -1)
            continue;
        cout<<m[res[j]]<<"-"<<m[res[j+1]]<<"|";
    }
    cout<<endl;
    vector<vector<int>> path;
    for (int j=0;j<res.size()-1;j+=2) {
        vector<int> comp;
        if (res[j]<res[j+1]) {
            for (int k=res[j];k<=res[j+1];k++) {
                bfbRes += to_string(m[k]);
                comp.push_back(m[k]);
                if(k<res[j+1]) bfbRes += ":";
                cout<<m[k]<<":";
            }
        }
        else {
            for (int k=res[j];k>=res[j+1];k--) {
                bfbRes += to_string(m[k]);
                comp.push_back(m[k]);
                if(k>res[j+1]) bfbRes += ":";
                cout<<m[k]<<":";
            }
        }
        if(j<res.size()-3&&res[j+2]!=-1) {
            bfbRes += "|";
            cout<<"|";
        }
        path.push_back(comp);
    }
    cout<<endl;
    // second integration
    Junction* trans = nullptr;
    for(Junction *junc: unusedSV) {
        if(junc->getSource()->getChrom() != junc->getTarget()->getChrom()) {// unused translocation
            trans = junc;
            break;
        }
    }
    if(trans != nullptr) {
        // find strands
        vector<bool> isPositive;
        int s1 = path[0].front(), s2 = path[0][1], s3 = path[0].back();
        if(abs(s1-s2) == 1) isPositive.push_back(s1 < s2);
        else isPositive.push_back(s2 < s3);
        for(int i = 3; i < res.size(); i += 2) isPositive.push_back(!isPositive.back());
        int len = path.size();
        int sID = trans->getSource()->getId(), eID = trans->getTarget()->getId();
        char sDir = trans->getSourceDir(), eDir = trans->getTargetDir();
        cout<<sID<<sDir<<" "<<eID<<eDir<<endl;
        bool isSource = false;
        for(int i = 0; i < len; i++) {
            if(find(path[i].begin(), path[i].end(), sID) != path[i].end()) {
                isSource = true;
                break;
            }
        }
        int id = isSource? sID:eID, dir = isSource? sDir:eDir;
        bool isFinishd = false;
        for(int i = 0; i <= len/2 && isFinishd == false; i++) {
            for(int j = 0; j < 2 && isFinishd == false; j++) {
                int idx = j==0? i : len-i-1;
                auto iter = find(path[idx].begin(), path[idx].end(), id);
                if(iter != path[idx].end()) {
                    if(isSource == bool(j)) {
                        if((isPositive[idx]==true && dir=='+') || (isPositive[idx]==false && dir=='-')) {
                            if(isSource == true) {
                                vector<int> head(path[idx].begin(), iter+1);
                                head.push_back(eID);
                                path[idx] = head;
                                path.erase(path.begin()+idx+1, path.end());
                            }
                            else {
                                vector<int> tail(iter, path[idx].end());
                                tail.insert(tail.begin(), sID);
                                path[idx] = tail;
                                path.erase(path.begin(), path.begin()+idx);
                            }
                            isFinishd = true;
                        }
                    }
                    else {
                        if((isPositive[idx]==true && dir=='-') || (isPositive[idx]==false && dir=='+')) {
                            if(isSource == false) {
                                vector<int> head(path[idx].begin(), iter+1);
                                head.push_back(sID);
                                path[idx] = head;
                                path.erase(path.begin()+idx+1, path.end());
                            }
                            else {
                                vector<int> tail(iter, path[idx].end());
                                tail.insert(tail.begin(), eID);
                                path[idx] = tail;
                                path.erase(path.begin(), path.begin()+idx);
                            }
                            isFinishd = true;
                        }
                    }
                }
            }
        }
        // if(idx == len) { // cut the tail
        //     while(idx >= len/2) {
        //         auto iter = find(path[idx].begin(), path[idx].end(), eID);
        //         if(iter != path[idx].end()) {
        //             if((path[idx].back()-path[idx].front()>0 && eDir=='+') ||
        //             (path[idx].back()-path[idx].front()<0 && eDir=='-')) {
        //                 cout<<"Find subpath\n";
        //                 vector<int> tail(iter, path[idx].end());
        //                 tail.insert(tail.begin(), sID);
        //                 path[idx] = tail;
        //                 path.erase(path.begin(), path.begin()+idx);
        //                 break;
        //             }
        //         }
        //         idx--;
        //     }
        // }
        cout<<"BFB path after the second integration:\n";
        for(int i = 0; i < len; i++) {
            for(int j: path[i]) cout<<j<<":";
            if(i < len-1) cout<<"|";
        }
        cout<<endl;
    }

    bfbRes += "\n";
    ofstream bfbPaths;
    bfbPaths.open("bfbPaths.txt",std::ios_base::app);
    bfbPaths<<bfbRes;
    bfbPaths.close();
}

void LocalGenomicMap::BFB_ILP(const char *lpFn, vector<vector<int>> &patterns, vector<vector<int>> &loops, map<string, int> &variableIdx, 
                        double** juncCN, vector<vector<int>> &components, const bool juncsInfo, const double maxError, const bool seqMode) {
    // initialize variables for ILP
    OsiClpSolverInterface *si = new OsiClpSolverInterface();
    int startSegID = patterns.front()[0], endSegID = patterns.back()[1];
    cout<<"start-end: "<<startSegID<<" "<<endSegID<<endl;
    vector<Segment *> segs;
    for (Segment * seg: *this->getGraph()->getSegments()) {
        if (startSegID<=seg->getId() && seg->getId()<=endSegID)
            segs.push_back(seg);
    }
    
    int numSegments = segs.size();
    int numElements = variableIdx.size(), numEpsilons = numSegments*3;
    int numVariables = numElements+numEpsilons, numPat = patterns.size(), numLoop = loops.size();
    int numConstrains = (numSegments*2*3 + 3*numLoop)+(seqMode? numPat*numPat*3 : numPat*numPat);

    double *objective = new double[numVariables];
    double *variableLowerBound = new double[numVariables];
    double *variableUpperBound = new double[numVariables];
    double *constrainLowerBound = new double[numConstrains];
    double *constrainUpperBound = new double[numConstrains];
    cout << "Declare done" << endl;

    CoinPackedMatrix *matrix = new CoinPackedMatrix(false, 0, 0);
    int idx = 0;//number of constrains/inequalities
    //inequality formula: constrains on segment and junction CNs
    for (int i=startSegID;i<=endSegID;i++) {
        //??p + ??2*l + e_i >= c_i and 
        //??p + ??2*l - e_i <= c_i
        CoinPackedVector constrain1, constrain2;
        for (int j=0;j<numPat;j++) {//??p
            if (patterns[j][0]<=i && i<=patterns[j][1]) {
                string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                constrain1.insert(variableIdx[key], 1);
                constrain2.insert(variableIdx[key], 1);
            }
        }
        for (int j=0;j<numLoop;j++) {//??2*l
            if (loops[j][0]<=i && i<=loops[j][1]) {
                string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                constrain1.insert(variableIdx[key], 2);
                constrain2.insert(variableIdx[key], 2);
            }
        }
        constrain1.insert(numElements+idx/2, 1);//e_i
        constrainLowerBound[idx] = segs[i-startSegID]->getWeight()->getCopyNum();
        constrainUpperBound[idx] = si->getInfinity();
        idx++;
        matrix->appendRow(constrain1);

        constrain2.insert(numElements+idx/2, -1);//-e_i
        constrainLowerBound[idx] = -1*si->getInfinity();
        constrainUpperBound[idx] = segs[i-startSegID]->getWeight()->getCopyNum();
        idx++;
        matrix->appendRow(constrain2);

        //??p + ??2*l + e_ij >= c_ij where j=i+1 and
        //??p + ??2*l - e_ij <= c_ij where j=i+1
        CoinPackedVector constrain3, constrain4;
        for (int j=0;j<numPat;j++) {//??p
            if (patterns[j][0]<=i && i<patterns[j][1]) {
                string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                constrain3.insert(variableIdx[key], 1);
                constrain4.insert(variableIdx[key], 1);
            }
        }
        for (int j=0;j<numLoop;j++) {//??2*l
            if (loops[j][0]<=i && i<loops[j][1]) {
                string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                constrain3.insert(variableIdx[key], 2);
                constrain4.insert(variableIdx[key], 2);
            }
        }
        constrain3.insert(numElements+idx/2, 1);//e_ij
        constrainLowerBound[idx] = juncCN[i][0];//c_ij where j=i+1
        constrainUpperBound[idx] = si->getInfinity();
        idx++;
        matrix->appendRow(constrain3);

        constrain4.insert(numElements+idx/2, -1);//-e_ij
        constrainLowerBound[idx] = -1*si->getInfinity();
        constrainUpperBound[idx] = juncCN[i][0];//c_ij where j=i+1
        idx++;
        matrix->appendRow(constrain4);

        //??l + ??(p1+p2)/2 + ??(p+l)/2 + ??(l1+l2)/2 + e_ii >= c_ii and
        //??l + ??(p1+p2)/2 + ??(p+l)/2 + ??(l1+l2)/2 - e_ii <= c_ii
        double *coef = new double[numElements];
        memset(coef, 0, sizeof(double)*numElements);
        CoinPackedVector constrain5, constrain6;
        for (int j=0;j<numLoop;j++) {//??l
            if (loops[j][0] == i || loops[j][1] == i) {
                string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                coef[variableIdx[key]] += 1;
            }
        }
        for (int j=0;j<numPat;j++) {//??(p1+p2)/2
            for (int k=0;k<numPat;k++) {
                if ((patterns[j][0] == i && patterns[k][0] == i) ||
                    (patterns[j][1] == i && patterns[k][1] == i)) {
                    int diff1 = patterns[j][0]-patterns[j][1], diff2 = patterns[k][0]-patterns[k][1];
                    if (abs(diff1) > abs(diff2)) {
                        string key1 = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                        coef[variableIdx[key1]] += 0.5;
                        string key2 = "p:"+to_string(patterns[k][0])+","+to_string(patterns[k][1]);
                        coef[variableIdx[key2]] += 0.5;
                    }
                }
            }
        }
        for (int j=0;j<numPat;j++) {//??(p+l)/2
            for (int k=0;k<numLoop;k++) {
                if ((patterns[j][0] == i && loops[k][0] == i) ||
                    (patterns[j][1] == i && loops[k][1] == i)) {
                    int diff1 = patterns[j][0]-patterns[j][1], diff2 = loops[k][0]-loops[k][1];
                    if (abs(diff1) > abs(diff2)) {
                        string key1 = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                        coef[variableIdx[key1]] += 0.5;
                        string key2 = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                        coef[variableIdx[key2]] += 0.5;
                    }
                }
            }
        }
        for (int j=0;j<numLoop;j++) {//??(l1+l2)/2 (loop2 is inserted into the middle of loop1)
            for (int k=0;k<numLoop;k++) {
                if ((loops[j][0] == i && loops[k][0] == i) ||
                    (loops[j][1] == i && loops[k][1] == i)) {
                    int diff1 = loops[j][0]-loops[j][1], diff2 = loops[k][0]-loops[k][1];
                    if (abs(diff1) > abs(diff2)) {
                        string key1 = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                        coef[variableIdx[key1]] += 0.5;
                        string key2 = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                        coef[variableIdx[key2]] += 0.5;
                    }
                }
            }
        }
        for(int i=0;i<numElements;i++) {//insert the variables with coefficient>0
            if (coef[i] > 0.1) {
                constrain5.insert(i, coef[i]);
                constrain6.insert(i, coef[i]);
            }
        }
        cout<<"Inversion CN: "<<i<<"-"<<juncCN[i][1]<<endl;
        constrain5.insert(numElements+idx/2, 1);//e_ii
        constrainLowerBound[idx] = juncCN[i][1];//c_ii
        constrainUpperBound[idx] = si->getInfinity();
        idx++;
        matrix->appendRow(constrain5);

        constrain6.insert(numElements+idx/2, -1);//-e_ii
        constrainLowerBound[idx] = -si->getInfinity();
        constrainUpperBound[idx] = juncCN[i][1];//c_ii
        idx++;
        matrix->appendRow(constrain6);
        delete [] coef;
    }
    //constrains on errors
    if(maxError >= 0) {
        cout<<"Number of error variables: "<<idx/2<<endl;
        CoinPackedVector errorConstrain;
        for(int i=2; i<idx/2; i+=3) errorConstrain.insert(numElements+i, 1);
        constrainLowerBound[idx] = 0;
        constrainUpperBound[idx] = maxError;
        idx++;
        matrix->appendRow(errorConstrain);
    }

    //inequality formula: constrains on patterns and loops

    //0<=p(a,b)+p(c,d)<=1 and 0<=p(a,b)+l(c,d)<=1 and 0<=l(a,b)+l(c,d)<=1: constrains on exclusiveness
    for (int i=0;i<numPat;i++) {         
        for (int j=i+1;j<numPat;j++) {
            if ((patterns[i][0] < patterns[j][0] && patterns[i][1] < patterns[j][1]) ||
                (patterns[i][0] > patterns[j][0] && patterns[i][1] > patterns[j][1])) {
                CoinPackedVector constrain7, constrain8, constrain9;
                string key1 = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);                               
                string key2 = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                constrain7.insert(variableIdx[key1], 1);
                constrain7.insert(variableIdx[key2], 1);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = 1;
                idx++;
                matrix->appendRow(constrain7);
                if(seqMode) {
                    key2 = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                    constrain8.insert(variableIdx[key1], 1);
                    constrain8.insert(variableIdx[key2], 1);
                    constrainLowerBound[idx] = 0;
                    constrainUpperBound[idx] = 1;
                    idx++;
                    matrix->appendRow(constrain8);  

                    key1 = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
                    constrain9.insert(variableIdx[key1], 1);
                    constrain9.insert(variableIdx[key2], 1);
                    constrainLowerBound[idx] = 0;
                    constrainUpperBound[idx] = 1;
                    idx++;
                    matrix->appendRow(constrain9);
                } 
            }
        }
    }

    //??p(a,b)-p(c,d)>=0 where p(a,b) is a parent pattern of p(c,d)
    for (int i=0;i<numPat;i++) {
        CoinPackedVector constrain8;
        bool flag = false;
        for (int j=0;j<numPat;j++) {
            if ((patterns[i][0] == patterns[j][0]) ||
                (patterns[i][1] == patterns[j][1])) {
                int diff1 = patterns[i][0]-patterns[i][1], diff2 = patterns[j][0]-patterns[j][1];
                if (abs(diff1)<abs(diff2)) {
                    flag = true;
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain8.insert(variableIdx[key], 1);
                }
            }
        }
        if (flag) {
            string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
            constrain8.insert(variableIdx[key], -1);
            constrainLowerBound[idx] = 0;
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain8);
        }
    }
    int scale = 2;
    //??p(x1,y1)+??l(x2,y2)-l(a,b)>=0 where l(a,b) pattern is a sub-pattern of p and l
    for (int i=0;i<numLoop;i++) {
        CoinPackedVector constrain9;
        bool flag = false;
        for (int j=0;j<numPat;j++) {
            if ((patterns[j][0] == loops[i][0])||
                (patterns[j][1] == loops[i][1])) {
                int diff1 = loops[i][0]-loops[i][1], diff2 = patterns[j][0]-patterns[j][1];
                if (abs(diff1)<abs(diff2)) {
                    flag = true;
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain9.insert(variableIdx[key], 1);
                }
            }
        }
        for (int k=0;k<numLoop;k++) { // l(a,b) can be inserted into the middle of l(x2,y2)
            if(seqMode) break;
            if ((loops[k][0] == loops[i][0]) ||
                (loops[k][1] == loops[i][1])) {
                int diff1 = loops[i][0]-loops[i][1], diff2 = loops[k][0]-loops[k][1];
                if (abs(diff1)<abs(diff2)) {
                    flag = true;
                    string key = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                    constrain9.insert(variableIdx[key], scale);// scaling for copy number of small loop 
                }
            }
        }
        if (flag) {
            string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
            constrain9.insert(variableIdx[key], -1);//scaling coefficient
            constrainLowerBound[idx] = 0;
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain9);
        }
    }

    //third-generation information
    if(components.size()>0&&juncsInfo) {
        CoinPackedVector constrain10;
        for(int i=0; i<components.size(); i++) {     
            int s = min(components[i].front(), components[i].back()),
                e = max(components[i].front(), components[i].back());
            if(s==startSegID&&e==endSegID) continue; // ignore end-to-end case, e.g. 1-6
            cout<<"Round: "<<i+1<<"\t"<<s<<"-"<<e<<endl;
            string key = "l:"+to_string(s)+","+to_string(e);
            constrain10.insert(variableIdx[key], 1);
            key = "p:"+to_string(s)+","+to_string(e);
            constrain10.insert(variableIdx[key], 1);            
        }
        constrainLowerBound[idx] = 1;
        constrainUpperBound[idx] = 5;
        idx++;
        matrix->appendRow(constrain10);
        //p(a,d)+l(a,d)-??p(a,b)-??l(c,d)>=0: constrains on exclusiveness
        for (int i=0;i<numPat;i++) {
            CoinPackedVector constrain7, constrain8;
            bool flag = false;
            for (int j=0;j<numPat;j++) {
                if (patterns[i][0] == patterns[j][0] && patterns[i][1] > patterns[j][1]) {                
                    flag = true;          
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain7.insert(variableIdx[key], -1);
                    key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                    constrain8.insert(variableIdx[key], -1);
                }
                if (patterns[i][0] < patterns[j][0] && patterns[i][1] == patterns[j][1]) {
                    flag = true;          
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain8.insert(variableIdx[key], -1);
                    key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                    constrain7.insert(variableIdx[key], -1);
                }
            }
            if (flag) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                constrain7.insert(variableIdx[key], 1);
                constrain8.insert(variableIdx[key], 1);
                key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
                constrain7.insert(variableIdx[key], 1);
                constrain8.insert(variableIdx[key], 1);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = si->getInfinity();
                idx++;
                matrix->appendRow(constrain7);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = si->getInfinity();
                idx++;
                matrix->appendRow(constrain8);
            }
        }
    }    

    cout<<"ILP formula done"<<endl;

    //constrains for variables
    double maxCN = 0;
    for (Segment* seg: *this->getGraph()->getSegments()) {
        maxCN += seg->getWeight()->getCopyNum();
    }
    for (int i=0;i<numPat;i++) {//p
        string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
        variableLowerBound[variableIdx[key]] = 0;
        variableUpperBound[variableIdx[key]] = 1;
    }
    //reference pattern
    // string key = "p:"+to_string(startSegID)+","+to_string(endSegID);
    // variableLowerBound[variableIdx[key]] = variableUpperBound[variableIdx[key]] = 1;
    for (int i=0;i<numLoop;i++) {//l
        string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
        variableLowerBound[variableIdx[key]] = 0;
        variableUpperBound[variableIdx[key]] = maxCN;
    }
    for (int i=0;i<numEpsilons;i++) {//?? (or e)
        variableLowerBound[numElements+i] = 0;
        variableUpperBound[numElements+i] = si->getInfinity();
    }
    cout<<"Variable constrains done"<<endl;
    //weights for all variables
    for (int i=0;i<numVariables;i++) {
        if (i<numElements)
            objective[i] = 0;
        else {
            objective[i] = 1;
        }
    }

    si->loadProblem(*matrix, variableLowerBound, variableUpperBound, objective, constrainLowerBound,
                    constrainUpperBound);
    // set integer for junction variables
    for (int i = 0; i < numElements; i++) {
        si->setInteger(i);
    }

    si->writeMps(lpFn);
    si->writeLp(lpFn);

}

void LocalGenomicMap::BFB_ILP_SC(const char *lpFn, vector<vector<int>> &patterns, vector<vector<int>> &loops, map<string, int> &variableIdx, 
                        vector<Graph*> graphs, const double maxError, vector<vector<int>> &evolution) {
    int numGraphs = graphs.size();
    OsiClpSolverInterface *si = new OsiClpSolverInterface();
    int startSegID = patterns.front()[0], endSegID = patterns.back()[1];
    cout<<"start-end: "<<startSegID<<" "<<endSegID<<endl;    
    
    int numSegments = endSegID-startSegID+1;
    int numElements = variableIdx.size()*numGraphs, numEpsilons = numSegments*3*numGraphs+(numGraphs*(numGraphs-1))/2*variableIdx.size();
    int numVariables = numElements+numEpsilons, numPat = patterns.size(), numLoop = loops.size();
    int numComp = numPat+numLoop;
    int numConstrains = numSegments*2*3 + 3*numLoop + numPat;
    if (numGraphs > 1) numConstrains += (numGraphs*(numGraphs-1))*variableIdx.size()*2;
    numConstrains *= numGraphs;
    cout<<"Num of elements: "<<numElements<<" Num of epsilons: "<<numEpsilons<<endl;
    cout<<"Num of variables: "<<numVariables<<" Num of constrains: "<<numConstrains<<endl;
    double *objective = new double[numVariables];
    double *variableLowerBound = new double[numVariables];
    double *variableUpperBound = new double[numVariables];
    double *constrainLowerBound = new double[numConstrains];
    double *constrainUpperBound = new double[numConstrains];
    cout << "Declare done" << endl;

    CoinPackedMatrix *matrix = new CoinPackedMatrix(false, 0, 0);
    int idx = 0;//number of constrains/inequalities
    for(int n=0; n<numGraphs; n++) {
        vector<Segment *> segs;
        for (Segment * seg: *graphs[n]->getSegments()) {
            if (startSegID<=seg->getId() && seg->getId()<=endSegID)
                segs.push_back(seg);
        }
        //calculate the loop/pattern index
        if(n > 0) {
            for (int i=0; i<numPat; i++) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                variableIdx[key] += numComp;
                key = "l:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                variableIdx[key] += numComp;
            }
        }
        //find copy number for both normal junctions and fold-back inversions
        vector<Junction *> inversions;
        double** juncCN = new double*[endSegID+1];
        this->getJuncCN(inversions, juncCN, *graphs[n], startSegID, endSegID);
        //inequality formula: constrains on segment and junction CNs
        for (int i=startSegID;i<=endSegID;i++) {
            //??p + ??2*l + e_i >= c_i and 
            //??p + ??2*l - e_i <= c_i
            CoinPackedVector constrain1, constrain2;
            for (int j=0;j<numPat;j++) {//??p
                if (patterns[j][0]<=i && i<=patterns[j][1]) {
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain1.insert(variableIdx[key], 1);
                    constrain2.insert(variableIdx[key], 1);
                }
            }
            for (int j=0;j<numLoop;j++) {//??2*l
                if (loops[j][0]<=i && i<=loops[j][1]) {
                    string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                    constrain1.insert(variableIdx[key], 2);
                    constrain2.insert(variableIdx[key], 2);
                }
            }
            constrain1.insert(numElements+idx/2, 1);//e_i
            constrainLowerBound[idx] = segs[i-startSegID]->getWeight()->getCopyNum();
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain1);

            constrain2.insert(numElements+idx/2, -1);//-e_i
            constrainLowerBound[idx] = -1*si->getInfinity();
            constrainUpperBound[idx] = segs[i-startSegID]->getWeight()->getCopyNum();
            idx++;
            matrix->appendRow(constrain2);

            //??p + ??2*l + e_ij >= c_ij where j=i+1 and
            //??p + ??2*l - e_ij <= c_ij where j=i+1
            CoinPackedVector constrain3, constrain4;
            for (int j=0;j<numPat;j++) {//??p
                if (patterns[j][0]<=i && i<patterns[j][1]) {
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain3.insert(variableIdx[key], 1);
                    constrain4.insert(variableIdx[key], 1);
                }
            }
            for (int j=0;j<numLoop;j++) {//??2*l
                if (loops[j][0]<=i && i<loops[j][1]) {
                    string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                    constrain3.insert(variableIdx[key], 2);
                    constrain4.insert(variableIdx[key], 2);
                }
            }
            constrain3.insert(numElements+idx/2, 1);//e_ij
            constrainLowerBound[idx] = juncCN[i][0];//c_ij where j=i+1
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain3);

            constrain4.insert(numElements+idx/2, -1);//-e_ij
            constrainLowerBound[idx] = -1*si->getInfinity();
            constrainUpperBound[idx] = juncCN[i][0];//c_ij where j=i+1
            idx++;
            matrix->appendRow(constrain4);

            //??l + ??(p1+p2)/2 + ??(p+l)/2 + ??(l1+l2)/2 + e_ii >= c_ii and
            //??l + ??(p1+p2)/2 + ??(p+l)/2 + ??(l1+l2)/2 - e_ii <= c_ii
            double *coef = new double[numElements];
            memset(coef, 0, sizeof(double)*numElements);
            CoinPackedVector constrain5, constrain6;
            for (int j=0;j<numLoop;j++) {//??l
                if (loops[j][0] == i || loops[j][1] == i) {
                    string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                    coef[variableIdx[key]] += 1;
                }
            }
            for (int j=0;j<numPat;j++) {//??(p1+p2)/2
                for (int k=0;k<numPat;k++) {
                    if ((patterns[j][0] == i && patterns[k][0] == i) ||
                        (patterns[j][1] == i && patterns[k][1] == i)) {
                        int diff1 = patterns[j][0]-patterns[j][1], diff2 = patterns[k][0]-patterns[k][1];
                        if (abs(diff1) > abs(diff2)) {
                            string key1 = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                            coef[variableIdx[key1]] += 0.5;
                            string key2 = "p:"+to_string(patterns[k][0])+","+to_string(patterns[k][1]);
                            coef[variableIdx[key2]] += 0.5;
                        }
                    }
                }
            }
            for (int j=0;j<numPat;j++) {//??(p+l)/2
                for (int k=0;k<numLoop;k++) {
                    if ((patterns[j][0] == i && loops[k][0] == i) ||
                        (patterns[j][1] == i && loops[k][1] == i)) {
                        int diff1 = patterns[j][0]-patterns[j][1], diff2 = loops[k][0]-loops[k][1];
                        if (abs(diff1) > abs(diff2)) {
                            string key1 = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                            coef[variableIdx[key1]] += 0.5;
                            string key2 = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                            coef[variableIdx[key2]] += 0.5;
                        }
                    }
                }
            }
            for (int j=0;j<numLoop;j++) {//??(l1+l2)/2 (loop2 is inserted into the middle of loop1)
                for (int k=0;k<numLoop;k++) {
                    if ((loops[j][0] == i && loops[k][0] == i) ||
                        (loops[j][1] == i && loops[k][1] == i)) {
                        int diff1 = loops[j][0]-loops[j][1], diff2 = loops[k][0]-loops[k][1];
                        if (abs(diff1) > abs(diff2)) {
                            string key1 = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                            coef[variableIdx[key1]] += 0.5;
                            string key2 = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                            coef[variableIdx[key2]] += 0.5;
                        }
                    }
                }
            }
            for(int i=0;i<numElements;i++) {//insert the variables with coefficient>0
                if (coef[i] > 0.1) {
                    constrain5.insert(i, coef[i]);
                    constrain6.insert(i, coef[i]);
                }
            }
            cout<<"Inversion CN: "<<i<<"-"<<juncCN[i][1]<<endl;
            constrain5.insert(numElements+idx/2, 1);//e_ii
            constrainLowerBound[idx] = juncCN[i][1];//c_ii
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain5);

            constrain6.insert(numElements+idx/2, -1);//-e_ii
            constrainLowerBound[idx] = -si->getInfinity();
            constrainUpperBound[idx] = juncCN[i][1];//c_ii
            idx++;
            matrix->appendRow(constrain6);
            delete [] coef;
        }
        //constrains on errors
        if(maxError >= 0) {
            cout<<"Number of error variables: "<<idx/2<<endl;
            CoinPackedVector errorConstrain;
            for(int i=2; i<idx/2; i+=3) errorConstrain.insert(numElements+i, 1);
            constrainLowerBound[idx] = 0;
            constrainUpperBound[idx] = maxError;
            idx++;
            matrix->appendRow(errorConstrain);
        }

        //inequality formula: constrains on patterns and loops    
        //0<=p(a,b)+??p(c,d)<=1: constrains on exclusiveness
        for (int i=0;i<numPat;i++) {
            CoinPackedVector constrain7;
            bool flag = false;
            for (int j=0;j<numPat;j++) {
                if ((patterns[i][0] < patterns[j][0] && patterns[i][1] < patterns[j][1]) ||
                    (patterns[i][0] > patterns[j][0] && patterns[i][1] > patterns[j][1])) {                
                    flag = true;          
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain7.insert(variableIdx[key], 1);                
                }
            }
            if (flag) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                constrain7.insert(variableIdx[key], 1);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = 1;
                idx++;
                matrix->appendRow(constrain7);
            }
        }

        //p(c,d)-??p(a,b)>=0 where p(a,b) is a sub-pattern of p(c,d)
        for (int i=0;i<numPat;i++) {
            CoinPackedVector constrain8;
            bool flag = false;
            for (int j=0;j<numPat;j++) {
                if ((patterns[i][0] == patterns[j][0]) ||
                    (patterns[i][1] == patterns[j][1])) {
                    int diff1 = patterns[i][0]-patterns[i][1], diff2 = patterns[j][0]-patterns[j][1];
                    if (abs(diff1)>abs(diff2)) {
                        flag = true;
                        string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                        constrain8.insert(variableIdx[key], -1);
                    }
                }
            }
            if (flag) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                constrain8.insert(variableIdx[key], 1);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = si->getInfinity();
                idx++;
                matrix->appendRow(constrain8);
            }
        }
        int scale = 2;
        //??p(x1,y1)+??l(x2,y2)-l(a,b)>=0 where l(a,b) pattern is a sub-pattern of p and l
        for (int i=0;i<numLoop;i++) {
            CoinPackedVector constrain9;
            bool flag = false;
            for (int j=0;j<numPat;j++) {
                if ((patterns[j][0] == loops[i][0])||
                    (patterns[j][1] == loops[i][1])) {
                    int diff1 = loops[i][0]-loops[i][1], diff2 = patterns[j][0]-patterns[j][1];
                    if (abs(diff1)<abs(diff2)) {
                        flag = true;
                        string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                        constrain9.insert(variableIdx[key], 1);
                    }
                }
            }
            for (int k=0;k<numLoop;k++) { // l(a,b) can be inserted into the middle of l(x2,y2)
                if ((loops[k][0] == loops[i][0]) ||
                    (loops[k][1] == loops[i][1])) {
                    int diff1 = loops[i][0]-loops[i][1], diff2 = loops[k][0]-loops[k][1];
                    if (abs(diff1)<abs(diff2)) {
                        flag = true;
                        string key = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                        constrain9.insert(variableIdx[key], 1);
                    }
                }
            }
            if (flag) {
                string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
                constrain9.insert(variableIdx[key], -1);//scaling coefficient
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = si->getInfinity();
                idx++;
                matrix->appendRow(constrain9);
            }
        }   

        cout<<"ILP formula done"<<endl;

        //constrains for variables
        double maxCN = 0;
        for (Segment* seg: segs) {
            maxCN += seg->getWeight()->getCopyNum();
        }
        for (int i=0;i<numPat;i++) {//p
            string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
            variableLowerBound[variableIdx[key]] = 0;
            variableUpperBound[variableIdx[key]] = 1;
        }
        //reference pattern
        // string key = "p:"+to_string(startSegID)+","+to_string(endSegID);
        // variableLowerBound[variableIdx[key]] = variableUpperBound[variableIdx[key]] = 1;
        for (int i=0;i<numLoop;i++) {//l
            string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
            variableLowerBound[variableIdx[key]] = 0;
            variableUpperBound[variableIdx[key]] = maxCN;
        }        
    }
    //constrains on common components
    cout<<"Number of constrains: "<<idx<<endl;
    int cnt = (numElements+numSegments*3*numGraphs)*2;//idx for remaining epsilons
    cout<<"Stat cnt: "<<cnt/2<<endl;
    for(auto iter: variableIdx) variableIdx[iter.first] = iter.second%numComp;
    for(int i=0; i<evolution.size(); i++) {
        for(int j: evolution[i]) {
            for(int k=0; k<numPat; k++) {
                CoinPackedVector constrain1, constrain2;
                string key = "p:"+to_string(patterns[k][0])+","+to_string(patterns[k][1]);
                constrain1.insert(variableIdx[key]+numComp*i, 1);
                constrain1.insert(variableIdx[key]+numComp*j, -1);
                constrain1.insert(cnt/2, 1);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = si->getInfinity();
                idx++, cnt++;
                matrix->appendRow(constrain1);
                constrain2.insert(variableIdx[key]+numComp*i, 1);
                constrain2.insert(variableIdx[key]+numComp*j, -1);
                constrain2.insert(cnt/2, -1);
                constrainLowerBound[idx] = -si->getInfinity();
                constrainUpperBound[idx] = 0;
                idx++, cnt++;
                matrix->appendRow(constrain2);
            }
            for (int k=0;k<numLoop;k++) {
                CoinPackedVector constrain1, constrain2;
                string key = "l:"+to_string(loops[k][0])+","+to_string(loops[k][1]);
                constrain1.insert(variableIdx[key]+numComp*i, 1);
                constrain1.insert(variableIdx[key]+numComp*j, -1);
                constrain1.insert(cnt/2, 1);
                constrainLowerBound[idx] = 0;
                constrainUpperBound[idx] = si->getInfinity();
                idx++, cnt++;
                matrix->appendRow(constrain1);
                constrain2.insert(variableIdx[key]+numComp*i, 1);
                constrain2.insert(variableIdx[key]+numComp*j, -1);
                constrain2.insert(cnt/2, -1);
                constrainLowerBound[idx] = -si->getInfinity();
                constrainUpperBound[idx] = 0;
                idx++, cnt++;
                matrix->appendRow(constrain2);
            }            
        }
    }
    cout<<"End cnt: "<<cnt/2<<endl;      
    for (int i=0;i<numEpsilons;i++) {//?? (or e)
        variableLowerBound[numElements+i] = 0;
        variableUpperBound[numElements+i] = si->getInfinity();
    }
    //weights for all variables
    for (int i=0;i<numVariables;i++) {
        if (i<numElements)
            objective[i] = 0;
        else   
            objective[i] = 1;
    }    
    cout<<"Variable constrains done"<<endl;
    si->loadProblem(*matrix, variableLowerBound, variableUpperBound, objective, constrainLowerBound,
                    constrainUpperBound);
    // set integer for junction variables
    for (int i = 0; i < numElements; i++) {
        si->setInteger(i);
    }
    si->writeMps(lpFn);
    si->writeLp(lpFn);
}

//Deal with complex case: concatenation
void LocalGenomicMap::bfbConcate(Junction *sv, bool edgeA, int pos1, int pos2, vector<vector<int>> bfbPaths, vector<int> &res) {
    Edge *e = edgeA? sv->getEdgeA():sv->getEdgeB();
    //find the breakage positions and merge bfb paths
    vector<Segment *> segs = *this->getGraph()->getSegments();
    bool found = false;
    int sourceID = e->getSource()->getId(), targetID = e->getTarget()->getId();
    char sourceDir = e->getSource()->getDir(), targetDir = e->getTarget()->getDir();
    int chr1 = segs[sourceID-1]->getChrId(), chr2 = segs[targetID-1]->getChrId();
    // cout<<(edgeA?"True":"False")<<endl;
    // cout<<sourceID<<sourceDir<<" -> "<<targetID<<targetDir<<endl;         
    for (int i=pos1; i<bfbPaths[chr1].size()-1; i+=2) {
        if ((bfbPaths[chr1][i] <= sourceID && sourceID <= bfbPaths[chr1][i+1] && sourceDir == '+') ||
            (bfbPaths[chr1][i] >= sourceID && sourceID >= bfbPaths[chr1][i+1] && sourceDir == '-')){                          
            pos1 = i;                        
            for (int j=pos2; j<bfbPaths[chr2].size()-1; j+=2) {
                if ((bfbPaths[chr2][j] <= targetID && targetID <= bfbPaths[chr2][j+1] && targetDir == '+') ||
                    (bfbPaths[chr2][j] >= targetID && targetID >= bfbPaths[chr2][j+1] && targetDir == '-')) {
                    pos2 = j;                               
                    found = true;                          
                    break;
                }
            }
        } 
        if (found)
            break;
    }
    if (found) {
        //cout<<chr1<<": "<<pos1<<" -> "<<chr2<<": "<<pos2<<endl;
        res.push_back(chr1), res.push_back(pos1), res.push_back(sourceID);
        res.push_back(chr2), res.push_back(pos2), res.push_back(targetID);
    }
    else 
        cout<<"Cannot concatenate the two bfb paths."<<endl;    
}

// deal with complex case: insertion
void LocalGenomicMap::bfbInsertion(vector<Junction *> &SVs, vector<vector<int>> bfbPaths, bool edgeA[], vector<int> &res) {
    int pos1 = 0, pos2 = 0;
    for (int i=0; i<SVs.size(); i++) {
        vector<int> copy = res;
        this->bfbConcate(SVs[i], edgeA[i], pos1, pos2, bfbPaths, res);
        int n = res.size();
        if (copy.size() == n) {
            cout<<"Cannot insert all the SVs into the main chromosome."<<endl;
            break;
        }
        pos1 = res[n-2];
        pos2 = i==SVs.size()-2? res[1] : 0;
    }    
}

void LocalGenomicMap::readComponents(vector<vector<int>>& res, const char *juncsFn) {
    ifstream inFile(juncsFn);
    string line;
    Graph* g = this->getGraph();
    if(juncsFn=="") return;
    while (getline(inFile, line)) {
        istringstream iss(line);
        vector<int> segs;
        vector<char> sign;
        string temp;
        while(iss >> temp) {
            segs.push_back(stoi(temp.substr(0,temp.length()-1)));
            sign.push_back(temp.back());
        }
        // process inversion, translocation, and component
        int lastIdx = 0;
        for(int i=1; i<segs.size(); i++) {
            if(g->getSegmentById(segs[lastIdx])->getPartition() != g->getSegmentById(segs[i])->getPartition() ||
                sign[i-1] != sign[i]) {// breakpoints                        
                // add component
                if(i-lastIdx>=2) {
                    vector<int> subset;
                    subset.assign(segs.begin()+lastIdx, segs.begin()+i);
                    sort(subset.begin(), subset.end());
                    res.push_back(subset);
                }
                // add junction (inversion or translocation)
                int sourceId = segs[i-1], targetId = segs[i];
                char sourceDir = sign[i-1], targetDir = sign[i];
                cout<<sourceId<<sourceDir<<" -> "<<targetId<<targetDir<<endl;
                double junCoverage = g->getAvgCoverage();
                Junction* junc1 = new Junction(g->getSegmentById(sourceId), g->getSegmentById(targetId), sourceDir, targetDir, junCoverage, 1, 1, false, true, false);
                Junction* junc2 = g->findJunction(junc1);
                if(junc2 == NULL) {
                    g->addJunction(sourceId, sourceDir, targetId, targetDir, junCoverage, 1, 1, false, true, false);
                }
                else {
                    junc2->getWeight()->setCopyNum(1);
                }
                lastIdx = i;
            }
        }
        if(segs.size()-lastIdx>=2) {
            // add trailing component
            vector<int> subset;
            subset.assign(segs.begin()+lastIdx, segs.end());
            sort(subset.begin(), subset.end());
            res.push_back(subset);
        }
    }
    // remove duplicates
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
}

/*
* old version
*/
VertexPath* LocalGenomicMap::findBFB(VertexPath* currPath, int n, set<Edge *>* visited, int error) {
    int len = currPath->size();
    if (len == n)
        return currPath;
    if (len > n)
        return NULL;
    VertexPath* BFBpath;
    Vertex* lastV = currPath->back();
    EdgePath* nextEdges = lastV->getEdgesAsSource();
    // record the original BFB pattern and visit set
    VertexPath initPath = *currPath;
    set<Edge *> initVisited = *visited;
    cout<<"Current path size: "<<currPath->size()<<endl;
    for (Edge *e: *nextEdges) {
        *currPath = initPath;
        *visited = initVisited;
        if (visited->find(e) == visited->end() && this->checkBFB(currPath, e->getTarget())) {
            visited->insert(e);
            // cout<<e->getTarget()->getId()<<e->getTarget()->getDir()<<" "<<&e<<endl;
            currPath->push_back(e->getTarget());
            BFBpath = findBFB(currPath, n, visited, error);
            if (BFBpath != NULL)
                return BFBpath;
        }
    }
    if (nextEdges->size()>0 && error>0) {// deal with special cases
        //randomness
        random_device rd;
        mt19937 g(rd());
        shuffle(nextEdges->begin(), nextEdges->end(), g);
        Edge* e = nextEdges->front();
        //initialization
        *currPath = initPath;
        *visited = initVisited;
        error--;
        cout<<"original: "<<visited->size()<<" "<<currPath->size()<<endl;
        // insertion
        if (visited->find(e) == visited->end()) {
            visited->insert(e);
            currPath->push_back(e->getTarget());
            BFBpath = findBFB(currPath, n, visited, error);
            if (BFBpath != NULL)
                return BFBpath;
        }
        // deletion
        *visited = initVisited;
        *currPath = initPath;
        cout<<"before deletion: "<<visited->size()<<" "<<currPath->size()<<endl;
        currPath->pop_back();
        BFBpath = findBFB(currPath, n, visited, error);
        if (BFBpath != NULL)
            return BFBpath;
        //replacement
        if (visited->find(e) == visited->end()) {
            *visited = initVisited;
            *currPath = initPath;
            cout<<"before replacement: "<<visited->size()<<" "<<currPath->size()<<endl;
            visited->insert(e);
            currPath->pop_back();
            currPath->push_back(e->getTarget());
            BFBpath = findBFB(currPath, n, visited, error);
            if (BFBpath != NULL)
                return BFBpath; 
        }
    }
    return NULL;
}

bool LocalGenomicMap::checkBFB(VertexPath* currPath, Vertex* v) {
    Vertex* lastV = currPath->back();
    cout<<v->getId()<<v->getDir()<<" ";
    int len = currPath->size();
    for (int i=1;i<len;i++) {
        if (this->isPalindrome(currPath, i)) {
            Vertex* preV = (i==len-1)?(*currPath)[i] : (*currPath)[i-1];
            if (preV->isReverse(v)) {
                return true;
            }
        }
    }
    return false;
}

bool LocalGenomicMap::isPalindrome(VertexPath* path, int start) {
    int len = path->size();
    int bound = (start+len)/2;
    for (int i=start;i<bound;i++){
        if (!(*path)[i]->isReverse((*path)[--len]))
            return false;
    }
    return true;
}

//Bipartite matching
void LocalGenomicMap::findMaxBPMatching(vector<Junction *> &juncs, vector<Junction *> &results) {
    //get source and target segments & construct connection matrix
    set<Segment *> sources, targets;
    int num = this->getGraph()->getSegments()->size()+1;
    bool** connection = new bool*[num];
    for (int i=0; i < num; i++) {
        connection[i] = new bool[num];
        memset(connection[i], false, sizeof(connection[i]));
    }
    for (Junction *junc: juncs){
        sources.insert(junc->getSource());
        targets.insert(junc->getTarget());
        connection[junc->getSource()->getId()][junc->getTarget()->getId()] = true;  
    }
    int* match = new int[num];
    memset(match, -1, sizeof(int)*num);

    cout<<"Search for the max BPM..."<<endl;
    int result = 0;
    bool* visited = new bool[num];
    for (Segment *s: sources) {
        memset(visited, false, sizeof(visited));
        cout<<"Check source: "<<s->getId()<<endl;
        if(this->bpm(connection, s, targets, visited, match))
            result++;
    }
    for (int i=0; i<num; i++) {
        cout<<match[i]<<" ";
        if (match[i] > 0) {
            for (Junction *junc: juncs) {
                if (junc->getTarget()->getId() == i && junc->getSource()->getId() == match[i])
                    results.push_back(junc);
            }
        }
    }
    cout<<"Max BP match: "<<result<<endl;
}

bool LocalGenomicMap::bpm(bool** connection, Segment *source, set<Segment *> &targets, bool visited[], int match[]) {
    for (Segment *t: targets) {
        if (connection[source->getId()][t->getId()] && !visited[t->getId()]) {
            visited[t->getId()] = true;
            //if the target is occupied, check other matches
            Segment *next; 
            if (match[t->getId()] > 0) {
                next = this->getGraph()->getSegmentById(match[t->getId()]);
            }
            if (match[t->getId()]<=0 || bpm(connection, next, targets, visited, match)) {
                match[t->getId()] = source->getId();
                cout<<source->getId()<<" -> "<<t->getId()<<endl;
                return true;
            }
        }
    }
    return false;
}

//Hierholzer???s Algorithm for finding an Euler circuit in a directed graph
void LocalGenomicMap::findCircuits(vector<vector<int>> adj) {
    unordered_map<int, int> edgeCount;//number of edges for each starting segment
    for (int i=1; i<=adj.size(); i++)
        edgeCount[i] = adj[i].size();
    
    stack<int> currPath;
    vector<int> circuit;
    currPath.push(1);
    int currSeg = 1;
    
    while (!currPath.empty()) {
        if (edgeCount[currSeg]) {
            currPath.push(currSeg);
            int nextSeg = adj[currSeg].back();
            edgeCount[currSeg]--;
            adj[currSeg].pop_back();
            currSeg = nextSeg;    
        }
        else {
            circuit.push_back(currSeg);
            currSeg = currPath.top();
            currPath.pop();
        }
    }

    for (int i=circuit.size()-1; i>=0; i--) {
        cout<< circuit[i];
        if (i > 0)
            cout<< " -> ";
        else
            cout<<endl;
    }
}

//Construct circuits/paths
void LocalGenomicMap::constructCircuits(vector<vector<int>> sv) {
    vector<vector<int>> res;
    res.push_back(sv[0]);
    sv.erase(sv.begin());    
    while (sv.size()) {
        bool finished = true; 
        int idx = 0;       
        for (vector<int> junc: sv) {
            if (junc[0] == res.back()[1]) {
                finished = false;
                res.push_back(junc);
                sv.erase(sv.begin()+idx);
            }
            else if (junc[1] == res.front()[0]) {
                finished = false;
                res.insert(res.begin(), junc);
                sv.erase(sv.begin()+idx);
            }
            idx++;
        }
        if (finished) break;
    }
    cout<<res[0][0]<<"->"<<res[0][1]<<" ";
    for (int i=1; i<res.size(); i++) {
        cout<<res[i][1]<<"->";
    }
    cout<<endl;
}