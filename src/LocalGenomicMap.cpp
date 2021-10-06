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
//            如果seg是病毒seg那么对任意第一个定点reachable就可以, 否则需要在自己的partition内reachable
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
// //                    TODO 两个check需要更改
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
//                    TODO 两个check需要更改
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

// Find BFB patterns
void LocalGenomicMap::combinations(int start, int end, int len, vector<vector<int>> &res, vector<int> temp) {
	if (len == 0) {
        res.push_back(temp);
        return;
    }
    for (int i=start;i<=end;i++) {
        temp.push_back(i);
        combinations(i+1, end, len-1, res, temp);
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

void LocalGenomicMap::constructDAG(vector<vector<int>> &adj, bool** mLoop, vector<vector<int>> &node2pat, 
                            vector<vector<int>> &node2loop, map<string, int> &variableIdx, int *elementCN) {
    vector<vector<int>> parents;
    for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) {
        if(elementCN[iter->second] > 0) {     
            vector<int> temp;
            //push an empty vector<int>
            adj.push_back(temp), parents.push_back(temp);
            //split the string
            string key = iter->first;
            temp.push_back(stoi(key.substr(2, key.find(",")-2)));
            temp.push_back(stoi(key.substr(key.find(",")+1)));
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

void LocalGenomicMap::printBFB(vector<vector<int>> &orders, vector<vector<int>> &node2pat, vector<vector<int>> &node2loop) {
    vector<int> res;
    for (vector<int> bfb: orders) {
            int i;
            res.clear();
            for (i=0;i<bfb.size();i++) {
                //cout<<bfb[i]+1<<" ";
                if (node2pat[bfb[i]].size()) {
                    if (res.size()==0) {
                        res.push_back(node2pat[bfb[i]][0]);
                        res.push_back(node2pat[bfb[i]][1]);
                    }
                    else if (res.back() == node2pat[bfb[i]][0]) {
                        res.push_back(node2pat[bfb[i]][0]);
                        res.push_back(node2pat[bfb[i]][1]);
                    }
                    else if (res.back() == node2pat[bfb[i]][1]) {
                        res.push_back(node2pat[bfb[i]][1]);
                        res.push_back(node2pat[bfb[i]][0]);
                    } 
                    else
                        break; 
                }
                else if (node2loop[bfb[i]].size()) {
                    // lgm->printLoop(node2pat, node2loop, mLoop, bfb[i]);
                    if (res.size()==0) {
                        res.push_back(node2loop[bfb[i]][0]);
                        res.push_back(node2loop[bfb[i]][1]);
                        res.push_back(node2loop[bfb[i]][1]);
                        res.push_back(node2loop[bfb[i]][0]);
                        continue;
                    }
                    int idx = res.size()-1;
                    while (idx>=0 && res[idx] != node2loop[bfb[i]][0] && 
                            res[idx] != node2loop[bfb[i]][1]) idx--;
                    if (idx<0)
                        break;
                    if (res[idx] == node2loop[bfb[i]][0]) {
                        res.insert(idx==res.size()-1?res.end():res.begin()+idx, {node2loop[bfb[i]][0], node2loop[bfb[i]][1], 
                            node2loop[bfb[i]][1], node2loop[bfb[i]][0]});
                    }
                    else if (res[idx] == node2loop[bfb[i]][1]) {
                        res.insert(idx==res.size()-1?res.end():res.begin()+idx, {node2loop[bfb[i]][1], node2loop[bfb[i]][0], 
                            node2loop[bfb[i]][0], node2loop[bfb[i]][1]});
                    }
                }            
            }
            if (i == bfb.size()) {
                cout<<"find a BFB path"<<endl;
                // for (int j=0;j<res.size();j++)
                //     cout<<res[j]<<" ";
                // cout<<endl;
                for (int j=0;j<res.size()-1;j+=2) {
                    if (res[j]<res[j+1]) {
                        for (int k=res[j];k<=res[j+1];k++)
                            cout<<k;
                    }
                    else {
                        for (int k=res[j];k>=res[j+1];k--)
                            cout<<k;
                    }
                    cout<<"|";
                }
                break;
            }                
        }
}

void LocalGenomicMap::BFB_ILP(const char *lpFn, vector<vector<int>> &patterns, 
                            vector<vector<int>> &loops, map<string, int> &variableIdx, double** juncCN) {
    OsiClpSolverInterface *si = new OsiClpSolverInterface();
    vector<Segment *> *segs = this->getGraph()->getSegments();
    int numSegments = segs->size();
    int numElements = variableIdx.size(), numEpsilons = numSegments*3;
    int numVariables = numElements+numEpsilons, numPat = patterns.size(), numLoop = loops.size();
    int numConstrains = numSegments*2*3 + 2*numPat + 2*numLoop;

    double *objective = new double[numVariables];
    double *variableLowerBound = new double[numVariables];
    double *variableUpperBound = new double[numVariables];
    double *constrainLowerBound = new double[numConstrains];
    double *constrainUpperBound = new double[numConstrains];
    cout << "Declare done" << endl;

    CoinPackedMatrix *matrix = new CoinPackedMatrix(false, 0, 0);
    int idx = 0;//number of constrains/inequalities
    //inequality formula: constrains on segment and junction CNs
    for (int i=1;i<=numSegments;i++) {
        //Σp + Σ2*l + e_i >= c_i and 
        //Σp + Σ2*l - e_i <= c_i
        CoinPackedVector constrain1, constrain2;
        for (int j=0;j<numPat;j++) {//Σp
            if (patterns[j][0]<=i && i<=patterns[j][1]) {
                string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                constrain1.insert(variableIdx[key], 1);
                constrain2.insert(variableIdx[key], 1);
            }
        }
        for (int j=0;j<numLoop;j++) {//Σ2*l
            if (loops[j][0]<=i && i<=loops[j][1]) {
                string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                constrain1.insert(variableIdx[key], 2);
                constrain2.insert(variableIdx[key], 2);
            }
        }
        constrain1.insert(numElements+idx/2, 1);//e_i
        constrainLowerBound[idx] = (*segs)[i-1]->getWeight()->getCorrectedCoverage();
        constrainUpperBound[idx] = si->getInfinity();
        idx++;
        matrix->appendRow(constrain1);

        constrain2.insert(numElements+idx/2, -1);//-e_i
        constrainLowerBound[idx] = -1*si->getInfinity();
        constrainUpperBound[idx] = (*segs)[i-1]->getWeight()->getCorrectedCoverage();
        idx++;
        matrix->appendRow(constrain2);

        //Σp + Σ2*l + e_ij >= c_ij where j=i+1 and
        //Σp + Σ2*l - e_ij <= c_ij where j=i+1
        CoinPackedVector constrain3, constrain4;
        for (int j=0;j<numPat;j++) {//Σp
            if (patterns[j][0]<=i && i<patterns[j][1]) {
                string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                constrain3.insert(variableIdx[key], 1);
                constrain4.insert(variableIdx[key], 1);
            }
        }
        for (int j=0;j<numLoop;j++) {//Σ2*l
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

        //Σl + Σ(p1+p2)/2 + Σ(p+l)/2 + Σ(l1+l2)/2 + e_ii >= c_ii and
        //Σl + Σ(p1+p2)/2 + Σ(p+l)/2 + Σ(l1+l2)/2 - e_ii <= c_ii
        double *coef = new double[numElements];
        memset(coef, 0, sizeof(double)*numElements);
        CoinPackedVector constrain5, constrain6;
        for (int j=0;j<numLoop;j++) {//Σl
            if (loops[j][0] == i || loops[j][1] == i) {
                string key = "l:"+to_string(loops[j][0])+","+to_string(loops[j][1]);
                coef[variableIdx[key]] += 1;
            }
        }
        for (int j=0;j<numPat;j++) {//Σ(p1+p2)/2
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
        for (int j=0;j<numPat;j++) {//Σ(p+l)/2
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
        for (int j=0;j<numLoop;j++) {//Σ(l1+l2)/2 (loop2 is inserted into the middle of loop1)
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
    //inequality formula: constrains on patterns and loops
    //Σp(c,d)-p(a,b)>=0 where p(a,b) is a sub-pattern of p(c,d)
    for (int i=0;i<numPat;i++) {
        CoinPackedVector constrain7;
        bool flag = false;
        for (int j=0;j<numPat;j++) {
            if ((patterns[i][0] == patterns[j][0]) ||
                (patterns[i][1] == patterns[j][1])) {
                int diff1 = patterns[i][0]-patterns[i][1], diff2 = patterns[j][0]-patterns[j][1];
                if (abs(diff1)<abs(diff2)) {
                    flag = true;
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain7.insert(variableIdx[key], 1);
                }
            }
        }
        if (flag) {
            string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
            constrain7.insert(variableIdx[key], -1);
            constrainLowerBound[idx] = 0;
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain7);
        }
    }
    //Σp(x1,y1)+Σl(x2,y2)-l(a,b)>=0 where l(a,b) pattern is a sub-pattern of p and l
    for (int i=0;i<numLoop;i++) {
        CoinPackedVector constrain8;
        bool flag = false;
        for (int j=0;j<numPat;j++) {
            if ((patterns[j][0] == loops[i][0])||
                (patterns[j][1] == loops[i][1])) {
                int diff1 = loops[i][0]-loops[i][1], diff2 = patterns[j][0]-patterns[j][1];
                if (abs(diff1)<abs(diff2)) {
                    flag = true;
                    string key = "p:"+to_string(patterns[j][0])+","+to_string(patterns[j][1]);
                    constrain8.insert(variableIdx[key], 1);
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
                    constrain8.insert(variableIdx[key], 1);
                }
            }
        }
        if (flag) {
            string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
            constrain8.insert(variableIdx[key], -1);
            constrainLowerBound[idx] = 0;
            constrainUpperBound[idx] = si->getInfinity();
            idx++;
            matrix->appendRow(constrain8);
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
    for (int i=0;i<numLoop;i++) {//l
        string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
        variableLowerBound[variableIdx[key]] = 0;
        variableUpperBound[variableIdx[key]] = maxCN;
    }
    for (int i=0;i<numEpsilons;i++) {//ε (or e)
        variableLowerBound[numElements+i] = 0;
        variableUpperBound[numElements+i] = si->getInfinity();
    }
    cout<<"Variable constrains done"<<endl;
    //weights for all variables
    for (int i=0;i<numVariables;i++) {
        if (i<numElements)
            objective[i] = 0;
        else   
            objective[i] = 1;
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

//Hierholzer’s Algorithm
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