#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <queue>
#include <iomanip>
#include <assert.h>

#include "Graph.hpp"
#include "Exceptions.hpp"

using namespace std;

Graph::Graph() {
    mPurity = -1;
    mAvgPloidy = -1;
    mAvgTumorPloidy = -1;

    mSegments = new vector<Segment *>();
    mJunctions = new vector<Junction *>();
    mSources = new vector<Segment *>();
    mSinks = new vector<Segment *>();
}

Graph::Graph(const char *aFilename) {
    mPurity = -1;
    mAvgPloidy = 0;
    mAvgTumorPloidy = -1;
    mAvgCoverageRaw = -1;
    mAvgVirusDP = -1;

    mSegments = new vector<Segment *>();
    mJunctions = new vector<Junction *>();
    mSources = new vector<Segment *>();
    mSinks = new vector<Segment *>();
    mAvgCoverages = new vector<double>();
    this->readGraph(aFilename);
}

Graph::~Graph() { ; }

string Graph::getSampleName() { return mSampleName; }

int Graph::getJunctionIndexByEdge(Edge *aEdge) {
    int idx = 0;
    while (true) {
        if ((*mJunctions)[idx] == aEdge->getJunction())
            return idx;
        idx++;
    }
    throw JunctionDoesNotExistException(aEdge);
}

int Graph::getExpectedPloidy() { return mExpectedPloidy; }
int Graph::getExpectedPloidy2(int i) {
    int * r = new int [8,2,2,2];
    return r[i];
}

double Graph::getPurity() { return mPurity; }

double Graph::getAvgPloidy() { return mAvgPloidy; }

double Graph::getAvgTumorPloidy() { return mAvgTumorPloidy; }

double Graph::getAvgCoverage() { return mAvgCoverage; }

int Graph::getVirusSegStart() { return mVirusSegStart;}

double Graph::getAvgRawCoverage() { return mAvgCoverageRaw; }

double Graph::getAvgCoverageJunc() { return mAvgCoverageJunc; }

double Graph::getAvgRawCoverageJunc() { return mAvgCoverageRawJunc; }

double Graph::getHaploidDepth() { return mHaploidDepth; }
double Graph::getRatio() {return mRatio;}

void Graph::setPurity(double aPurity) { mPurity = aPurity; }

void Graph::setAvgPloidy(double aAvgPloidy) { mAvgPloidy = aAvgPloidy; }

void Graph::setAvgTumorPloidy(double aAvgTumorPloidy) { mAvgTumorPloidy = aAvgTumorPloidy; }

void Graph::setAvgCoverage(double aAvgCoverage) { mAvgCoverage = aAvgCoverage; }

void Graph::setAvgRawCoverage(double aAvgRawCoverage) { mAvgCoverageRaw = aAvgRawCoverage; }

Segment *Graph::getFirstSource() { return (*mSources)[0]; }

Segment *Graph::getFirstSink() { return (*mSinks)[0]; }

vector<Segment *> *Graph::getSegments() { return mSegments; }

vector<Junction *> *Graph::getJunctions() { return mJunctions; }

// functionality
void Graph::readGraph(const char *aFilename) {
    ifstream graphFile(aFilename);
    if (!graphFile) {
        cerr << "Cannot open file " << aFilename << endl;
        exit(1);
    }

    cout << "Reading graph..." << endl;
    char line[8192];
    char *token;
    vector<int> sourceIds, sinkIds;
//    int sourceId, sinkId;
    while (!graphFile.eof()) {
        graphFile.getline(line, 8192);

        char *line_p = line;
        while (*line_p != '\0') {
            if (*line_p != '\t' && *line_p != ' ') {
                break;
            }
            line_p++;
        }
        if (*line_p == '#') {
            continue;
        }
        // cout << "\"" << line << "\"" << endl;

        token = strtok(line, " \t");
        if (token == NULL) {
            continue;
        }
        if (strcmp(token, "SAMPLE_NAME") == 0) {
            mSampleName = string(strtok(NULL, " "));
        } else if (strcmp(token, "AVG_CHR_SEG_DP") == 0) {
            token = strtok(NULL, " ");
            token = strtok(token, ",");
            while (token != NULL) {
                mAvgCoverages->push_back(atof(token));
                token = strtok(NULL, ",");
            }
        } else if (strcmp(token, "AVG_WHOLE_HOST_DP") == 0) {
            mAvgCoverageRaw = atof(strtok(NULL, " "));
        } else if (strcmp(token, "AVG_VIRUS_SEG_DP") == 0) {
            mAvgVirusDP = atof(strtok(NULL, " "));
        } else if (strcmp(token, "VIRUS_START") == 0) {
            mVirusSegStart = atoi(strtok(NULL, " "));
        } else if (strcmp(token, "AVG_JUNC_DP") == 0) {
            mAvgCoverageJunc = atof(strtok(NULL, " "));
            mAvgCoverageRawJunc = mAvgCoverageJunc;
        } else if (strcmp(token, "PURITY") == 0) {
            mPurity = atof(strtok(NULL, " "));
        } else if (strcmp(token, "AVG_TUMOR_PLOIDY") == 0) {
            mAvgTumorPloidy = atof(strtok(NULL, " "));
        } else if (strcmp(token, "AVG_PLOIDY") == 0) {
            mAvgPloidy = atof(strtok(NULL, " "));
        } else if (strcmp(token, "PLOIDY") == 0) {
            token = strtok(NULL, " ");
            mPloidy = string(token);
            mExpectedPloidy = atoi(strtok(token, "m"));
        } else if (strcmp(token, "SOURCE") == 0) {
            token = strtok(NULL, " ");
            token = strtok(token, ",");
            while (token != NULL) {
                sourceIds.push_back(atoi(token));
                token = strtok(NULL, ",");
            }
        } else if (strcmp(token, "SINK") == 0) {
            token = strtok(NULL, " ");
            token = strtok(token, ",");
            while (token != NULL) {
                sinkIds.push_back(atoi(token));
                token = strtok(NULL, ",");
            }
        } else if (strcmp(token, "SEG") == 0) {
            char *node = strtok(NULL, " ");
            double segCoverage = max(atof(strtok(NULL, " ")), 0.0);
            double segCopy = atof(strtok(NULL, " "));

            strtok(node, ":");
            int segId = atoi(strtok(NULL, ":"));
            string chrom = strtok(NULL, ":");
            int start = atoi(strtok(NULL, ":"));
            int end = atoi(strtok(NULL, ":"));
            double segCred = 1.0;
            int chrId;
            for (int i=0; i<sourceIds.size(); i++) {
                if (sourceIds[i]<=segId && segId<=sinkIds[i])
                    chrId = i;
            }
            this->addSegment(segId, chrId, chrom, start, end, segCoverage, segCred, segCopy);
        } else if (strcmp(token, "JUNC") == 0) {
            char *sourceNode = strtok(NULL, " ");
            char *targetNode = strtok(NULL, " ");
            double junCoverage = atof(strtok(NULL, " "));
            double junCopy = atof(strtok(NULL, " "));
            token = strtok(NULL, " ");
            // cout << token << endl;
            bool isInferred = (token[0] == 'I') ? true : false;
            token = strtok(NULL, " ");
            // cout << token << endl;
            bool isBounded = (token[0] == 'B') ? true : false;
            // cout << (isBounded ? "yes" : "no") << endl;
            if (junCoverage <= 0 && junCopy <= 0) continue;

            strtok(sourceNode, ":");
            int sourceId = atoi(strtok(NULL, ":"));
            char sourceDir = strtok(NULL, ":")[0];

            strtok(targetNode, ":");
            int targetId = atoi(strtok(NULL, ":"));
            char targetDir = strtok(NULL, ":")[0];

            double junCred = 1.0;

            this->addJunction(sourceId, sourceDir, targetId, targetDir, junCoverage, junCred, junCopy, isInferred,
                              isBounded, false);
        }
    }
    assert(sourceIds.size() == sinkIds.size());

    for (int i = 0 ;i < sourceIds.size(); i++) {
        mSources->push_back(this->getSegmentById(sourceIds[i]));
        mSinks->push_back(this->getSegmentById(sinkIds[i]));
    }
    mInferredBegin = mJunctions->end();
//
//    mSource = this->getSegmentById(sourceId);
//    mSink = this->getSegmentById(sinkId);
}

void Graph::writeGraph(const char *aFilename) {
    ofstream fout(aFilename);
    fout << "SAMPLE_NAME " << "TEST" << endl
         << "AVG_SEG_DP " << mAvgCoverage << endl
         << "AVG_JUNC_DP " << mAvgCoverageJunc << endl
         << "PURITY " << mPurity << endl
         << "AVG_PLOIDY " << mAvgPloidy << endl
         << "PLOIDY " << mPloidy << endl
         << "SOURCE " << getSourcesIds() << endl
         << "SINK " << getSinksIds() << endl;
    cout << "write seg" << endl;
    for (Segment *seg : *mSegments) {
        fout << "SEG " << "H:" << seg->getId() << ":" << seg->getChrom()
             << ":" << seg->getStart() << ":" << seg->getEnd() << " "
             << seg->getWeight()->getCoverage() << " "
             << seg->getWeight()->getCopyNum() << " "
             << (seg->hasLowerBoundLimit() ? "B" : "U") << endl;
    }
    for (Junction *junc : *mJunctions) {
        Edge *e = junc->getEdgeA();
        fout << "JUNC " << "H:" << e->getSource()->getId() << ":" << e->getSource()->getDir() << " "
             << "H:" << e->getTarget()->getId() << ":" << e->getTarget()->getDir() << " "
             << junc->getWeight()->getCoverage() << " "
             << junc->getWeight()->getCopyNum() << " "
             << (junc->isInferred() ? "I" : "U") << " " << (junc->hasLowerBoundLimit() ? "B" : "U") << endl;
    }
    fout.close();
}

void Graph::writePartialGraph(vector<Segment *> *segs, const char *aFilename) {
    ofstream fout(aFilename);
    fout << "SAMPLE_NAME " << mSampleName << endl
         << "AVG_SEG_DP " << mAvgCoverage << endl
         << "AVG_JUNC_DP " << mAvgCoverageJunc << endl
         << "PURITY " << mPurity << endl
         << "AVG_PLOIDY " << mAvgPloidy << endl
         << "PLOIDY " << mPloidy << endl
         << "SOURCE " << getSourcesIds() << endl
         << "SINK " << getSinksIds() << endl;
    cout << "write seg" << endl;
    for (Segment *seg : *segs) {
        fout << "SEG " << "H:" << seg->getId() << ":" << seg->getChrom()
             << ":" << seg->getStart() << ":" << seg->getEnd() << " "
             << seg->getWeight()->getCoverage() << " "
             << seg->getWeight()->getCopyNum() << " "
             << (seg->hasLowerBoundLimit() ? "B" : "U") << endl;
    }
    for (Junction *junc : *mJunctions) {
        bool sourceFlag = false;
        bool targetFlag = false;
        for (Segment *seg : *segs) {
            if (seg->getId() == junc->getSource()->getId()) sourceFlag = true;
            if (seg->getId() == junc->getTarget()->getId()) targetFlag = true;
        }
        if (sourceFlag && targetFlag) {
            Edge *e = junc->getEdgeA();
            fout << "JUNC " << "H:" << e->getSource()->getId() << ":" << e->getSource()->getDir() << " "
                 << "H:" << e->getTarget()->getId() << ":" << e->getTarget()->getDir() << " "
                 << junc->getWeight()->getCoverage() << " "
                 << junc->getWeight()->getCopyNum() << " "
                 << (junc->isInferred() ? "I" : "U") << " " << (junc->hasLowerBoundLimit() ? "B" : "U") << endl;
        }
    }
    fout.close();
}


void Graph::checkOrphan() {
    for (Segment *seg : *mSegments) {
        seg->checkOrphan();
    }
}

void Graph::calculateHapDepth() {
    /* haploid depth = average depth / average ploidy
     * average ploidy = purity * average tumor ploidy + (1 - purity) * 2
     * purity * segment copy * haploytye depth + (1 - purity) * 2 * haploid depth = segment depth
    */

    if (mAvgPloidy < 0) {
        if (mAvgTumorPloidy < 0) {
            cerr
                    << "input error: there is no ploidy information provided. There must be at least one of \"AVG_PLOIDY\" and \"AVG_TUMOR_PLOIDY\"."
                    << endl;
            exit(1);
        } else {
            if (mPurity < 0) {
                cerr << "input error: no purity information provided." << endl;
                exit(1);
            } else {
                mAvgPloidy = mPurity * mAvgTumorPloidy + (1 - mPurity) * 2;
            }
        }
    } else {
        if (mAvgTumorPloidy >= 0) {
            cout << "WARN: both AVG_PLOIDY and AVG_TUMOR_PLOIDY are provided, ";

            if (mPurity < 0) {
                cout << "WARN: no purity information provided, use the given AVG_PLOIDY" << endl;
            } else {
                double ratio = 1 - (mPurity * mAvgTumorPloidy)/((mPurity * mAvgTumorPloidy) + (1 - mPurity) * 2 );
                double avgPloidy = mPurity * mAvgTumorPloidy + (1 - mPurity) * 2;
                mRatio = ratio;
                if (abs(mAvgPloidy - avgPloidy) <= 0.1) {
                    cout
                            << "calculated AVG_PLOIDY using AVG_TUMOR_PLOIDY is close enough to the given AVG_PLOIDY, use the given AVG_PLOIDY"
                            << endl;
                } else {
                    mAvgPloidy = avgPloidy;
                    cout
                            << "calculated AVG_PLOIDY using AVG_TUMOR_PLOIDY is distinguishable from the given AVG_PLOIDY, use the calculated AVG_PLOIDY"
                            << endl;
                }
            }
        } else {
            cout << "WARN: only AVG_PLOIDY is given, use that" << endl;
        }
    }

    // mHaploidDepth = mAvgCoverageRaw / mAvgPloidy;
    // mHaploidDepthJunc = mAvgCoverageRawJunc / mAvgPloidy;
    mHaploidDepth = mAvgCoverageRaw * mPurity / mAvgPloidy;
    mHaploidDepthJunc = mHaploidDepth;
    cout << "Average ploidy: " << mAvgPloidy << endl
         << "Haploid depth: " << mHaploidDepth << endl;
    mAvgCoverage = mAvgPloidy * mHaploidDepth;
    mAvgCoverageJunc = mAvgPloidy * mHaploidDepthJunc;

}

void Graph::calculateCopyNum() {
    double ratio = getRatio();
    double hDP = getHaploidDepth();
    for (Segment *seg : *mSegments) {
//        TODO 判断virus seg
        double segCopy;
        if(seg->getId() >= mVirusSegStart) {
            segCopy = seg->getWeight()->getCoverage() / mAvgCoverageRaw * 2;
        } else {
            double depthT = seg->getWeight()->getCoverage() - mAvgCoverageRaw * ratio;
            seg->getWeight()->setCorrectedCoverage(depthT);
            segCopy = depthT/hDP;
        }
        // seg->getWeight()->setAdjustedCoverage(max(segAdjustedCoverage, 0.0));
        // seg->getWeight()->setCoverage(seg->getWeight()->getAdjustedCoverage());
        seg->getWeight()->setCopyNum(max(segCopy, 0.0));
    }

    for (Junction *junc : *mJunctions) {
//        double juncCopy = junc->getWeight()->getCoverage() / mHaploidDepthJunc;
        double juncCopy;
        if (junc->isInferred()) {
            cout << mHaploidDepth << endl;
//            junc->getWeight()->setCoverage(hDP);
        }
//        juncCopy = junc->getWeight()->getCoverage() / mHaploidDepth;
        double depthT = junc->getWeight()->getCoverage() - mAvgCoverageRaw * ratio;
        juncCopy = depthT/hDP;
        junc->getWeight()->setCorrectedCoverage(depthT);
        // junc->getWeight()->setAdjustedCoverage(max(juncAdjustedCoverage, 0.0));
        // junc->getWeight()->setCoverage(junc->getWeight()->getAdjustedCoverage());
        junc->getWeight()->setCopyNum(max(juncCopy, 0.0));
    }
}

void Graph::restoreCopy() {
    for (Segment *seg : *mSegments) {
        seg->restoreCopy();
    }
    for (Junction *junc : *mJunctions) {
        junc->restoreCopy();
    }
}

void Graph::backupCopy() {
    for (Segment *seg : *mSegments) {
        seg->backupCopy();
    }
    for (Junction *junc : *mJunctions) {
        junc->backupCopy();
    }
}

void Graph::resetVertexVisitFlag() {
    for (Segment *seg : *mSegments) {
        seg->getPositiveVertex()->resetVisited();
        seg->getNegativeVertex()->resetVisited();
        // seg->getPositiveVertex()->clearShortestPrevVertex();
        // seg->getNegativeVertex()->clearShortestPrevVertex();
        // seg->getPositiveVertex()->clearShortestPrevEdge();
        // seg->getNegativeVertex()->clearShortestPrevEdge();
    }
}

void Graph::resetJunctionVisitFlag() {
    for (Junction *junc : *mJunctions) {
        junc->getEdgeA()->resetVisited();
        junc->getEdgeB()->resetVisited();
    }
}

void Graph::resetShortestPrevEdge() {
    for (Segment *seg : *mSegments) {
        seg->getPositiveVertex()->setShortestPrevEdge(NULL);
        seg->getNegativeVertex()->setShortestPrevEdge(NULL);
    }
}

void Graph::checkLowerBound() {
    this->checkOrphan();
    for (Segment *seg : *mSegments) {
//        if(isTarget) {
//            seg->setHasLowerBoundLimit();
//            continue;
//        }
        // seg->checkLowerBound();
        if (seg->isOrphan()) {
            if (seg->getWeight()->getCoverage() <= 0.25 * mAvgCoverage) {
                seg->resetHasLowerBoundLimit();
                // } else {
                //     seg->resetHasLowerBoundLimit();
            }
        } else {
            seg->setHasLowerBoundLimit();
        }
    }
    for (Junction *junc : *mJunctions) {
        if (junc->getWeight()->getCoverage() > 0.25 * mAvgCoverageJunc) {
            junc->setHasLowerBoundLimit();
        } else {
            junc->resetHasLowerBoundLimit();
        }
    }
    // for (Junction * junc : *mJunctions) {
    //     junc->checkLowerBound();
    // }
}

bool Graph::isCopyExhaustive() {
    for (Segment *seg : *mSegments) {
        if (seg->hasCopy()) {
            return false;
        }
    }
    return true;
}

bool Graph::doesJunctionExist(Junction *aJunction) {
    vector<string> aJuncInfo = aJunction->getInfo();
    for (Junction *junc : *mJunctions) {
        vector<string> juncInfo = junc->getInfo();
        if ((juncInfo[0] == aJuncInfo[0] && juncInfo[1] == aJuncInfo[1]) ||
            (juncInfo[0] == aJuncInfo[1] && juncInfo[1] == aJuncInfo[0])) {
            return true;
        }
    }
    return false;
}

Segment *Graph::getSegmentById(int aSegId) {
    for (Segment *seg : *mSegments) {
        if (seg->getId() == aSegId) {
            return seg;
        }
    }
    throw SegmentDoesNotExistException(aSegId);
}

Segment *Graph::getSegmentByChromStart(string aChrom, int aStart) {
    for (Segment *seg : *mSegments) {
        if (seg->getChrom() == aChrom && seg->getStart() == aStart) {
            return seg;
        }
    }
    throw SegmentDoesNotExistException(aStart);
}

Segment *Graph::getSegmentByChromEnd(string aChrom, int aEnd) {
    for (Segment *seg : *mSegments) {
        if (seg->getChrom() == aChrom && seg->getEnd() == aEnd) {
            return seg;
        }
    }
    throw SegmentDoesNotExistException(aEnd);
}

Segment *Graph::addSegment(int aId, int chr, string aChrom, int aStart, int aEnd, double aCoverage, double aCredibility, double aCopy) {
    Segment *seg = new Segment(aId, chr, aChrom, aStart, aEnd, aCoverage, aCredibility, aCopy);
    mSegments->push_back(seg);
    return seg;
}

Junction *Graph::addJunction(Vertex *aSource, Vertex *aTarget, double aCoverage, double aCredibility, double aCopy,
                             bool aInferred, bool aIsBounded, bool aIsSourceSinkJunction) {
    Segment *sourceSeg = aSource->getSegment();
    Segment *targetSeg = aTarget->getSegment();

    if (!sourceSeg->hasLowerBoundLimit() || !targetSeg->hasLowerBoundLimit()) {
        // cout << (sourceSeg->hasLowerBoundLimit() ? "s yes" : "s no") << endl;
        // cout << (targetSeg->hasLowerBoundLimit() ? "t yes" : "t no") << endl;
        return NULL;
    }

    Junction *junc = new Junction(sourceSeg, targetSeg, aSource->getDir(), aTarget->getDir(), aCoverage, aCredibility,
                                  aCopy, aInferred, aIsBounded, aIsSourceSinkJunction);
    if (this->doesJunctionExist(junc)) {
        throw DuplicateJunctionException(junc);
        delete junc;
    }
    junc->insertEdgesToVertices();
    mJunctions->push_back(junc);
    // if (aIsSourceSinkJunction) {
    //     mJunctions->push_back(junc);
    // } else {
    //     if (mJunctions->size() > 0) {
    //         mJunctions->insert(mJunctions->end() - 1, junc);
    //     } else {
    //         mJunctions->push_back(junc);
    //     }
    // }
    return junc;
    // TODO
    // if mInferred, append to a vector storing inferred junctions
}

Junction *Graph::addJunction(int aSourceId, char aSourceDir, int aTargetId, char aTargetDir, double aCoverage,
                             double aCredibility, double aCopy, bool aInferred, bool aIsBounded,
                             bool aIsSourceSinkJunction) {
    Segment *sourceSeg = this->getSegmentById(aSourceId);
    Segment *targetSeg = this->getSegmentById(aTargetId);

    if (!sourceSeg->hasLowerBoundLimit() || !targetSeg->hasLowerBoundLimit()) {
        // cout << "dd" << endl;
        return NULL;
    }

    Junction *junc = new Junction(sourceSeg, targetSeg, aSourceDir, aTargetDir, aCoverage, aCredibility, aCopy,
                                  aInferred, aIsBounded, aIsSourceSinkJunction);
    if (this->doesJunctionExist(junc)) {
        // throw DuplicateJunctionException(junc);
        return junc;
    }
    junc->insertEdgesToVertices();
    mJunctions->push_back(junc);
    // if (aIsSourceSinkJunction) {
    //     mJunctions->push_back(junc);
    // } else {
    //     if (mJunctions->size() > 0) {
    //         mJunctions->insert(mJunctions->end() - 1, junc);
    //     } else {
    //         mJunctions->push_back(junc);
    //     }
    // }
    return junc;
    // TODO
    // if mInferred, append to a vector storing inferred junctions
}

Vertex *Graph::getNextVertexById(Vertex *aSourceVertex) {
    Segment *nextSegment;
    if (aSourceVertex->getDir() == '+') {
        nextSegment = this->getSegmentById(aSourceVertex->getId() + 1);
        // if (nextSegment == NULL) {
        //     throw NotReachableException(aSourceVertex);
        // }
        return nextSegment->getPositiveVertex();
    } else {
        nextSegment = this->getSegmentById(aSourceVertex->getId() - 1);
        // if (nextSegment == NULL) {
        //     throw NotReachableException(aSourceVertex);
        // }
        return nextSegment->getNegativeVertex();
    }
    // return nextSegment->getNegativeVertex();
}

Vertex *Graph::getPrevVertexById(Vertex *aTargetVertex) {
    Segment *prevSegment;
    if (aTargetVertex->getDir() == '+') {
        prevSegment = this->getSegmentById(aTargetVertex->getId() - 1);
        // if (prevSegment == NULL) {
        //     throw NotReachableException(aTargetVertex);
        // }
        return prevSegment->getPositiveVertex();
    } else {
        prevSegment = this->getSegmentById(aTargetVertex->getId() + 1);
        // if (prevSegment == NULL) {
        //     throw NotReachableException(aTargetVertex);
        // }
        return prevSegment->getNegativeVertex();
    }
    // return prevSegment->getNegativeVertex();
}

int Graph::BFS(Vertex *aStartVertex, Vertex *aTargetVertex) {
    queue<Vertex *> vertexQueue;
    vertexQueue.push(aStartVertex);

    while (!vertexQueue.empty()) {
        Vertex *currentVertex = vertexQueue.front();
        vertexQueue.pop();

        for (Edge *e : *(currentVertex->getEdgesAsSource())) {
            // if (e->hasCopy()) {
            Vertex *nextVertex = e->getTarget();
            if (!nextVertex->isVisited()) {
                nextVertex->setVisited();
                nextVertex->setShortestPrevEdge(e);
                vertexQueue.push(nextVertex);
            }
            // }
        }
    }
    this->resetVertexVisitFlag();

    Edge *prevEdge;
    Vertex *prevVertex;
    Vertex *currentVertex = aTargetVertex;
    int found;
    while (true) {
        prevEdge = currentVertex->getShortestPrevEdge();
        if (prevEdge == NULL) {
            found = -1;
            break;
            // return -1;  // no path from aStartVertex to aTargetVertex
        }

        prevVertex = prevEdge->getSource();
        if (prevVertex == aStartVertex) {
            found = 0;
            break;
            // return 0;
        }
        currentVertex = prevVertex;
    }
    this->resetShortestPrevEdge();
    return found;
}

// int Graph::findShortestPath(Vertex * aStartVertex, Vertex * aTargetVertex) {
//     Edge * prevEdge;
//     Vertex * prevVertex;
//     Vertex * currentVertex = aTargetVertex;
//     while (true) {
//         prevEdge = currentVertex->getShortestPrevEdge();
//         if (prevEdge == NULL) {
//             return -1;  // no path from aStartVertex to aTargetVertex
//         }

//         prevVertex = prevEdge->getFirstSource();
//         if (prevVertex == aStartVertex) {
//             return 0;
//         }
//         currentVertex = prevVertex;
//     }
// }

void Graph::print() {
    cout << "``````````````````````````````````````````````````````````````````````````````````````````" << endl;
    cout << "Ploidy: " << mExpectedPloidy << endl;
    cout << "Avg ploidy: " << mAvgPloidy << endl;
    cout << "Avg coverage: " << mAvgCoverage << endl;
    cout << "Haploid coverage: " << mHaploidDepth << endl;
    cout << "Haploid junc coverage: " << mHaploidDepthJunc << endl;
    cout << "Source: " << getSourcesIds() << endl;
    cout << "Sink: " << getSinksIds() << endl;
    cout << "Segments: " << mSegments->size() << endl;
    cout << "Junctions: " << mJunctions->size() << endl;
    cout << "``````````````````````````````````````````````````````````````````````````````````````````" << endl;
    for (Segment *seg : *mSegments) {
        cout << fixed << setprecision(4)
             << "SEG" << "\t"
             << seg->getId() << "\t"
             << seg->getChrom() << "\t"
             << seg->getStart() << "\t"
             << seg->getEnd() << "\t"
             << seg->getWeight()->getCoverage() << "\t"
             // << seg->getWeight()->getOriginalCoverage() << "\t"
             // << seg->getWeight()->getAdjustedCoverage() << "\t"
             << seg->getWeight()->getCoverage() / mHaploidDepth << "\t"
             << "\033[1;31m" << seg->getWeight()->getCopyNum() << "\033[0m" << "\t"
             << seg->getCredibility() << "\t"
             << (seg->isOrphan() ? "OO" : "ONO") << "\t"
             << (seg->getPositiveVertex()->isOrphan() ? "OO" : "ONO") << "\t"
             << (seg->getNegativeVertex()->isOrphan() ? "OO" : "ONO") << "\t"
             << (seg->isOrphan(false) ? "O" : "NO") << "\t"
             << (seg->getPositiveVertex()->isOrphan(false) ? "O" : "NO") << "\t"
             << (seg->getNegativeVertex()->isOrphan(false) ? "O" : "NO") << "\t"
             << (seg->isDeadEnd() ? "D" : "ND") << "\t"
             << (seg->hasLowerBoundLimit() ? "LB" : "NLB") << "\t"
             << ((seg->getPositiveVertex()->getShortestPrevEdge() == NULL) ? "NULL"
                                                                           : seg->getPositiveVertex()->getShortestPrevEdge()->getSource()->getInfo())
             << "\t"
             << ((seg->getNegativeVertex()->getShortestPrevEdge() == NULL) ? "NULL"
                                                                           : seg->getNegativeVertex()->getShortestPrevEdge()->getSource()->getInfo())
             << endl;
    }

    int c = 0;
    for (Junction *junc : *mJunctions) {
        vector<string> info = junc->getInfo();
        cout << left << fixed << setprecision(4)
             << "JUNC" << "\t"
             << c << "\t"
             << info[0] << "\t" << info[1] << "\t"
             << junc->getWeight()->getCoverage() << "\t"
             // << junc->getWeight()->getOriginalCoverage() << "\t"
             // << junc->getWeight()->getAdjustedCoverage() << "\t"
             << junc->getWeight()->getCoverage() / mHaploidDepth << "\t"
             << "\033[1;31m" << junc->getWeight()->getCopyNum() << "\033[0m" << "\t"
             << junc->getCredibility() << "\t"
             << (junc->isInferred() ? "I" : "NI") << "\t"
             << (junc->hasLowerBoundLimit() ? "LB" : "NLB") << endl;
        c++;
    }
    cout << "``````````````````````````````````````````````````````````````````````````````````````````" << endl;
}

vector<Segment *> *Graph::getMSources() const {
    return mSources;
}

vector<Segment *> *Graph::getMSinks() const {
    return mSinks;
}

string Graph::getSourcesIds() {
    string r;
    for (auto s : *mSources) {
        r += std::to_string(s->getId());
        r += ',';
    }
    return r;
}

string Graph::getSinksIds() {
    string r="";
    for (auto s : *mSinks) {
        r += std::to_string(s->getId());
        r += ',';
    }
    return r;
}