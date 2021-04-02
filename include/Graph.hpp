#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>

#include "Segment.hpp"
#include "Junction.hpp"

using namespace std;

class Segment;
class Junction;

class Graph {
    protected:
        string mSampleName;
        
        int mExpectedPloidy;     // expected ploidy

        double mHaploidDepth;
        double mHaploidDepthJunc;
        double mPurity;   // tumor purity
        double mAvgPloidy;   // the ploidy of whole sequencing region, including tumor if possible
        double mAvgTumorPloidy;   // the ploidy of tumor
        double mAvgCoverage;
        double mAvgCoverageRaw;
        double mAvgCoverageJunc;
        double mAvgCoverageRawJunc;

        string mPloidy;

        Segment * mSource;
        Segment * mSink;

        vector<Segment *> * mSegments;
        vector<Junction *> * mJunctions;
        vector<Junction *>::iterator mInferredBegin;

    public:
        Graph();
        Graph(const char * aFilename);
        ~Graph();

        string getSampleName();

        int getJunctionIndexByEdge(Edge * aEdge);

        int getExpectedPloidy();
        double getPurity();
        double getAvgPloidy();
        double getAvgTumorPloidy();
        double getAvgCoverage();
        double getAvgRawCoverage();
        double getAvgCoverageJunc();
        double getAvgRawCoverageJunc();
        double getHaploidDepth();
        void setPurity(double aPurity);
        void setAvgPloidy(double aAvgPloidy);
        void setAvgTumorPloidy(double aAvgTumorPloidy);
        void setAvgCoverage(double aAvgCoverage);
        void setAvgRawCoverage(double aAvgCoverageRaw);

        Segment * getSource();
        Segment * getSink();
        
        vector<Segment *> * getSegments();
        vector<Junction *> * getJunctions();

        /* functionality */
        
        bool isCopyExhaustive();
        bool doesJunctionExist(Junction * aJunction);

        void readGraph(const char * aFilename);
        void writeGraph(const char * aFilename);
        void checkOrphan();
        void calculateHapDepth();
        void calculateCopyNum();
        void restoreCopy();
        void backupCopy();
        void resetVertexVisitFlag();
        void resetJunctionVisitFlag();
        void resetShortestPrevEdge();
        void checkLowerBound();

        Segment * getSegmentById(int aSegId);
        Segment * getSegmentByChromStart(string aChrom, int aStart);
        Segment * getSegmentByChromEnd(string aChrom, int aEnd);
        Segment * addSegment(int aId, string aChrom, int aStart, int aEnd, double aCoverage, double aCredibility, double aCopy);
        Junction * addJunction(Vertex * aSource, Vertex * aTarget, double aCoverage, double aCredibility, double aCopy, bool aInferred, bool aIsBounded, bool aIsSourceSinkJunction);
        Junction * addJunction(int aSourceId, char aSourceDir, int aTargetId, char aTargetDir, double aCoverage, double aCredibility, double aCopy, bool aInferred, bool aIsBounded, bool aIsSourceSinkJunction);
        Vertex * getPrevVertexById(Vertex * aTargetVertex);
        Vertex * getNextVertexById(Vertex * aSourceVertex);

        int BFS(Vertex * aStartVertex, Vertex * aTargetVertex);
        // int findShortestPath(Vertex * aStartVertex, Vertex * aTargetVertex);

        void print();
};

#endif
