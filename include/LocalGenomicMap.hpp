#ifndef _LOCAL_GENOMIC_MAP_H_
#define _LOCAL_GENOMIC_MAP_H_

#include "Graph.hpp"
#include "JunctionDB.hpp"

class Graph;
class JunctionDB;

typedef vector<Vertex *> VertexPath;
typedef vector<Edge *> EdgePath;

typedef struct TrieNode {
    bool mIsEnd;
    vector<TrieNode *> mChildren;
    vector<int> mIndexes;
} trie_node_t;

class LocalGenomicMap {
    protected:
        Graph * mGraph;
        
        vector<VertexPath *> * mRawPathVertices;
        vector<EdgePath> * mRawPathEdges;

        vector<VertexPath *> * mCircuits;
        vector<VertexPath *> * mHaploids;

        vector<VertexPath *> * mLongFrags;

        bool usingLong;
        bool usingHic;

        double ** hicMatrix;
        double ** decreaseMatrix;
    public:
        LocalGenomicMap(Graph * aGraph);
        ~LocalGenomicMap();
        
        /* getter and setter */
        Graph * getGraph();
        void setGraph(Graph * aGraph);
        bool isUsingLong();
        bool isUsingHic();
        void setUsingLong(bool v);
        void setUsingHic(bool );
        
        /* Long fragments for TGS */
        vector<VertexPath *> *get_long_frags();
        void read_long_frags(const char *fn);
        void merge_long_frags(VertexPath *aFrag);
        bool equal_frags(vector<VertexPath *> *f1, vector<VertexPath *> *f2);
        VertexPath *get_complement(VertexPath *path);

//        hic
        void read_hic_matrix(const char *fn);

    /* functionality */
        // double normalizeCoverage(double value, double mean, double std);
        vector<double> scaleILPCoef(vector<double> aCovs);
        int balancerILP(const char * lpFn);

        void addAllJuncsFromDB(JunctionDB *aJuncDB);

        double inferCoverage(Vertex * aSource, Vertex * aTarget);
        double weightedCredibility(Vertex * aVertex, bool aIsSource);
        double inferCredibility(Vertex * aSource, Vertex * aTarget);

        int strandCross(Vertex * aVertex);

        bool doesPathExists(Vertex * aStartVertex, Vertex * aEndVertex);

        void connectSourceSink();
        void countReachability(int * nReachability, VertexPath & notReachableVertices, Vertex * targetVertex);
        void findWithReachability(int * nReachability, VertexPath & notReachableVertices, int reachability);
        void reconnectNegativeToPositive(JunctionDB * aJuncDB, bool verbose = false);
        Vertex * getMostReachable(VertexPath & notReachableVertices, Vertex * targetVertex);
        Vertex * getFirstStrandCrossVertex(VertexPath & notReachableVertices);
        bool adjustReachability(VertexPath & notReachableVertices, Vertex * targetVertex, JunctionDB * aJuncDB, bool verbose = false);
        bool isReachable();
        void checkReachability(JunctionDB * aJuncDB, bool verbose = false);
        vector<Segment*> * extractReachableGraph(bool verbose = false);
        void addNormalJunctions();
        void checkInferredJunctionCredibility();
        void clearSegmentJunctionCredibility(Segment * aSegment);

        int findCircuit(Vertex * aVertex, VertexPath & pathVertices, EdgePath & pathEdges);
        void traverse(Vertex * startVertex, JunctionDB * aJuncDB);
        Edge * traverseNextEdge(Vertex * aStartVertex, VertexPath* vp, JunctionDB * aJuncDB);
        Edge * traverseWithHic(VertexPath* vp);
        double calculateHicInteraction(VertexPath* vp, Vertex* currentVertex);
//        decrease hic matrix when traversing
        void decreaseHicMatrix(VertexPath* vp, Edge* e);
        void decreaseHicInteraction(Vertex* v1, Vertex* v2);
        void traverseGraph(JunctionDB * aJuncDB);            // traverse all junctions and segments first
//        Traverse long path first
        void traverseWithLong(Vertex * aStartVertex, JunctionDB * aJuncDB);
        Vertex * traverseLongPath(Vertex * aStartVertex, VertexPath* vPath, bool skip_first);
//        find if long path in current graph (return how far the path can be traversed in current graph)
        int longPathLenInGraph(VertexPath* longPath);
        void isCircuitSimple(VertexPath * circuit, pair<int, int> & notSimpleIdx);
        void allCircuitsSimple(vector< tuple<int, int, int> > & notSimpleIdx);
        void extractCircuits();          // after traversing the whole graph, extract circuits from found paths
        void sortCircuits();
        
        void writeCircuits(const char * outFn);

        void generateHaploids();
        void writeHaploids(const char * outFn);

        Vertex * selectPrevVertex(Vertex * currentVertex, Vertex * targetVertex, JunctionDB * aJuncDB);
        Vertex * selectNextVertex(Vertex * currentVertex, Vertex * targetVertex, JunctionDB * aJuncDB);
        Edge * selectPrevEdge(Vertex * aTargetVertex, bool isTraversing = false);
        Edge * selectNextEdge(Vertex * aSourceVertex, bool isTraversing = false);

        /* print functions */
        void print(EdgePath& path);
        void print(VertexPath& path);
        void printCircuits();
        void printHaploids();
};

#endif
