#ifndef _LOCAL_GENOMIC_MAP_H_
#define _LOCAL_GENOMIC_MAP_H_

#include <unordered_map>
#include <map>
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
//TODO, There are some functions need to be reconstruct into util
class LocalGenomicMap {
protected:
    Graph *mGraph;

    vector<VertexPath *> *mRawPathVertices;
    vector<EdgePath> *mRawPathEdges;

    vector<VertexPath *> *mCircuits;
    unordered_map<int, vector<VertexPath *> *> *traversedCircuits;
    unordered_map<int, vector<VertexPath *> *> *dividedCircuits;
    unordered_map<int, vector<VertexPath *> *> *dividedHaploids;
    vector<VertexPath *> *mHaploids;

    unordered_map<int, vector<VertexPath *> *> *mLongFrags;
//    vector<VertexPath *>* mTraversedPath;

    bool usingLong;
    bool usingHic;

    double **hicMatrix;
    double **decreaseMatrix;

public:
    LocalGenomicMap(Graph *aGraph);

    ~LocalGenomicMap();

    /* getter and setter */
    Graph *getGraph();

    void setGraph(Graph *aGraph);

    bool isUsingLong();

    bool isUsingHic();

    void setUsingLong(bool v);

    void setUsingHic(bool);

    /* Long fragments for TGS */
    unordered_map<int, vector<VertexPath *> *> *get_long_frags();

    void read_long_frags(const char *fn);

    void merge_long_frags(vector<VertexPath *> * partitionLong, VertexPath *aFrag);

    bool equal_frags(vector<VertexPath *> *f1, vector<VertexPath *> *f2);

    VertexPath *get_complement(VertexPath *path);

//        hic
    void read_hic_matrix(const char *fn);

    /* functionality */
    // double normalizeCoverage(double value, double mean, double std);
    vector<double> scaleILPCoef(vector<double> aCovs);

    int balancerILP(const char *lpFn);

    void addAllJuncsFromDB(JunctionDB *aJuncDB);

    double inferCoverage(Vertex *aSource, Vertex *aTarget);

    double weightedCredibility(Vertex *aVertex, bool aIsSource);

    double inferCredibility(Vertex *aSource, Vertex *aTarget);

    int strandCross(Vertex *aVertex);

    bool doesPathExists(Vertex *aStartVertex, Vertex *aEndVertex);

    void connectSourceSink();

    void countReachability(int *nReachability, VertexPath &notReachableVertices, Vertex *targetVertex);

    void findWithReachability(int *nReachability, VertexPath &notReachableVertices, int reachability);

    void reconnectNegativeToPositive(JunctionDB *aJuncDB, bool verbose = false);

    Vertex *getMostReachable(VertexPath &notReachableVertices, Vertex *targetVertex);

    Vertex *getFirstStrandCrossVertex(VertexPath &notReachableVertices);

    bool adjustReachability(VertexPath &notReachableVertices, Vertex *targetVertex, JunctionDB *aJuncDB,
                            bool verbose = false);

    bool isReachable();

    void checkReachability(JunctionDB *aJuncDB, bool verbose = false);

    vector<Segment *> *extractReachableGraph(bool verbose = false);

    void addNormalJunctions();

    void checkInferredJunctionCredibility();

    void clearSegmentJunctionCredibility(Segment *aSegment);

    int findCircuit(Vertex *aVertex, VertexPath &pathVertices, EdgePath &pathEdges);

    void traverse(Vertex *startVertex, JunctionDB *aJuncDB);

    Edge *traverseNextEdge(Vertex *aStartVertex, VertexPath *vp, JunctionDB *aJuncDB);

    Edge *traverseNextEdgeByPartition(Vertex *aStartVertex, VertexPath *vp, JunctionDB *aJuncDB, int *partitionStart,
                                      int *partitionEnd);

    Edge *traverseWithHic(VertexPath *vp);

    pair<int, int> findPartition(int id);
    pair<int, int> findLongPathPartition(VertexPath* vp);

    double calculateHicInteraction(VertexPath *vp, Vertex *currentVertex);

//        decrease hic matrix when traversing
    void decreaseHicMatrix(VertexPath *vp, Edge *e);

    void decreaseHicInteraction(Vertex *v1, Vertex *v2);

    void traverseGraphByPartition(JunctionDB *aJuncDB);

    void traverseGraph(JunctionDB *aJuncDB);            // traverse all junctions and segments first
//        Traverse long path first
    void traverseWithLong(Vertex *aStartVertex, JunctionDB *aJuncDB);

    Vertex *traverseLongPath(Vertex *aStartVertex, VertexPath *vPath, int *partitionStart, int *partitionEnd);

//        find if long path in current graph (return how far the path can be traversed in current graph)
    int longPathLenInGraph(VertexPath *longPath);

    bool checkPartition(int id, int* partitionStart, int* partitionEnd);
    bool checkCommon(int id, int* partitionStart, int* partitionEnd);
    bool vReachable(bool isBackwardSourceReachable, bool isForwardSinkReachable, bool isBackwardSinkReachable, bool isForwardSourceReachable);

    void isCircuitSimple(VertexPath *circuit, pair<int, int> &notSimpleIdx);

    void allCircuitsSimple(vector<tuple<int, int, int> > &notSimpleIdx);

    void extractCircuits();          // after traversing the whole graph, extract circuits from found paths
    void sortCircuits();
    void divideCircuits();

    void writeCircuits(const char *outFn);
    void writeTraversedPath(const char *outFn);

    void generateHaploids();

    void writeHaploids(const char *outFn);

    Vertex *selectPrevVertex(Vertex *currentVertex, Vertex *targetVertex, JunctionDB *aJuncDB);

    Vertex *selectNextVertex(Vertex *currentVertex, Vertex *targetVertex, JunctionDB *aJuncDB);

    Edge *selectPrevEdge(Vertex *aTargetVertex, bool isTraversing = false);

    Edge *selectNextEdge(Vertex *aSourceVertex, bool isTraversing = false);
    Edge *selectNextEdgeByPartition(int partitionStart, int partitionEnd,Vertex *aSourceVertex, bool isTraversing = false);


    /* print functions */
    void print(EdgePath &path);

    void print(VertexPath &path);

    void printCircuits();

    void printHaploids();

    /* BFB functions*/
    void combinations(int start, int end, int len, vector<vector<int>> &per, vector<int> temp);
    void constructDAG(vector<vector<int>> &adj, vector<vector<int>> &node2pat, vector<vector<int>> &node2loop, map<string, int> &variableIdx, int *elementCN);
    void allTopologicalOrders(vector<int> &res, bool visited[], int num, int indeg[], vector<vector<int>> &adj, vector<vector<int>> &orders);
    void printBFB(VertexPath *path);
    void imperfectFBI(VertexPath *bkpPath, unordered_map<int, Junction*> &inversions);
    void getBFB(vector<vector<int>> &orders, vector<vector<int>> &node2pat, vector<vector<int>> &node2loop, VertexPath *path, unordered_map<int, Junction*> &inversions, bool isReversed, bool printAll);
    void getIndelBias(int startSegID, int endSegID);
    void indelBFB(VertexPath *path, int startSegID, int endSegID);
    void virusBFB(VertexPath *path, unordered_map<Segment*, Segment*> &originalSegs, vector<Junction *> unusedSV);
    void readBFBProps(string &mainChr, int &insMode, vector<string> &insChr, int &conMode, vector<string> &conChr, vector<int> &startSegs, const char *lhRawFn);
    void getJuncCN(unordered_map<int, Junction*> &inversions, double** juncCN, Graph &graph, int startSegID, int endSegID);
    void translocationBFB(vector<VertexPath*> paths, VertexPath *res, string mainChr);
    void insertBeforeBFB(Graph*& g, vector<string>& insChr, unordered_map<Segment*, Segment*>& originalSegs, vector<Junction *>& unusedSV);
    void concatBeforeBFB(Graph*& g, vector<string>& conChr, unordered_map<Segment*, Segment*>& originalSegs, vector<Junction *>& unusedSV);
    void BFB_ILP(const char *lpFn, vector<vector<int>> &patterns, vector<vector<int>> &loops, map<string, int> &variableIdx, double** juncCN, vector<vector<int>> &components, bool juncsInfo, int bias);
    void BFB_ILP_SC(const char *lpFn, vector<vector<int>> &patterns, vector<vector<int>> &loops, map<string, int> &variableIdx, vector<Graph*> graphs, vector<vector<int>> &evolution);
    void readComponents(vector<vector<int>>& res, unordered_map<Segment*, Segment*>& originalSegs, const char *juncsFn);

};

#endif
