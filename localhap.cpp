#include <iostream>
#include <fstream>
#include <cstring>
#include <set>
#include <map>
#include <string>
#include <stack>

#include "Graph.hpp"
#include "Exceptions.hpp"
#include "LocalGenomicMap.hpp"
#include "JunctionDB.hpp"
#include "cxxopts.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    cxxopts::Options options("localhap", "Local Haplotype constructer");

    options.add_options()
            ("op", "Operate: bfb", cxxopts::value<std::string>())
            // ("out_lh", "Checked local hap input file, lh format", cxxopts::value<std::string>())
            // ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            // ("hic_matrix", "Segment Hic matrix file, only for solve", cxxopts::value<std::string>())
            // ("tgs_order", "Segment tgs local order file, only for solve", cxxopts::value<std::string>())
            // ("hap", "Haplotype out file, only for solve", cxxopts::value<std::string>())
            // ("traversed", "traversed path out file, only for solve", cxxopts::value<std::string>())
            // ("circuits", "Circuits out file, only for solve", cxxopts::value<std::string>())
            
            // BFB
            ("in_lh", "Input .lh file (required)", cxxopts::value<std::string>())
            ("lp_prefix", "ILP output file prefix", cxxopts::value<std::string>())
            ("juncdb", "Input .junc file with linkage information from linked/long reads", cxxopts::value<std::string>()->default_value(""))
            ("junc_info", "Whether use linked/long reads information in ILP (Default: false)", cxxopts::value<bool>()->default_value("false"))
            ("reversed", "Find BFB paths starting from the negative strand (Default: false)", cxxopts::value<bool>()->default_value("false"))
            ("max_error", "The maximal acceptable error rate (real number)", cxxopts::value<double>()->default_value("-1"))
            ("seq_mode", "Resolve a sequential BFB path without nested loops (Default: false)", cxxopts::value<bool>()->default_value("false"))
            // ("edges", "Edges that indicate the evolution of single-cell data", cxxopts::value<std::string>()->default_value(""))
            ("help", "Print usage");
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    result["op"].as<std::string>().c_str();
    std::cout << result["op"].as<std::string>().c_str() << std::endl;

    if (strcmp(result["op"].as<std::string>().c_str(), "bfb") == 0) {
        const char *lhRawFn = result["in_lh"].as<std::string>().c_str();// path to .lh file
        const char *lpFn = result["lp_prefix"].as<std::string>().c_str();// sample name
        const char *juncsFn = result["juncdb"].as<std::string>().c_str();// extra junction information from TGS data
        const bool juncsInfo = result["junc_info"].as<bool>();// indicate whether add extra junction iformation into ILP constrains
        const bool isReversed = result["reversed"].as<bool>();// the reference strand: true - forward; flase - backward
        const double maxError = result["max_error"].as<double>();// upper boundary of ILP error
        const bool seqMode = result["seq_mode"].as<bool>();// indicate whether use sequential mode (no small loop inserted in big loop)
        // string edges = result["edges"].as<std::string>();// relationship among sub-clones in single-cell data e.g. 1:2,1:3
        string edges = "";
       
        /* get multiple .lh file names for single-cell data */
        Graph *g = new Graph(lhRawFn);
        g->calculateHapDepth();
        g->calculateCopyNum();
        LocalGenomicMap *lgm = new LocalGenomicMap(g);
        
        // read options from input
        string mainChr;
        int insMode = 0, conMode = 0;// 0: pre-BFB insertion/concatenation, 1: post-BFB insertion/concatenation
        vector<string> insChr, conChr;
        vector<int> startSegs;// starting segments for insertions
        lgm->readBFBProps(mainChr, insMode, insChr, conMode, conChr, startSegs, lhRawFn);//read properties
        unordered_map<Segment*, Segment*> originalSegs;// mapping rearranged segments into original ones
        // Insertion Mode 1: pre-BFB insertion
        vector<Junction *> unusedSV;
        if(insMode == 1) {
            lgm->insertBeforeBFB(g, insChr, originalSegs, unusedSV);
            delete lgm;
            lgm = new LocalGenomicMap(g);
        }
        // Concatenation Mode 1: pre-BFB concatenation
        else if(conMode == 1) {
            lgm->concatBeforeBFB(g, conChr, originalSegs, unusedSV);
            delete lgm;
            lgm = new LocalGenomicMap(g);
        }

        vector<Segment *> sources = *g->getMSources();
        vector<Segment *> sinks = *g->getMSinks();
        vector<Segment *> segs = *g->getSegments();
        // set segment partitions
        for(int i=0; i<sources.size(); i++) {
            for(int j=sources[i]->getId(); j<=sinks[i]->getId(); j++) {
                segs[j-1]->setPartition(i);
            }
        }

        // get information of third-generation data
        vector<vector<int>> components;
        lgm->readComponents(components, juncsFn);// third-generation data information

        vector<vector<int>> bfbPaths;
        vector<VertexPath*> paths;
        //record target CN of segments
        vector<int> targetCN(g->getSegments()->size(),0);
        // construct bfb path on each chromosome
        for(int n = 0; n < sinks.size(); n++) {
            // enumerate all the patterns and loops
            int startID = sources[n]->getId();
            int endID = sinks[n]->getId();
            vector<vector<int>> patterns, loops;
            vector<int> temp;
            lgm->combinations(startID,endID,2,patterns,temp);
            temp.clear();
            lgm->combinations(startID,endID,2,loops,temp);

            // construct mapping from pattern/loop to index
            map<string, int> variableIdx;
            for (int i=0;i<patterns.size();i++) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                variableIdx[key] = i;
                // cout<<variableIdx[key]<<" "<<key<<endl;
            }
            for (int i=0;i<loops.size();i++) {
                string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
                variableIdx[key] = i+patterns.size();
                // cout<<variableIdx[key]<<" "<<key<<endl;
            }
            int numComp = variableIdx.size();

            // find copy number for both normal junctions and fold-back inversions
            unordered_map<int, Junction*> inversions;
            double** juncCN = new double*[endID+1];
            lgm->getJuncCN(inversions, juncCN, *g, startID, endID);

            // compute bias from imperct FBI and intra-chromosomal SV
            int bias = 1;
            for(int i = startID; i <= endID; i++) {
                if(juncCN[i][1] > 0) {
                    if(inversions[i]->getSource() != inversions[i]->getTarget()) bias += int(juncCN[i][1])%2;
                }
            }
            lgm->getIndelBias(startID, endID);
            
            // check if there is any fold-back inversion
            double inversionCNSum = 0;
            for (int i=0;i<=endID;i++) {
                inversionCNSum += juncCN[i][1];
            }
            //copy number of patterns and loops
            int* elementCN = new int[numComp];
            memset(elementCN, 0, numComp*sizeof(int));
            //pick components in the segment interval
            vector<vector<int>> validComponents;
            for(int i=0; i<components.size(); i++) {
                if(g->getSegmentById(components[i][0])->getPartition() == n) 
                    validComponents.push_back(components[i]);
            }

            if (abs(inversionCNSum)<0.000001&&validComponents.size()==0) {// no fold-back inversion
                VertexPath *temp = new VertexPath();
                for(int i = startID; i <= endID; i++) temp->push_back(g->getSegmentById(i)->getPositiveVertex());
                lgm->printBFB(temp);
                paths.push_back(temp);
                continue;
            }
                        
            // construct ILP and generate .lp file for cbc
            lgm->BFB_ILP(lpFn, patterns, loops, variableIdx, juncCN, validComponents, juncsInfo, maxError, seqMode, bias);
            // reset variableIdx (for single-cell data)
            for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) 
                variableIdx[iter->first] = variableIdx[iter->first]%numComp;

            // run cbc under the directory containing test.lp
            string str = "cbc "+string(lpFn) +".lp solve solu "+string(lpFn)+".sol";
            const char *cmd = str.c_str();
            system(cmd);
            // return 0;
            // read patterns and loops from .sol file
            str = "./" + string(lpFn)+".sol";
            const char *solDir = str.c_str();
            ifstream solFile(solDir);
            if (!solFile) {
                cerr << "ILP error: cannot open file " << solDir << endl;
                exit(1);
            }
            
            string element, cn;
            bool infeasible = false;
            while (solFile >> element) {
                if(element == "Infeasible") {
                    infeasible = true;
                    break;
                }
                if (element[0] == 'x') {                
                    int x = stoi(element.substr(1));
                    if (x < numComp) {// exclude epsilons
                        solFile >> cn;
                        int copynum = stoi(cn);
                        elementCN[x] = copynum;
                    }
                }
            }
            if(infeasible) {
                VertexPath *temp = new VertexPath();
                for(int i = startID; i <= endID; i++) temp->push_back(g->getSegmentById(i)->getPositiveVertex());
                lgm->printBFB(temp);
                cout<<"ILP is unsolvable.\n";
                paths.push_back(temp);
                continue;
            }
            //compute target CN of segments based on loop/pattern
            for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) {
                if(elementCN[iter->second] > 0) {
                    string key = iter->first;
                    cout<<"X"<<variableIdx[key]<<" "+key<<" CN: "<<elementCN[iter->second]<<endl;
                    int idx1 = stoi(key.substr(2, key.find(",")-2)), idx2 = stoi(key.substr(key.find(",")+1));
                    for(int i=idx1-1; i<idx2; i++) {
                        if(key[0]=='p') targetCN[i] += elementCN[iter->second];
                        else targetCN[i] += elementCN[iter->second]*2;
                    }
                }
            }
            
            // construct BFB DAG and find all topological orders
            vector<vector<int>> adj, node2pat, node2loop;
            lgm->constructDAG(adj, node2pat, node2loop, variableIdx, elementCN);
            int num = adj.size();
            bool *visited = new bool[num];
            int *indeg = new int[num];
            for (int i = 0; i < num; i++) {
                visited[i] = false;
                indeg[i] = 0;
            }
            // set up indegree
            int cnt = 0;
            for (int i = 0; i < num; i++) {
                for (auto next = adj[i].begin(); next != adj[i].end(); next++) {
                    indeg[*next]++;                
                }
            }
            // find all topological orders in BFB DAG
            vector<int> res;
            vector<vector<int>> orders;
            lgm->allTopologicalOrders(res, visited, num, indeg, adj, orders);
            for(vector<int> bfb: orders) {
                for(int i: bfb) cout<<i<<" ";
                cout<<endl;
            }
            // get one valid bfb path
            VertexPath *path = new VertexPath();
            lgm->getBFB(orders, node2pat, node2loop, path, inversions, isReversed);//get a valid BFB path
            lgm->indelBFB(path, startID, endID);
            if(insMode == 1 || conMode == 1) {
                // restore path segment into original segments
                for(auto iter = path->begin(); iter != path->end(); iter++) {
                    if((*iter)->getDir() == '+') 
                        *iter = originalSegs[(*iter)->getSegment()]->getPositiveVertex();
                    else
                        *iter = originalSegs[(*iter)->getSegment()]->getNegativeVertex();
                }
                lgm->printBFB(path);
            }
            paths.push_back(path);
        }
        // BFB-TRX mode
        VertexPath *res = new VertexPath();
        if(insMode == 2 || conMode == 2)  lgm->translocationBFB(paths, res, mainChr);

    } else if (strcmp(result["op"].as<std::string>().c_str(), "sc_bfb") == 0) {
       const char *lhRawFn = result["in_lh"].as<std::string>().c_str();// path to .lh file
        const char *lpFn = result["lp_prefix"].as<std::string>().c_str();// sample name
        const char *juncsFn = result["juncdb"].as<std::string>().c_str();// extra junction information from TGS data
        const bool juncsInfo = result["junc_info"].as<bool>();// indicate whether add extra junction iformation into ILP constrains
        const bool isReversed = result["reversed"].as<bool>();// the reference strand: true - forward; flase - backward
        const double maxError = result["max_error"].as<double>();// upper boundary of ILP error
        const bool seqMode = result["seq_mode"].as<bool>();// indicate whether use sequential mode (no small loop inserted in big loop)
        // string edges = result["edges"].as<std::string>();// relationship among sub-clones in single-cell data e.g. 1:2,1:3
        string edges = "";
       
        /* get multiple .lh file names for single-cell data */
        vector<Graph*> graphs;
        vector<LocalGenomicMap*> lgms;
        unordered_map<string, int> graphIdx;
        char* lhNames = (char*)lhRawFn;
        char* lhFn;
        while (lhFn = strtok_r(lhNames, ",", &lhNames)) {            
            Graph *g = new Graph(lhFn);
            graphIdx[lhFn] = graphs.size();
            graphs.push_back(g);
            g->calculateHapDepth();
            g->calculateCopyNum();
            lgms.push_back(new LocalGenomicMap(g));
        }

        // construct an adjacency list for single-cell evolution
        vector<vector<int>> evolution(graphs.size(), vector<int>());
        if(edges.length() > 0) {
            size_t pos = 0;
            string token;
            while ((pos = edges.find(",")) != std::string::npos) {
                token = edges.substr(0, pos);
                edges.erase(0, pos + 1);
                pos = edges.find(":");
                evolution[graphIdx[token.substr(0, pos)]].push_back(graphIdx[token.substr(pos+1)]);            
            }
            pos = edges.find(":");
            evolution[graphIdx[edges.substr(0, pos)]].push_back(graphIdx[edges.substr(pos+1)]);
        }
        else {// default: each pair of sub-clones share some similar patterns
            for(int i=0; i<evolution.size(); i++) {
                for(int j=i+1; j<evolution.size(); j++) evolution[i].push_back(j);
            }
        }

        /*graph data structure for single .lh file */
        int numGraphs = graphs.size();
        Graph *g = graphs[0];
        g->calculateHapDepth();
        g->calculateCopyNum();
        LocalGenomicMap *lgm = new LocalGenomicMap(g);

        // read options from input
        string mainChr;
        int insMode = 0, conMode = 0;// 0: pre-BFB insertion/concatenation, 1: post-BFB insertion/concatenation
        vector<string> insChr, conChr;
        vector<int> startSegs;// starting segments for insertions
        lgm->readBFBProps(mainChr, insMode, insChr, conMode, conChr, startSegs, lhRawFn);//read properties

        vector<Segment *> sources = *g->getMSources();
        vector<Segment *> sinks = *g->getMSinks();
        vector<Segment *> segs = *g->getSegments();
        // set segment partitions
        for(int i=0; i<sources.size(); i++) {
            for(int j=sources[i]->getId(); j<=sinks[i]->getId(); j++) {
                segs[j-1]->setPartition(i);
            }
        }

        vector<vector<int>> bfbPaths;
        vector<VertexPath*> paths;
        //record target CN of segments
        vector<int> targetCN(g->getSegments()->size(),0);
        // construct bfb path on each chromosome
        for (int n=0; n<g->getMSources()->size(); n++) {
            // enumerate all the patterns and loops
            int startID = sources[n]->getId();
            int endID = sinks[n]->getId();
            vector<vector<int>> patterns, loops;
            vector<int> temp;
            lgm->combinations(startID,endID,2,patterns,temp);
            temp.clear();
            lgm->combinations(startID,endID,2,loops,temp);

            // construct mapping from pattern/loop to index
            map<string, int> variableIdx;
            for (int i=0;i<patterns.size();i++) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                variableIdx[key] = i;
                // cout<<variableIdx[key]<<" "<<key<<endl;
            }
            for (int i=0;i<loops.size();i++) {
                string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
                variableIdx[key] = i+patterns.size();
                // cout<<variableIdx[key]<<" "<<key<<endl;
            }
            int numComp = variableIdx.size();

            // find copy number for both normal junctions and fold-back inversions
            unordered_map<int, Junction*> inversions;
            double** juncCN = new double*[endID+1];
            lgm->getJuncCN(inversions, juncCN, *g, startID, endID);

            // compute bias from imperct FBI and intra-chromosomal SV
            int bias = 1;
            for(int i = startID; i <= endID; i++) {
                if(juncCN[i][1] > 0) {
                    if(inversions[i]->getSource() != inversions[i]->getTarget()) bias += int(juncCN[i][1])%2;
                }
            }
            lgm->getIndelBias(startID, endID);
            
            // check if there is any fold-back inversion
            double inversionCNSum = 0;
            for (int i=0;i<=endID;i++) {
                inversionCNSum += juncCN[i][1];
            }
            //copy number of patterns and loops
            int* elementCN = new int[numComp*numGraphs];
            memset(elementCN, 0, numComp*numGraphs*sizeof(int));

            if (abs(inversionCNSum)<0.000001) {// no fold-back inversion
                VertexPath *temp = new VertexPath();
                for(int i = startID; i <= endID; i++) temp->push_back(g->getSegmentById(i)->getPositiveVertex());
                lgm->printBFB(temp);
                paths.push_back(temp);
                continue;
            }
                        
            // construct ILP and generate .lp file for cbc
            lgm->BFB_ILP_SC(lpFn, patterns, loops, variableIdx, graphs, maxError, evolution);
            // reset variableIdx (for single-cell data)
            for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) 
                variableIdx[iter->first] = variableIdx[iter->first]%numComp;

            // run cbc under the directory containing test.lp
            string str = "cbc "+string(lpFn) +".lp solve solu "+string(lpFn)+".sol";
            const char *cmd = str.c_str();
            system(cmd);
            // return 0;
            // read patterns and loops from .sol file
            str = "./" + string(lpFn)+".sol";
            const char *solDir = str.c_str();
            ifstream solFile(solDir);
            if (!solFile) {
                cerr << "ILP error: cannot open file " << solDir << endl;
                exit(1);
            }
            
            string element, cn;
            bool infeasible = false;
            while (solFile >> element) {
                if(element == "Infeasible") {
                    infeasible = true;
                    break;
                }
                if (element[0] == 'x') {                
                    int x = stoi(element.substr(1));
                    if (x < numComp*numGraphs) {// exclude epsilons
                        solFile >> cn;
                        int copynum = stoi(cn);
                        elementCN[x] = copynum;
                    }
                }
            }
            if(infeasible) {
                VertexPath *temp = new VertexPath();
                for(int i = startID; i <= endID; i++) temp->push_back(g->getSegmentById(i)->getPositiveVertex());
                lgm->printBFB(temp);
                cout<<"ILP is unsolvable.\n";
                paths.push_back(temp);
                continue;
            }
            //compute target CN of segments based on loop/pattern
            for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) {
                if(elementCN[iter->second] > 0) {
                    string key = iter->first;
                    cout<<"X"<<variableIdx[key]<<" "+key<<" CN: "<<elementCN[iter->second]<<endl;
                    int idx1 = stoi(key.substr(2, key.find(",")-2)), idx2 = stoi(key.substr(key.find(",")+1));
                    for(int i=idx1-1; i<idx2; i++) {
                        if(key[0]=='p') targetCN[i] += elementCN[iter->second];
                        else targetCN[i] += elementCN[iter->second]*2;
                    }
                }
            }
            
            // construct BFB DAG and find all topological orders
            for (int k = 0; k < numGraphs; k++) {
                vector<vector<int>> adj, node2pat, node2loop;
                for (int i = 0;i < numComp/2; i++) {
                    string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                    variableIdx[key] = i + k*numComp;
                    key = "l:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                    variableIdx[key] = i + numComp/2 + k*numComp;
                }
                lgms[k]->constructDAG(adj, node2pat, node2loop, variableIdx, elementCN);
                int num = adj.size();
                bool *visited = new bool[num];
                int *indeg = new int[num];
                for (int i = 0; i < num; i++) {
                    visited[i] = false;
                    indeg[i] = 0;
                }
                // set up indegree
                int cnt = 0;
                for (int i = 0; i < num; i++) {
                    for (auto next = adj[i].begin(); next != adj[i].end(); next++) {
                        indeg[*next]++;                
                    }
                }
                // find all topological orders in BFB DAG
                vector<int> res;
                vector<vector<int>> orders;
                lgms[k]->allTopologicalOrders(res, visited, num, indeg, adj, orders);
                for(vector<int> bfb: orders) {
                    for(int i: bfb) cout<<i<<" ";
                    cout<<endl;
                }
                // get one valid bfb path
                VertexPath *path = new VertexPath();
                lgms[k]->getBFB(orders, node2pat, node2loop, path, inversions, isReversed);//get a valid BFB path
                lgms[k]->indelBFB(path, startID, endID);
                paths.push_back(path);
            }
        }
    }
}
