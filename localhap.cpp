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
            ("op", "operate: check or solve", cxxopts::value<std::string>())
            ("juncdb", "Junction database", cxxopts::value<std::string>()->default_value(""))
            ("in_lh", "Input lh file", cxxopts::value<std::string>())
            ("out_lh", "Checked local hap input file, lh format", cxxopts::value<std::string>())
            ("lp_prefix", "ILP out file prefix, only for check", cxxopts::value<std::string>())
            ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("hic_matrix", "Segment Hic matrix file, only for solve", cxxopts::value<std::string>())
            ("tgs_order", "Segment tgs local order file, only for solve", cxxopts::value<std::string>())
            ("hap", "Haplotype out file, only for solve", cxxopts::value<std::string>())
            ("traversed", "traversed path out file, only for solve", cxxopts::value<std::string>())
            ("circuits", "Circuits out file, only for solve", cxxopts::value<std::string>())
            ("help", "Print usage")
            // BFB
            ("junc_info", "Use extra junction information", cxxopts::value<bool>()->default_value("false"))
            ("max_error", "The maximal acceptable rate", cxxopts::value<double>()->default_value("-1"))
            ("seq_mode", "Resolve a sequential bfb path without nested loops", cxxopts::value<bool>()->default_value("false"))
            ("edges", "Edges that indicate the evolution in single-cell data", cxxopts::value<std::string>()->default_value(""));
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    result["op"].as<std::string>().c_str();
    std::cout << result["op"].as<std::string>().c_str() << std::endl;

    if (strcmp(result["op"].as<std::string>().c_str(), "check") == 0) {
        const char *juncdbFn_sample = result["juncdb"].as<std::string>().c_str();
        const char *lhRawFn = result["in_lh"].as<std::string>().c_str();
        const char *lhCheckedFn = result["out_lh"].as<std::string>().c_str();
        const char *lpFn = result["lp_prefix"].as<std::string>().c_str();
        bool verbose = result["verbose"].as<bool>();
        try {
            auto *db_sample = new JunctionDB(juncdbFn_sample);
            db_sample->sortRecordEntry();
//            db_sample->print();

            Graph *g = new Graph(lhRawFn);
            g->calculateHapDepth();
            g->checkLowerBound();
            g->print();
            LocalGenomicMap *lgm = new LocalGenomicMap(g);
            cout << "checking reachability" << endl;
            g->print();
            lgm->checkReachability(db_sample, verbose);

            g->calculateCopyNum();
            g->checkOrphan();
            g->print();

            // lgm->checkInferredJunctionCredibility();
            cout << "connect ss" << endl;
            try {
                lgm->connectSourceSink();
            } catch (DuplicateJunctionException &e) {
                cout << e.what() << endl;
            }
//            g->calculateCopyNum();
            g->print();
            cout << "write" << endl;
            g->writeGraph(lhCheckedFn);
            cout << "write done" << endl;

            int feasible = lgm->balancerILP(lpFn);

            cout << "Done" << endl;
        } catch (DuplicateJunctionException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (SegmentDoesNotExistException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (JunctionDoesNotExistException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (BackwardReachSourceNegativeException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (BackwardReachSinkPositiveException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (ForwardReachSinkNegativeException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (ForwardReachSourcePositiveException &e) {
            cout << e.what() << endl;
            return 1;
            // } catch (...) {
            //     cout << "Logic error." << endl;
            //     return 1 ;
        }
        return 0;
    } else if (strcmp(result["op"].as<std::string>().c_str(), "solve") == 0) {
        const char *juncdbFn = result["juncdb"].as<std::string>().c_str();
        const char *lhFn = result["in_lh"].as<std::string>().c_str();
        // const char * lpOutFn = argv[4];
        const char *circuitsFn = result["circuits"].as<std::string>().c_str();
        const char *traversedFn = result["traversed"].as<std::string>().c_str();
        const char *hapFn = result["hap"].as<std::string>().c_str();
        bool verbose = result["verbose"].as<bool>();
        const char *longFragFn;
        const char *hicMatrix;
        try {
            JunctionDB *db = new JunctionDB(juncdbFn);
            db->sortRecordEntry();
            // db->writeDB("/home/tanbowen/LocalHap/test/test.db");
            // db->print();

            Graph *g = new Graph(lhFn);
            // g->checkOrphan();
            // g->print();
            // g->calculateCopyNum();
            // g->checkLowerBound();
            // g->print();
            // g->calculateCopyNum();
            g->print();
            LocalGenomicMap *lgm = new LocalGenomicMap(g);
            if (result.count("tgs_order")) {
                longFragFn = result["tgs_order"].as<std::string>().c_str();
                lgm->read_long_frags(longFragFn);
                lgm->setUsingLong(true);
                for (auto frag : *lgm->get_long_frags()) {
                    for(VertexPath* suvp : *frag.second) {
                        for (Vertex *v: *suvp) {
                            cout << v->getInfo() << " ";
                        }
                    }
                    cout << endl;
                }
            }
            if (result.count("hic_matrix")) {
                hicMatrix = result["hic_matrix"].as<std::string>().c_str();

                lgm->read_hic_matrix(hicMatrix);
            }
//             lgm->read_hic_matrix(hicMatrix);
//             return 0;
            // TODO some problems on checking reachability
            // especially with inversion
            // cout << "add normal" << endl;
            // lgm->addNormalJunctions();
            // g->print();
            // lgm->addNormalJunctions();
            // lgm->checkReachability(db, verbose);
            // lgm->addAllJuncsFromDB(db);
            // g->checkOrphan();
            // g->print();

            // lgm->checkInferredJunctionCredibility();
            // g->print();

            /* integer linear programming */
            try {
//                lgm->connectSourceSink();
            } catch (DuplicateJunctionException &e) {
                cout << e.what() << endl;
            }
            g->print();
            // int feasible = lgm->balancerILP(lpOutFn);
            // if (feasible >= 0) cout << "feasible" << endl;
            // if (feasible < 0) return 1;
            // cout << "After balancing: >>>>>>>>>>>>" << endl;
            // g->print();

            // delete g->getJunctions()->back();
            // g->getJunctions()->pop_back();
            // g->print();

            cout << "Identifying circuits..." << endl;
            lgm->traverseGraph(db);
            // g->print();
            // lgm->printCircuits();
            // cout << "Extracting circuits..." << endl;
            lgm->writeTraversedPath(traversedFn);
            lgm->extractCircuits();
            lgm->sortCircuits();
            lgm->divideCircuits();
            // g->print();
            if (verbose) {
               lgm->printCircuits();
            }

           lgm->writeCircuits(circuitsFn);

            cout << "Generating haploids..." << endl;
            g->print();
//            lgm->generateHaploids();
            // // cout << "Estimated haploids: " << endl;
            // // lgm->printHaploids();
//            lgm->writeHaploids(hapFn);
            // cout << "Done" << endl;
        } catch (DuplicateJunctionException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (SegmentDoesNotExistException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (JunctionDoesNotExistException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (ILPBalancerInfeasibleException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (BackwardReachSourceNegativeException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (BackwardReachSinkPositiveException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (ForwardReachSinkNegativeException &e) {
            cout << e.what() << endl;
            return 1;
        } catch (ForwardReachSourcePositiveException &e) {
            cout << e.what() << endl;
            return 1;
            // } catch (...) {
            //     cout << "Logic error." << endl;
            //     return 1 ;
        }
    }else if (strcmp(result["op"].as<std::string>().c_str(), "bfb") == 0) {
        const char *lhRawFn = result["in_lh"].as<std::string>().c_str();
        const char *lpFn = result["lp_prefix"].as<std::string>().c_str();
        const char *juncsFn = result["juncdb"].as<std::string>().c_str();
        const bool juncsInfo = result["junc_info"].as<bool>();//whether add extra junction iformation into ILP constrains
        const double maxError = result["max_error"].as<double>();
        const bool seqMode = result["seq_mode"].as<bool>();//whether use sequential mode (no small loop inserted in big loop)
        string edges = result["edges"].as<std::string>();

        // string bfbRes = "\n"+string(lhRawFn)+" "+string(juncsFn)+"\n";
        // ofstream outString;
        // outString.open("bfbPaths.txt",std::ios_base::app);
        // outString<<bfbRes;
        // outString.close();

        // get multiple .lh file names for single cell data
        vector<Graph*> graphs;
        vector<LocalGenomicMap*> lgms;
        unordered_map<string, int> graphIdx;
        // cout<<string(lhRawFn)<<endl;
        char* lhNames = (char*)lhRawFn;
        char* lhFn;
        while (lhFn = strtok_r(lhNames, ",", &lhNames)) {
            cout<<lhFn<<endl;            
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
        else {
            for(int i=0; i<evolution.size(); i++) {
                for(int j=i+1; j<evolution.size(); j++) evolution[i].push_back(j);
            }
        }
        for(int i=0; i<evolution.size(); i++) {
            cout<<i<<": ";
            for(int j: evolution[i]) cout<<j<<" ";
            cout<<endl;
        }
        // graph data structure for single .lh file
        int numGraphs = graphs.size();
        Graph *g = graphs[0];
        g->calculateHapDepth();
        g->calculateCopyNum();
        LocalGenomicMap *lgm = new LocalGenomicMap(g);
        //read options from input
        string mainChr;
        int insMode = 0, conMode = 0;//0: pre-BFB insertion/concatenation, 1: post-BFB insertion/concatenation
        vector<string> insChr, conChr;
        vector<int> startSegs;//starting segments for insertions
        int virusInfo[3] = {0};
        lgm->readBFBProps(mainChr, insMode, insChr, conMode, conChr, virusInfo, startSegs, lhRawFn);//read properties
        unordered_map<int, int> originalSegs;
        //Insertion Mode 1: pre-BFB insertion        
        if(insMode == 1) {
            lgm->insertBeforeBFB(g, insChr, originalSegs);
            delete lgm;
            lgm = new LocalGenomicMap(g);
        }
        //Concatenation Mode 1: pre-BFB concatenation
        else if(conMode == 1) {
            lgm->concatBeforeBFB(g, conChr, originalSegs);
            delete lgm;
            lgm = new LocalGenomicMap(g);
        }

        vector<Segment *> sources = *g->getMSources();
        vector<Segment *> sinks = *g->getMSinks();

        //construct segment intervals
        unordered_map<int,int> intervals;
        for(int i=0; i<sources.size(); i++) {
            for(int j=sources[i]->getId(); j<=sinks[i]->getId(); j++) {
                cout<<j<<"-"<<i<<endl;
                intervals.insert(pair<int,int>(j,i));
            }
        }

        //get information of third-generation data
        vector<vector<int>> components;
        lgm->readComponents(components, juncsFn, intervals);//third-generation data information
        for(int i=0; i<components.size(); i++) {
            for(int j=0; j<components[i].size(); j++) {
                cout<<components[i][j]<<" ";
            }
            cout<<endl;
        }
        
        //record target CN of segments
        vector<int> targetCN(g->getSegments()->size(),0);

        vector<vector<int>> bfbPaths;
        //construct bfb path on each chromosome
        for (int n=0; n<g->getMSources()->size(); n++) {
            //enumerate all the patterns and loops
            int startID = sources[n]->getId();
            int endID = sinks[n]->getId();
            vector<vector<int>> patterns, loops;
            vector<int> temp;
            lgm->combinations(startID,endID,2,patterns,temp);
            temp.clear();
            lgm->combinations(startID,endID,2,loops,temp);

            //construct mapping from pattern/loop to index
            map<string, int> variableIdx;
            for (int i=0;i<patterns.size();i++) {
                string key = "p:"+to_string(patterns[i][0])+","+to_string(patterns[i][1]);
                variableIdx[key] = i;
                cout<<variableIdx[key]<<" "<<key<<endl;
            }
            for (int i=0;i<loops.size();i++) {
                string key = "l:"+to_string(loops[i][0])+","+to_string(loops[i][1]);
                variableIdx[key] = i+patterns.size();
                cout<<variableIdx[key]<<" "<<key<<endl;
            }
            int numComp = variableIdx.size();

            //find copy number for both normal junctions and fold-back inversions
            vector<Junction *> inversions;
            double** juncCN = new double*[endID+1];
            lgm->getJuncCN(inversions, juncCN, *g, startID, endID);
            //check if there is any fold-back inversion
            cout<<"Junction CN"<<endl;
            double inversionCNSum = 0;
            for (int i=0;i<=endID;i++) {
                inversionCNSum += juncCN[i][1];
                cout<<i<<","<<i+1<<" "<<juncCN[i][0]<<"\t"<<i<<","<<i<<" "<<juncCN[i][1]<<endl;
            }
            //copy number of patterns and loops
            int* elementCN = new int[numComp*numGraphs];
            memset(elementCN, 0, numComp*numGraphs*sizeof(int));
            //pick components in the segment interval
            vector<vector<int>> validComponents;
            for(int i=0; i<components.size(); i++) {
                if(intervals[components[i][0]]==n) validComponents.push_back(components[i]);
            }

            if (abs(inversionCNSum)<0.000001&&validComponents.size()==0) {//no fold-back inversion
                for(int i=startID-1; i<endID; i++) targetCN[i] += 2;//compute target CN of segments
                vector<int> temp({startID, endID, endID, startID});
                lgm->editInversions(temp, inversions, juncCN, elementCN, variableIdx);
                bfbPaths.push_back(temp);
                continue;
            }
                        
            //construct ILP and generate .lp file for cbc
            if (numGraphs == 1) lgm->BFB_ILP(lpFn, patterns, loops, variableIdx, juncCN, validComponents, juncsInfo, maxError, seqMode);
            else lgm->BFB_ILP_SC(lpFn, patterns, loops, variableIdx, graphs, maxError, evolution);
            //reset variableIdx
            for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) 
                variableIdx[iter->first] = variableIdx[iter->first]%numComp;

            //run cbc under the directory containing test.lp
            string str = "cbc "+string(lpFn) +".lp solve solu "+string(lpFn)+".sol";
            const char *cmd = str.c_str();
            system(cmd);
            //read patterns and loops from .sol
            str = "./" + string(lpFn)+".sol";
            const char *solDir = str.c_str();
            ifstream solFile(solDir);
            if (!solFile) {
                cerr << "Cannot open file " << solDir << endl;
                exit(1);
            }
            
            string element, cn;
            bool infeasible = false;
            double cn_error = 0;
            ofstream errorString;
            errorString.open("./result.txt",std::ios_base::app);
            while (solFile >> element) {
                if(element == "Infeasible") {
                    infeasible = true;
                    break;
                }
                // if(element == "value") {
                //     solFile >> element;
                //     errorString<<element<<" ";
                // }
                if (element[0] == 'x') {                
                    int x = stoi(element.substr(1));
                    if (x < numComp*numGraphs) {//exclude epsilons
                        solFile >> cn;
                        int copynum = stoi(cn);
                        elementCN[x] = copynum;
                    }
                    else {//epsilons (errors)
                        if((x-variableIdx.size())%3==0) {
                            double seg_error;
                            solFile >> seg_error;
                            cn_error += seg_error;
                        }
                    }
                }
            }
            if(infeasible) {
                for(int i=startID-1; i<endID; i++) targetCN[i] += 2;//compute target CN of segments
                vector<int> temp({startID, endID, endID, startID});
                lgm->editInversions(temp, inversions, juncCN, elementCN, variableIdx);
                bfbPaths.push_back(temp);
                continue;
            }
            //compute target CN of segments based loop/pattern
            for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) {
                if(elementCN[iter->second] > 0) {
                    string key = iter->first;
                    int idx1 = stoi(key.substr(2, key.find(",")-2)), idx2 = stoi(key.substr(key.find(",")+1));
                    for(int i=idx1-1; i<idx2; i++) {
                        if(key[0]=='p') targetCN[i] += elementCN[iter->second];
                        else targetCN[i] += elementCN[iter->second]*2;
                    }
                }
            }

            // vector<int> segCN(endID, 0);
            // for (auto iter=variableIdx.begin();iter!=variableIdx.end();iter++) {
            //     if(elementCN[iter->second] > 0) {
            //         string key = iter->first;
            //         int idx1 = stoi(key.substr(2, key.find(",")-2)), idx2 = stoi(key.substr(key.find(",")+1));
            //         for(int i=idx1-1; i<idx2; i++) {
            //             if(key[0]=='p') segCN[i-1] += elementCN[iter->second];
            //             else segCN[i-1] += elementCN[iter->second]*2;
            //         }
            //     }
            // }
            // int cn_diff = 0;                        
            // for(int i=startID-1; i<endID; i++) {
            //     cn_diff += (*g->getSegments())[i]->getWeight()->getCopyNum()-segCN[i];
            // }
            
            // output errors
            // errorString<<cn_error<<endl;
            // errorString.close();
            // exit(0);
            
            //construct BFB DAG and find all topological orders
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
                //set up indegree
                int cnt = 0;
                for (int i = 0; i < num; i++) {
                    for (auto next = adj[i].begin(); next != adj[i].end(); next++) {
                        indeg[*next]++;                
                    }
                    cout<<i+1<<": ";
                    for (int j=0;j<adj[i].size();j++) {
                        cout<<adj[i][j]+1<<" ";
                    }
                    cout<<endl;
                }
                //find all topological orders in BFB DAG
                vector<int> res;
                vector<vector<int>> orders;
                lgms[k]->allTopologicalOrders(res, visited, num, indeg, adj, orders);
                cout<<"All topological orders: "<<endl;
                for (vector<int> bfb: orders) {
                    for (int i=0;i<bfb.size();i++)
                        cout<<bfb[i]+1<<" ";
                    cout<<endl;
                }
                //get one valid bfb path
                vector<int> path;
                lgms[k]->getBFB(orders, node2pat, node2loop, path);//get a valid BFB path          
                //output the text for visualization
                if(numGraphs==1) { 
                    lgms[k]->editInversions(path, inversions, juncCN, elementCN, variableIdx);//edit the imperfect fold-back inversions (with deletion)
                    if(insMode==1 || conMode==1) lgms[k]->printOriginalBFB(path, originalSegs);
                }
                else {
                    lgms[k]->getJuncCN(inversions, juncCN, *graphs[k], startID, endID);
                    lgms[k]->editInversions(path, inversions, juncCN, elementCN, variableIdx);
                    // lgms[k]->printBFB(path);
                }
                bfbPaths.push_back(path);
            }
        }              
        //output target CN
        ofstream targetCNString;
        targetCNString.open("./target_cn.txt",std::ios_base::app);
        double diffCN = 0;
        for(int i=0; i<targetCN.size(); i++) {
            diffCN += targetCN[i]-(*g->getSegments())[i]->getWeight()->getCopyNum();
            // targetCNString<<string(lhRawFn)<<"\t"<<i+1<<"\t"<< (*g->getSegments())[i]->getChrom()<<":"<<(*g->getSegments())[i]->getStart()<<"-"<<
            // (*g->getSegments())[i]->getEnd()<<"\t"<<(*g->getSegments())[i]->getWeight()->getCopyNum() <<"\t"<<targetCN[i]<<"\t"<<string(juncsFn) <<endl;
        }
        targetCNString<<string(lhRawFn)<<"\t"<<diffCN<<"\tnumSeg: "<<targetCN.size()<<endl;
        targetCNString.close();
        
        if(numGraphs > 1) exit(0);

        //Insertion Mode 2: post-BFB insertion     
        if(insMode == 2) lgm->insertAfterBFB(insChr, mainChr, startSegs, bfbPaths);
        //Concatenation Mode 2: post-BFB Concatenation 
        if(conMode == 2) lgm->concatAfterBFB(conChr, bfbPaths);
        
        //deal with insertion of virus
        if(virusInfo[0] != 0) {
            for(vector<int> path: bfbPaths) {
                cout<<"bfb path with virus: "<<endl;
                for(int i = 1; i<path.size(); i+=2) {
                    cout<<i<<": "<<virusInfo[0]<<" "<<virusInfo[1]<<" "<<virusInfo[2]<<endl;
                    if(path[i-1]<=virusInfo[1]&&virusInfo[2]<=path[i]) {
                        path.insert(path.begin()+i, {virusInfo[1],-1,-1,virusInfo[0],virusInfo[0],-1,-1,virusInfo[2]});
                        i += 8;
                    }
                    else if(path[i]<=virusInfo[1]&&virusInfo[2]<=path[i-1]) {
                        path.insert(path.begin()+i, {virusInfo[2],-1,-1,virusInfo[0],virusInfo[0],-1,-1,virusInfo[1]});
                        i += 8;
                    }
                    if(i<path.size()-1&& (path[i]==virusInfo[2]&&path[i+1]==virusInfo[1] 
                        || path[i]==virusInfo[1]&&path[i+1]==virusInfo[2])) {
                        path.insert(path.begin()+i+1, {-1,-1,virusInfo[0],virusInfo[0],-1,-1});
                        i += 6;
                    }                
                }    
                lgm->printBFB(path);
            }
        }

        //print the result
        // for (int i=0;i<bfbPaths.size();i++)
        //     lgm->printBFB(bfbPaths[i]);

    } else if (strcmp(result["op"].as<std::string>().c_str(), "bpm") == 0) {
        const char *lhRawFn = result["in_lh"].as<std::string>().c_str();
        Graph *g = new Graph(lhRawFn);
        LocalGenomicMap *lgm = new LocalGenomicMap(g);
        //Extract SV and normal edges and divide segments into two parts
        vector<Junction *> sv, normal;
        vector<Segment *> sources, targets;
        int num =  g->getSegments()->size();
        for (Junction *junc: *(g->getJunctions())){
            if (junc->getSource()->getId()%num +1  == junc->getTarget()->getId() &&
                junc->getSourceDir()=='+' && junc->getTargetDir()=='+') {                
                normal.push_back(junc);
            }
            else {
                sv.push_back(junc);
            }                
        }
        //Find the maximum bipartite matching
        vector<Junction *> selectedJunc;
        lgm->findMaxBPMatching(sv, selectedJunc);
        // for (Junction *junc: selectedJunc) {
        //     cout<<junc->getSource()->getId()<<" -> "<<junc->getTarget()->getId()<<endl;
        // }

        //Find the Eulerian Circuit with the selected and normal junctions
        vector<vector<int>> adj; //adjacency list
        for (int i=0; i<=num; i++) {
            vector<int> temp;
            adj.push_back(temp);
        }
        cout<<"SV: "<<endl;
        for (Junction* junc: selectedJunc) {
            cout<<junc->getSource()->getId()<<"->"<<junc->getTarget()->getId()<<endl;
            adj[junc->getSource()->getId()].push_back(junc->getTarget()->getId());
        }
        cout<<"Normal: "<<endl;
        for (Junction* junc: normal) {
            cout<<junc->getSource()->getId()<<"->"<<junc->getTarget()->getId()<<endl;
            adj[junc->getSource()->getId()].push_back(junc->getTarget()->getId());
        }

        lgm->findCircuits(adj);
        //lgm->constructCircuits(adj);

        //Traverse the graph with selected junc
        // JunctionDB *db = new JunctionDB(selectedJunc);
        // db->sortRecordEntry();
    
    } else if (strcmp(argv[1], "split") == 0) {
        const char *lhFn = argv[2];
        const char *chrom = argv[3];
        const char *lhOut = argv[4];
        bool verbose = (strcmp(argv[6], "--verbose") == 0) ? true : false;
        Graph *g = new Graph(lhFn);
        g->calculateHapDepth();
        g->print();

        LocalGenomicMap *lgm = new LocalGenomicMap(g);

        auto segs = lgm->extractReachableGraph(verbose);
        g->writePartialGraph(segs, lhOut);
        g->print();
        cout << "write" << endl;
        g->writeGraph(lhOut);
        cout << "write done" << endl;

    }
}
