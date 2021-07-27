#include <iostream>
#include <cstring>

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
            ("juncdb", "Junction database", cxxopts::value<std::string>())
            ("in_lh", "Input lh file", cxxopts::value<std::string>())
            ("out_lh", "Checked local hap input file, lh format", cxxopts::value<std::string>())
            ("lp_prefix", "ILP out file prefix, only for check", cxxopts::value<std::string>())
            ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("hic_matrix", "Segment Hic matrix file, only for solve", cxxopts::value<std::string>())
            ("tgs_order", "Segment tgs local order file, only for solve", cxxopts::value<std::string>())
            ("hap", "Haplotype out file, only for solve", cxxopts::value<std::string>())
            ("circuits", "Circuits out file, only for solve", cxxopts::value<std::string>())
            ("help", "Print usage");
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
                for (VertexPath *frag : *(lgm->get_long_frags())) {
                    for (Vertex *v: *frag) {
                        cout << v->getInfo() << " ";
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
            lgm->generateHaploids();
            // // cout << "Estimated haploids: " << endl;
            // // lgm->printHaploids();
            lgm->writeHaploids(hapFn);
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
