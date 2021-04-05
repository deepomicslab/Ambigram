#include <iostream>
#include <cstring>

#include "Graph.hpp"
#include "Exceptions.hpp"
#include "LocalGenomicMap.hpp"
#include "JunctionDB.hpp"

using namespace std;

int main(int argc, char * argv[]) {
    if (strcmp(argv[1], "check") == 0) {
        const char * juncdbFn_sample = argv[2];
        // const char * juncdbFn_ref = argv[3];
        const char * lhRawFn = argv[3];
        const char * lhCheckedFn = argv[4];
        const char * lpFn = argv[5];
        bool verbose = (strcmp(argv[6], "--verbose") == 0) ? true : false;
        try {
            JunctionDB * db_sample = new JunctionDB(juncdbFn_sample);
            db_sample->sortRecordEntry();
            db_sample->print();

            Graph * g = new Graph(lhRawFn); 
            g->calculateHapDepth();
            g->checkLowerBound();
            g->print();

            LocalGenomicMap * lgm = new LocalGenomicMap(g);
            
            // lgm->addNormalJunctions();
            lgm->addAllJuncsFromDB(db_sample);
            // lgm->reconnectNegativeToPositive(db_sample, true);
            cout << "checking reachability" << endl;
            g->print();
            // exit(1);
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
            g->calculateCopyNum();
            g->print();
            cout << "write" << endl;
            g->writeGraph(lhCheckedFn);
            cout << "write done" << endl;

            int feasible = lgm->balancerILP(lpFn);

            cout << "Done" << endl;
        } catch (DuplicateJunctionException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (SegmentDoesNotExistException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (JunctionDoesNotExistException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (BackwardReachSourceNegativeException & e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (BackwardReachSinkPositiveException & e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (ForwardReachSinkNegativeException & e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (ForwardReachSourcePositiveException & e) {
            cout << e.what() << endl;
            return 1 ;
            // } catch (...) {
            //     cout << "Logic error." << endl;
            //     return 1 ;
        }
        return 0;
    } else if (strcmp(argv[1], "solve") == 0) {
        const char * juncdbFn = argv[2];
        const char * lhFn = argv[3];
        // const char * lpOutFn = argv[4];
        const char * circuitsFn = argv[4];
        const char * hapFn = argv[5];
        bool verbose = (strcmp(argv[6], "--verbose") == 0) ? true : false;
        // const char * longFragFn = argv[7];
        try {
            JunctionDB * db = new JunctionDB(juncdbFn);
            db->sortRecordEntry();
            // db->writeDB("/home/tanbowen/LocalHap/test/test.db");
            // db->print();

            Graph * g = new Graph(lhFn);
            // g->checkOrphan();
            // g->print();
            // g->calculateCopyNum();
            // g->checkLowerBound();
            // g->print();
            // g->calculateCopyNum();
            g->print();

            LocalGenomicMap * lgm = new LocalGenomicMap(g);
            
            // lgm->read_long_frags(longFragFn);
            // for (VertexPath *frag : *(lgm->get_long_frags())) {
            //     for (Vertex * v: *frag) {
            //         cout << v->getInfo() << " ";
            //     }
            //     cout << endl;
            // }
            // return 0;
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
                lgm->connectSourceSink();
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
        } catch (DuplicateJunctionException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (SegmentDoesNotExistException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (JunctionDoesNotExistException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (ILPBalancerInfeasibleException& e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (BackwardReachSourceNegativeException & e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (BackwardReachSinkPositiveException & e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (ForwardReachSinkNegativeException & e) {
            cout << e.what() << endl;
            return 1 ;
        } catch (ForwardReachSourcePositiveException & e) {
            cout << e.what() << endl;
            return 1 ;
            // } catch (...) {
            //     cout << "Logic error." << endl;
            //     return 1 ;
        }
        return 0;
    }
}
