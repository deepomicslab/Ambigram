#include <iostream>
#include "Exceptions.hpp"

using namespace std;

/* duplicate junction */
DuplicateJunctionException::DuplicateJunctionException(Junction *aJunction) {
    mJunction = aJunction;
    whatMsg = "DuplicateJunctionException: " + mJunction->getInfo()[0] + ", " + mJunction->getInfo()[1];
}

const char *DuplicateJunctionException::what() const throw() {
    return whatMsg.c_str();
}

/* cannot find segment */
SegmentDoesNotExistException::SegmentDoesNotExistException(int aSegId) {
    mSegId = aSegId;
    whatMsg = "SegmentDoesNotExistException: Segment with ID " + to_string(mSegId) + " does not exist";
}

const char *SegmentDoesNotExistException::what() const throw() {
    return whatMsg.c_str();
}

/* cannot find junction */
JunctionDoesNotExistException::JunctionDoesNotExistException(Edge *aEdge) {
    mEdge = aEdge;
    whatMsg = "JunctionDoesNotExistException: Junction with edge " + mEdge->getInfo() + " does not exist";
}

const char *JunctionDoesNotExistException::what() const throw() {
    return whatMsg.c_str();
}

/* backward reach source- */
BackwardReachSourceNegativeException::BackwardReachSourceNegativeException(Vertex *aStartVertex) {
    mStartVertex = aStartVertex;
    whatMsg = "BackwardReachSourceNegativeException: " + mStartVertex->getInfo() + " can be reached from source-.";
}

const char *BackwardReachSourceNegativeException::what() const throw() {
    return whatMsg.c_str();
}

/* backward reach sink+ */
BackwardReachSinkPositiveException::BackwardReachSinkPositiveException(Vertex *aStartVertex) {
    mStartVertex = aStartVertex;
    whatMsg = "BackwardReachSinkPositiveException: " + mStartVertex->getInfo() + " can be reached from sink+.";
}

const char *BackwardReachSinkPositiveException::what() const throw() {
    return whatMsg.c_str();
}

/* forward reach sink- */
ForwardReachSinkNegativeException::ForwardReachSinkNegativeException(Vertex *aStartVertex) {
    mStartVertex = aStartVertex;
    whatMsg = "ForwardReachSinkNegativeException: " + mStartVertex->getInfo() + " can reach sink-.";
}

const char *ForwardReachSinkNegativeException::what() const throw() {
    return whatMsg.c_str();
}

/* forward reach source+ */
ForwardReachSourcePositiveException::ForwardReachSourcePositiveException(Vertex *aStartVertex) {
    mStartVertex = aStartVertex;
    whatMsg = "ForwardReachSourcePositiveException: " + mStartVertex->getInfo() + " can reach source+.";
}

const char *ForwardReachSourcePositiveException::what() const throw() {
    return whatMsg.c_str();
}

/* duplicate record in database */
DuplicateRecordException::DuplicateRecordException(int aSourceId, char aSourceDir, int aTargetId, char aTargetDir) {
    mSourceId = aSourceId;
    mSourceDir = aSourceDir;
    mTargetId = aTargetId;
    mTargetDir = aTargetDir;

    whatMsg = "Duplicate junction record in database: "
              + to_string(mSourceId) + mSourceDir
              + "=>"
              + to_string(mTargetId) + mTargetDir;
}

const char *DuplicateRecordException::what() const throw() {
    return whatMsg.c_str();
}

/* ILP is infeasible */
const char *ILPBalancerInfeasibleException::what() const throw() {
    return "ILPBalancerInfeasibleException:: The balancer for the given graph is proven infeasible.";
}
