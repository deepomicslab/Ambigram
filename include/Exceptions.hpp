#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

#include "Junction.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"

using namespace std;

class Junction;

class Vertex;

class Edge;

class DuplicateJunctionException : public exception {
private:
    Junction *mJunction;
    string whatMsg;
public:
    DuplicateJunctionException(Junction *aJunction);

    virtual const char *what() const throw();
};

class SegmentDoesNotExistException : public exception {
private:
    int mSegId;
    string whatMsg;
public:
    SegmentDoesNotExistException(int aSegId);

    virtual const char *what() const throw();
};

class JunctionDoesNotExistException : public exception {
private:
    Edge *mEdge;
    string whatMsg;
public:
    JunctionDoesNotExistException(Edge *aEdge);

    virtual const char *what() const throw();
};

class BackwardReachSourceNegativeException : public exception {
private:
    Vertex *mStartVertex;
    string whatMsg;
public:
    BackwardReachSourceNegativeException(Vertex *aStartVertex);

    virtual const char *what() const throw();
};

class BackwardReachSinkPositiveException : public exception {
private:
    Vertex *mStartVertex;
    string whatMsg;
public:
    BackwardReachSinkPositiveException(Vertex *aStartVertex);

    virtual const char *what() const throw();
};

class ForwardReachSinkNegativeException : public exception {
private:
    Vertex *mStartVertex;
    string whatMsg;
public:
    ForwardReachSinkNegativeException(Vertex *aStartVertex);

    virtual const char *what() const throw();
};

class ForwardReachSourcePositiveException : public exception {
private:
    Vertex *mStartVertex;
    string whatMsg;
public:
    ForwardReachSourcePositiveException(Vertex *aStartVertex);

    virtual const char *what() const throw();
};

class DuplicateRecordException : public exception {
private:
    int mSourceId, mTargetId;
    char mSourceDir, mTargetDir;
    string whatMsg;

public:
    DuplicateRecordException(int aSourceId, char aSourceDir, int aTargetId, char aTargetDir);

    virtual const char *what() const throw();
};

class ILPBalancerInfeasibleException : public exception {
public:
    virtual const char *what() const throw();
};


#endif
