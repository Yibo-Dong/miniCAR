#ifndef SIMPSOLVER_H
#define SIMPSOLVER_H
#include "SimpSolver.h"
#include "SolverTypes.h"
#include "assert.h"

class SSLV: public Glucose::SimpSolver{
public:
    SSLV(){}
    inline Glucose::ClauseIterator clausesBegin() const { return Glucose::ClauseIterator(ca, &clauses[0]); }
    inline Glucose::ClauseIterator clausesEnd  () const { return Glucose::ClauseIterator(ca, &clauses[clauses.size()]); }
    inline Glucose::TrailIterator  trailBegin  () const { return Glucose::TrailIterator(&trail[0]); }
    inline Glucose::TrailIterator  trailEnd    () const { 
        return Glucose::TrailIterator(&trail[decisionLevel() == 0 ? trail.size() : trail_lim[0]]); }
};


#endif