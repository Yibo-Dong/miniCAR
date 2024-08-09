#ifndef SIMPSOLVER_H
#define SIMPSOLVER_H
#include "SimpSolver.h"
#include "SolverTypes.h"
#include "assert.h"

class SSLV: public Minisat::SimpSolver{
public:
    SSLV(){}
    inline Minisat::ClauseIterator clausesBegin() const { return Minisat::ClauseIterator(ca, &clauses[0]); }
    inline Minisat::ClauseIterator clausesEnd  () const { return Minisat::ClauseIterator(ca, &clauses[clauses.size()]); }
    inline Minisat::TrailIterator  trailBegin  () const { return Minisat::TrailIterator(&trail[0]); }
    inline Minisat::TrailIterator  trailEnd    () const { 
        return Minisat::TrailIterator(&trail[decisionLevel() == 0 ? trail.size() : trail_lim[0]]); }

};


#endif