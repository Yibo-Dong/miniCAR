#include "implysolver.h"
#include <algorithm>
using namespace std;

#ifdef MINISAT
using namespace Minisat;
#else
using namespace Glucose;
#endif

namespace car
{
    std::map<int, int> counter;
    std::map<int, std::map<int, double>> moms;
    std::map<int, std::shared_ptr<ImplySolver>> ImplySolver::imp_solvers;

    bool ImplySolver::is_blocked(State *s, int level)
    {
        // not so good to those with large latches scales.
        std::shared_ptr<ImplySolver> solver = getSolver(level);
        solver->set_assumption(s->get_latches());

        /* note, we set the budget of solver in this level.
           we do not actually need the solution. We just want to check whether it is blocked by any UC.
           if blocked, it's clear, that to propagate the state will already make the assignment of the UC's clause to be UNSAT.
        */
        solver->Solver::setPropBudget(s->num_latches_);
        return solver->Solver::solve_()==l_False;
    };

    void ImplySolver::add_uc(const std::vector<int> &uc, int level)
    {
        std::shared_ptr<ImplySolver> solver = getSolver(level);
        // add negate of uc as a clause.
        // no flag. forever added.
        vec<Lit> lits(uc.size());
        int index = 0;
        for(size_t id : uc)
            lits[index++] = solver->SAT_lit(-id);
        solver->addClause(lits);
    }
} // namespace car