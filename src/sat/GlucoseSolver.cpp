#include "SATInterface.h"
#include "Solver.h"
#include "SolverTypes.h"
using namespace std;
using namespace Glucose;
namespace car{
    Glucose::Lit GlucoseSolver::SAT_lit(int id)
    {
        assert(id != 0); // only +var and -var.
        int var = abs(id) - 1;
        while (nVars() <= var)
        {
            newVar();
        }
        return mkLit(var, (id < 0));
    }

    int GlucoseSolver::lit_id(Glucose::Lit l) const
    {
        if (sign(l))
            return -(var(l) + 1);
        else
            return var(l) + 1;
    }

    void GlucoseSolver::add_clause(const Cube &cl)
    {
        int index = 0;
        vec<Lit> lits(cl.size());
        for (int id : cl)
            lits[index++] = SAT_lit(id);
        addClause(lits);
    }

    void GlucoseSolver::clear_assumption()
    {
        assumptions.clear();
    }

    void GlucoseSolver::push_assumption(const Cube &a)
    {
        for(auto lit:a)
            assumptions.push(SAT_lit(lit));
    }

    bool GlucoseSolver::solve_assumption()
    {
        return solve_() == l_True;
    }

    Cube GlucoseSolver::get_uc()
    {
        Cube reason(conflict.size(), 0);
        for (int k = 0; k < conflict.size(); ++k)
        {
            reason[k] = -lit_id(conflict[k]);
        }
        return reason;
    }

    Cube GlucoseSolver::get_model()
    {
        Cube res(nVars(), 0);
        for (int i = 0; i < nVars(); ++i)
        {
            res[i] = (model[i] == l_True) ? i + 1 : (-(i + 1));
        }
        return res;
    }

    void GlucoseSolver::print_clauses(std::ostream &out)
    {
        out << "clauses in SAT solver: \n";
        for (int i = 0; i < clauses.size(); ++i)
        {
            Glucose::Clause &c = ca[clauses[i]];
            for (int j = 0; j < c.size(); ++j)
                out << lit_id(c[j]) << " ";
            out << "0 " << endl;
        }
    }

    void GlucoseSolver::print_assumption(std::ostream &out)
    {
        out << "assumptions in SAT solver: \n";
        if (!assumptions.size())
            out << " Empty ";
        for (int i = 0; i < assumptions.size(); ++i)
            out << lit_id(assumptions[i]) << " ";
        out << endl;
    }
}