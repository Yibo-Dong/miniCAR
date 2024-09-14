#include "SATInterface.h"
#include "cadical.hpp"
using namespace std;
using namespace CaDiCaL;
namespace car
{

    void CaDiCaLSolver::add_clause(const Cube &cl)
    {
        for (int id : cl)
            Solver::add(id);
        Solver::add(0); // end of adding a clause.
    }

    void CaDiCaLSolver::clear_assumption()
    {
        asm_store.clear();
    }

    void CaDiCaLSolver::push_assumption(const Cube &a)
    {
        asm_store.insert(asm_store.end(), a.begin(), a.end());
    }

    bool CaDiCaLSolver::solve_assumption()
    {
        for (auto lit : asm_store)
            Solver::assume(lit);
        int res = Solver::solve();
        if (res == 10)
            return true;
        else if (res == 20)
            return false;
        else
            assert(false && "should no happen\n");
    }

    Cube CaDiCaLSolver::get_uc()
    {
        Cube reason;
        for (int i = asm_store.size() - 1; i >= 0; --i)
        {
            const auto &lit = asm_store[i];
            if (Solver::failed(lit))
                reason.emplace_back(lit);
        }
        return reason;
    }

    Cube CaDiCaLSolver::get_model()
    {
        // FIXME: we don't really need so many literals. Shall we shrink here instead of later?
        Cube res(model_size, 0);
        for (int i = 1; i <= model_size; ++i)
        {
            if (Solver::val(i) > 0)
                res[i-1] = i;
            else
                res[i-1] = -i;
        }

        return res;
    }

    void CaDiCaLSolver::print_clauses(std::ostream &out)
    {
        // TODO: fill in here
    }

    void CaDiCaLSolver::print_assumption(std::ostream &out)
    {
        out << "assumptions in SAT solver: \n";
        if (!asm_store.size())
            out << " Empty ";
        for (int i = 0; i < asm_store.size(); ++i)
            out << asm_store[i] << " ";
        out << endl;
    }
}