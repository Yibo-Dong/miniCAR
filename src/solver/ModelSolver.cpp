#include "ModelSolver.h"
#include <vector>
using std::vector;

namespace car
{
    int FlagManager::MFlagOf(const int frame_level)
    {
        assert(frame_level >= 0);
        while (size_t(frame_level) >= MFlags.size())
        {
            int flag = max_flag++;
            MFlags.push_back(flag);
        }
        return MFlags.at(frame_level);
    }

    // backward mapping from Mflag to level.
    int FlagManager::MLevelOf(const int Mflag)
    {
        int abs_flag = abs(Mflag);
        for (size_t i = 0; i < MFlags.size(); ++i)
        {
            if (abs_flag == MFlags[i])
                return i;
        }
        return NOT_M_FLAG; // does not exists.
    }

    int FlagManager::getNewPTFlag()
    {
        int flag = max_flag++;
        PTFlags.push_back(flag);
        return flag;
    }

    int FlagManager::PFlagOf(const int frame_level)
    {
        assert(frame_level >= 0);
        while (size_t(frame_level) >= PFlags.size())
        {
            int flag = max_flag++;
            PFlags.push_back(flag);
        }
        return PFlags.at(frame_level);
    }

    // backward mapping from Mflag to level.
    int FlagManager::PLevelOf(const int PFlag)
    {
        int abs_flag = abs(PFlag);
        for (size_t i = 0; i < PFlags.size(); ++i)
        {
            if (abs_flag == PFlags[i])
                return i;
        }
        return NOT_P_FLAG; // does not exists.
    }

    int FlagManager::LevelOf(const int flag, SolverFlagEnum ftype)
    {
        switch (ftype)
        {
        case SolverFlag_Main:
            return MLevelOf(flag);
        case SolverFlag_Prop:
            return PLevelOf(flag);
        default:
        case SolverFlag_PropTemp:
            assert(false && "temp prop flag has no level\n");
        }
    }

    int FlagManager::FlagOf(const int level, SolverFlagEnum ftype)
    {
        switch (ftype)
        {
        case SolverFlag_Main:
            return MFlagOf(level);
        case SolverFlag_Prop:
            return PFlagOf(level);
        case SolverFlag_PropTemp:
            return -PTFlags.back();
        default:
            assert(false && "ftype not supproted\n");
        }
    }

    void ModelSolver::loadTrans()
    {
        for (const auto &cls : m->cls_)
            SATslv->add_clause(cls);
    }

    void ModelSolver::loadTransSimp()
    {
// block it temporarily
#if 0 
    int max_id = m->max_id();
    SSLV *sslv = new SSLV();
    for (size_t i = 0; i <= max_id; ++i)
    {
        sslv->newVar();
    }
    // minus 1 : we skip 0 in Aiger.
    for (size_t i = 1; i <= m->num_inputs(); ++i)
    {
        Var id = i - 1;
        sslv->setFrozen(id, true);
    }
    for (size_t i = m->num_inputs() + 1; i <= m->num_inputs() + m->num_latches(); ++i)
    {
        Var id = i;
        sslv->setFrozen(id - 1, true);
        sslv->setFrozen(abs(m->prime(id)) - 1, true);
    }
    // true && false
    sslv->setFrozen(abs(m->output(0)) - 1, true);
    sslv->setFrozen(abs(m->true_id()) - 1, true);
    sslv->setFrozen(abs(m->false_id()) - 1, true);
    for (int i = 0; i < m->size(); ++i)
    {
        auto &cl = m->element(i);
        vec<Lit> lits(cl.size());
        int index = 0;
        for (int id : cl)
        {
            assert(id != 0);
            int var = abs(id) - 1;
            Lit l = (id > 0) ? mkLit(var) : ~mkLit(var);
            lits[index++] = l;
        }
        bool res = sslv->addClause(lits);
    }
    sslv->eliminate(true);

    while (nVars() < sslv->nVars())
        newVar();

    for (Glucose::ClauseIterator c = sslv->clausesBegin(); c != sslv->clausesEnd(); ++c)
    {
        const Glucose::Clause &cls = *c;
        Glucose::vec<Glucose::Lit> cls_;
        for (int i = 0; i < cls.size(); ++i)
        {
            cls_.push(cls[i]);
        }
        SATslv->addClause(cls_);
    }

    for (auto c = sslv->trailBegin(); c != sslv->trailEnd(); ++c)
    {
        SATslv->addClause(*c);
    }
    delete sslv;
#endif
    }

    void ModelSolver::addNotCubeToLevel(const Cube &cu, int level, SolverFlagEnum ftype, bool primed)
    {
        assert(level > -1);
        int flag = FlagManager::FlagOf(level, ftype);
        vector<int> to_push = {-flag};
        to_push.reserve(1 + cu.size());

        // add negation of each lit
        for (size_t i = 0; i < cu.size(); ++i)
        {
            if (primed)
                to_push.push_back(-m->prime(cu[i]));
            else
                to_push.push_back(-cu[i]);
        }

        SATslv->add_clause(to_push);
    }

    bool ModelSolver::zeroStepReachable(const Cube &from, int target)
    {
        SATslv->clear_assumption();
        SATslv->push_assumption({target});
        SATslv->push_assumption(from);
        return SATslv->solve_assumption();
    }

    bool ModelSolver::oneStepReachable(const Cube &from, int level, bool forwardT)
    {
        assert(forwardT);
        SATslv->clear_assumption();
        activateLevel(level, SolverFlag_Main); // FIXME
        SATslv->push_assumption(from);

        return SATslv->solve_assumption();
    }

    void ModelSolver::activateLevel(int level, SolverFlagEnum ftype)
    {
        assert(level > -1);
        vector<int> to_push;
        auto flag = FlagManager::FlagOf(level, ftype);
        // activate one
        to_push.push_back(flag);
        // deactivate others
        switch (ftype)
        {
        case SolverFlag_Main:
            for (auto flg : FlagManager::MFlags)
            {
                if (flg != flag)
                    to_push.push_back(-flg);
            }
            break;

        case SolverFlag_Prop:
            for (auto flg : FlagManager::PFlags)
            {
                if (flg != flag)
                    to_push.push_back(-flg);
            }
            break;

        case SolverFlag_PropTemp:
        default:
            assert(false && "flag type not supported.");
        }
        SATslv->push_assumption(to_push);
    }

    State *ModelSolver::getState(bool shrink_to_previous)
    {
        Assignment model = SATslv->get_model();
        if (shrink_to_previous)
            shrinkModelToPrevious(model);
        Assignment inputs(model.begin(), model.begin() + m->num_inputs());
        Assignment latches(model.begin() + m->num_inputs(), model.begin() + m->num_inputs() + m->num_latches());
        State *s = new State(inputs, latches);
        return s;
    }

    void ModelSolver::shrinkModelToPrevious(Assignment &model)
    {
        Assignment res = model;
        int num_inputs = m->num_inputs();
        int num_latches = m->num_latches();

        for (int i = num_inputs + 1; i <= num_inputs + num_latches; i++)
        {
            int p = m->prime(i);
            assert(p != 0);

            int abs_p = abs(p);
            assert(model.size() >= size_t(abs_p));
            int val = model[abs_p - 1];
            if (p == val)
                res[i - 1] = i;
            else
                res[i - 1] = -i;
        }
        model = res;
    }
    Cube ModelSolver::getUCofLatch(bool shrink_to_previous)
    {
        Cube conflict = SATslv->get_uc();
        assert(!conflict.empty());

        if (shrink_to_previous)
            m->shrink_to_previous_vars(conflict);
        else
            m->shrink_to_latch_vars(conflict);

        return conflict;
    }
}