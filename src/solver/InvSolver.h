#ifndef NEW_INV_SOLVER_H
#define NEW_INV_SOLVER_H

#include "definition.h"
#include "SATInterface.h"
#include <unordered_map>

namespace car
{
    // simple rewriten. Rethink this please.
    class InvSolver
    {
    private:
        int max_flag;
        SATSolver *SATslv;
        std::unordered_map<int, int> pos_flags, neg_flags; // level -> flag

    public:
        InvSolver(const Problem *m)
        {
            SATslv = SATSolverFactory::create_solver("Glucose");
            loadConstraints(m);
            max_flag = m->max_id() + 1;
        }
        ~InvSolver() { delete SATslv; }

        inline void loadConstraints(const Problem *m)
        {
            for (int i = 0; i < m->outputs_start(); ++i)
                SATslv->add_clause(m->element(i));
        }

        bool invFound(const OSequence &O, int fresh_level)
        {
            // TODO: we could reuse it, so we do not care the fresh level when loading clauses.
            updateClauses(O);
            for (size_t i = fresh_level; O.size() > i + 1; ++i)
                if (invFoundAt(i))
                    return true;
            return false;
        }

    private:
        void updateClauses(const OSequence &O)
        {
            // if empty, no need to load.
            if (1 > O.size())
                return;
            // O[0], will only be used as neg.
            updateFrameNegative(O[0], 0);
            for (size_t i = 1; O.size() > i; ++i)
            {
                updateFrameNegative(O[i], i);
                updateFramePositive(O[i], i);
            }
        }

        // Inv(k) ::=
        // |= O[k+1] \in  Union(O[0], ... O[k]).
        // |= O[k+1] -> ( O[0] \/ ... O[k])
        // UNSAT( !(O[k+1] -> ( O[0] \/ ... O[k])) )
        // UNSAT( !(!O[k+1] \/ O[0] \/ ... O[k] ) )
        // UNSAT( O[k+1] /\ !O[k] /\ ... /\ !O[0])
        bool invFoundAt(int level)
        {
            SATslv->clear_assumption();
            // deactivating
            for (size_t i = level + 1;; ++i)
            {
                if (neg_flags.find(i) != neg_flags.end())
                    SATslv->push_assumption({-neg_flags[i]});
                else
                    break;
            }
            for (size_t i = 1; i <= level; ++i) // no pos 0.
            {
                SATslv->push_assumption({-pos_flags[i]});
            }
            for (size_t i = level + 2;; ++i)
            {
                if (pos_flags.find(i) != pos_flags.end())
                    SATslv->push_assumption({-pos_flags[i]});
                else
                    break;
            }

            // activating
            for (size_t i = 0; i <= level; ++i)
            {
                assert(neg_flags.find(i) != neg_flags.end());
                SATslv->push_assumption({neg_flags[i]});
            }
            assert(pos_flags.find(level + 1) != pos_flags.end());
            SATslv->push_assumption({pos_flags[level + 1]});

            return !SATslv->solve_assumption();
        }

        // OFrame ::= !uc0 /\ !uc1 /\ ... /\ !uck
        void updateFramePositive(const OFrame &f, int level)
        {
            if (pos_flags.find(level) == pos_flags.end())
                pos_flags[level] = max_flag++;
            int frame_flag = pos_flags[level];
            for (size_t i = 0; f.size() > i; ++i)
            {
                Cube vec = f[i];
                for (auto &lit : vec)
                {
                    lit = -lit;
                }
                vec.push_back(-frame_flag);
                SATslv->add_clause(vec);
            }
        }

        // !OFrame ::= uc0 \/ uc1 \/ ... \/ uck
        // uc_flag0 <-> uc0
        // ...
        // uc_flag0 \/ uc_flag1 \/ ... \/ uc_flagk

        // uc_flag0 -> uc0  /\ uc0 -> uc_flag0
        // !uc_flag0 \/ uc0  /\ !uc0 \/ uc_flag0
        void updateFrameNegative(const OFrame &f, int level)
        {
            if (neg_flags.find(level) == neg_flags.end())
                neg_flags[level] = max_flag++;
            int frame_flag = neg_flags[level];

            Cube uc_flags = {-frame_flag};
            for (size_t i = 0; f.size() > i; ++i)
            {
                auto &uc = f[i];
                int uc_flag = max_flag++;
                uc_flags.push_back(uc_flag);

                // !uc_flag0 \/ uc0
                for (auto &lit : uc)
                {
                    SATslv->add_clause({-frame_flag, -uc_flag, lit});
                }

                // !uc0 \/ uc_flag0
                Cube vec = uc;
                for (auto &lit : vec)
                    lit = -lit;
                vec.push_back(uc_flag);
                vec.push_back(-frame_flag);
                SATslv->add_clause(vec);
            }
            // uc_flag0 \/ uc_flag1 \/ ... \/ uc_flagk
            SATslv->add_clause(uc_flags);
        }
    };

}

#endif