#ifndef ImplySolver_H
#define ImplySolver_H

/**
 * @file implysolver.h
 * @author yibodong (prodongf@gmail.com)
 * @brief A variant of solver, for checking implication.
 * @version 0.1.0
 * 
 * 
 */

#include "carsolver.h"
#include "definition.h"
#include "statistics.h"
#include <map>
#include <memory>
#include <vector>

namespace car
{
    extern Statistics CARStats;
    class ImplySolver : public CARSolver
    {
    public:
        explicit ImplySolver(int lvl) : level(lvl){};
        ~ImplySolver(){};

        static std::shared_ptr<ImplySolver> getSolver(int lvl)
        {
            if (imp_solvers.find(lvl) == imp_solvers.end())
            {
                // If solver for the level doesn't exist, create a new one
                imp_solvers[lvl] = std::make_shared<ImplySolver>(lvl);
            }
            return imp_solvers[lvl];
        }

        static bool is_blocked(State *s, int level);

        static void add_uc(const std::vector<int> &uc, int level);

        static inline void print_sz() {
            std::cerr<<"start printing sz"<<std::endl;
            for(auto& pr: imp_solvers) { std::cerr<<pr.first<<" : "<<pr.second->nClauses() <<std::endl; } 
            std::cerr<<"end printing sz"<<std::endl;
        }

        static void reset_all()
        {
            for(auto& pr: imp_solvers) { 
                pr.second = nullptr;
            }
            imp_solvers.clear();
        }

    private:
        
        inline void set_assumption(const std::vector<int> &state)
        {
            assumptions.clear();
            for (const int &var : state)
            {
                assumptions.push(SAT_lit(var));
            }
        }
        // which level it is for.
        int level;
        // Static map to store solvers for each level
        static std::map<int, std::shared_ptr<ImplySolver>> imp_solvers;
    };

    class UCMask
    {
        bool isState = false;

    public:
        std::vector<uint64_t> value;
        std::vector<uint64_t> mask; // uc only.
    public:
        UCMask(const Cube &uc, int nlatches, int startPos) : isState(false)
        {
            // 64 bit per Uint64.
            size_t nUInts = (nlatches + 0x3f) >> 6;
            value.resize(nUInts, 0);
            mask.resize(nUInts, 0);

            for (const auto &lit : uc)
            {
                // Which latch it is.
                int offset = abs(lit) - startPos;
                unsigned char Value = lit > 0 ? 0b1 : 0b0;
                assert(offset >= 0);
                // each UInt contains 64 lits. Mod 64 to see which UInt it is at.
                size_t whichUInt = offset >> 6;
                // % 64, get its position within this UInt. From Low To High.
                size_t IndexInUInt = (offset & 0x3f);
                uint64_t toadd = ((uint64_t)1 << IndexInUInt);
                mask[whichUInt] |= toadd;
                if (lit > 0)
                {
                    value[whichUInt] |= toadd;
                }
            }
        }

        UCMask(const Cube &state) : isState(true)
        {
            assert(!state.empty());
            int startPos = abs(state[0]);

            size_t nUInts = (state.size() + 0x3f) >> 6;
            value.resize(nUInts, 0);

            for (const auto &lit : state)
            {
                // Which latch it is.
                int offset = abs(lit) - startPos;
                unsigned char Value = lit > 0 ? 0b1 : 0b0;
                assert(offset >= 0);
                // each UInt contains 64 lits. Mod 64 to see which UInt it is at.
                size_t whichUInt = offset >> 6;
                // % 64, get its position within this UInt. From Low To High.
                size_t IndexInUInt = (offset & 0x3f);
                if (lit > 0)
                {
                    uint64_t toadd = ((uint64_t)1 << IndexInUInt);
                    value[whichUInt] |= toadd;
                }
            }
        }

        bool imply(const UCMask &another) const
        {
            assert(true == isState && false == another.isState);
            for (size_t i = 0; i < another.mask.size(); ++i)
            {
                if (((value[i] ^ another.value[i]) & another.mask[i]))
                    return false;
            }
            return true;
        }
    };

} // namespace car

#endif // !ImplySolver