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

} // namespace car

#endif // !ImplySolver