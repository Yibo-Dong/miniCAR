#ifndef BMCCHECKER_H

#define BMCCHECKER_H

#include "definition.h"
#include "mainsolver.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
using namespace car; // This is strange! remember to change later pls!


/**
 * @file carChecker.h
 * @author yibodong (prodongf@gmail.com)
 * @brief A checker, using Bounded Model Checking algorithm.
 * @version 0.1.0
 * 
 * 
 */

namespace bmc{
    // BMC is just: keep unrolling until reach the bad state.(Sure, may not stop on safe cases.)
    class BMCChecker
    {
        public:
            explicit BMCChecker(Problem *model) : model_(model), lev(0)
            {
                // FIXME: should not be true, ture here. change it.
                solver = new MainSolver(model,false, true, true, 1);
                // the initial state
                init = new State(model->init());
                cex.clear();
            }

            BMCChecker(const BMCChecker&)=delete;
            
            const BMCChecker& operator=(const BMCChecker&)=delete;

            ~BMCChecker()
            {
                delete init;
                delete solver;
                for(auto s: cex)
                    if(s)
                        delete s;
                cex.clear();

            }
            // FIXME: rewrite trivial check
            void check();
            void checkLimited(int limit);

            void printEvidence(std::ostream& res_file);

            inline int getLevel(){return lev;}
        private:
            Problem* model_;
            MainSolver *solver;
            int lev = 0;
            State *init;
            std::vector<car::State *> cex;
    };
}





#endif