#ifndef BMCCHECKER_H

#define BMCCHECKER_H

#include "data_structure.h"
#include "mainsolver.h"
#include "model.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
using namespace car; // This is strange! remember to change later pls!

namespace bmc{
    // BMC is just: keep unrolling until reach the bad state.(Sure, may not stop on safe cases.)
    class BMCChecker
    {
        public:
            BMCChecker(Model *model) : model_(model), lev(0)
            {
                solver = new MainSolver(model, false, 1, true);
                // the initial state
                init = new State(model->init());
                cex.clear();
            }

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
            Model* model_;
            MainSolver *solver;
            int lev = 0;
            State *init;
            std::vector<car::State *> cex;
    };
}





#endif