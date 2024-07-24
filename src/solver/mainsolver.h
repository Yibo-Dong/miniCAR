/*
	Copyright (C) 2018, Jianwen Li (lijwen2748@gmail.com), Iowa State University

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef NEW_MAIN_SOLVER_H
#define NEW_MAIN_SOLVER_H

#include "carsolver.h"
#include "definition.h"
#include "statistics.h"
#include <vector>
#include <assert.h>
#include <iostream>

namespace car
{
	class MainSolver : public CARSolver
	{
    public:
        const Problem *_model;
        const bool rotate_is_on;
        const bool uc_no_sort;
		int max_flag;
		int unroll_level;

	public:
		MainSolver(Problem *m, bool rotate_is_on, bool uc_no_sort);
		MainSolver(Problem* m, bool rotate_is_on, bool uc_no_sort, int unroll_level) ;
		~MainSolver() {}

        /**
         * @brief check SAT(_s_ /\ ~p).
         * used in:
         * * initial condition
         * * last check in backward CAR(forward search), where ~p is represented by an apprixmation.
         */
		void set_assumption(const Assignment &st, const int bad);

        /**
         * @brief check SAT(_s_  /\ T  /\ O_l' ) in backward CAR, or
         *              SAT(_s_' /\ T  /\ O_l  ) in forward CAR.
         */
		void set_assumption(OSequence *O, State *s, const int frame_level, const bool forward);
		// same as ↑, with `prefers` flattened and placed in the front.
        void set_assumption(OSequence *O, State *s, const int frame_level, const bool forward, const std::vector<Cube>& prefers);
		
        // if assumptions are already set, just solve.
        inline bool solve_with_assumption() { return CARSolver::solve_assumption(); }
        // otherwise, set assumption, then solve.
        bool solve_with_assumption(const Assignment &assump, const int bad);

        // here we could like to reuse the data-structure for phase-saving within Minisat, in order to guide it when SAT.
        // FIXME: cause false SAFE now.
        inline void set_expectation(const std::vector<int>& expectations, const bool forward)
        {
            assert(!forward && "temporarily not avaiable for forward-car");
            for(auto& id : expectations)
            {
                bool sgn = (id > 0 ? true : false);
                int var = sgn ? id : -id;
                // see whether they are state literals
                assert(_model->latch_var(var));
                // get the prime version of the literals
                auto id_prime = _model->prime(var);
                // it shall already be encoded.
                assert(nVars() >= id_prime);
                
                // set the polarity. It will work by default.
                if(forward)
                    setPolarity(var, sgn);
                else
                    setPolarity(id_prime, sgn);
            }
            
        }

		// @note: can only be used if the unroll level is 1. For unroll level>1, use `get_states`.
		State* get_state(const bool forward);

		Assignment get_state_full_assignment(const bool forward);

		// this version is used for bad check only
		Cube get_conflict_no_bad(const int bad);
		Cube get_conflict(const bool forward);
        Cube get_conflict_another(const bool forward,int option, int nth);

		inline void add_cube_negate(Cube &cu)
		{
			CARSolver::add_cube_negate(cu);
		}

		void add_new_frame(const OFrame& frame, const int frame_level, OSequence *O, const bool forward);

		void add_clause_from_cube(const Cube &cu, const int frame_level, OSequence *O, const bool forward);

		void shrink_model(Assignment &model);

	public:
		std::vector<int> flags;

		inline int flag_of(OSequence *o, const int frame_level)
		{
			assert(frame_level >= 0);
			while(frame_level >= flags.size())
			{
				flags.push_back(max_flag++);
			}
			return flags.at(frame_level);
		}
		


    public: 
        // unroll for BMC:
        void enable_level(int level);
		void enable_max();

		inline int lits_per_round() const
		{
			// flag for each level, with the first layer's flag remaining unused
			return _model->max_id() + 1;
		}
		
        // NOTE: this is deprecated for now, because we only pop one frame flag at present.
		void unroll();

        void get_states(std::vector<State*>& states, const bool forward);
	};

}

#endif
