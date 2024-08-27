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

/**
 * @file mainsolver.h
 * @author yibodong (prodongf@gmail.com), jianwenli(lijwen2748@gmail.com)
 * @brief The main solver.
 * @version 0.1.0
 * 
 * 
 */

#ifndef NEW_MAIN_SOLVER_H
#define NEW_MAIN_SOLVER_H

#include "SimpSolver.h"
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
        bool reverseT;
        int bad;

	public:
        void loadSimpCNF();
		MainSolver(Problem *m, bool forward, bool rotate_is_on, bool uc_no_sort, bool simp);
		~MainSolver() {}

        /**
         * @brief check SAT(_s_ /\ ~p).
         * used in:
         * * initial condition
         * * last check in backward CAR(forward search), where ~p is represented by an apprixmation.
         */
        bool badcheck(const Assignment &st, const int bad);

        /// @brief  push the activation & deactivation flags. 
        /// @param  level
        void push_flags_M(const int);
        void push_to_assumption(const Cube&);
        void push_to_assumption_primed(const Cube&);
        /**
         * @brief check SAT(_s_  /\ T  /\ O_l' ) in backward CAR, or
         *              SAT(_s_' /\ T  /\ O_l  ) in forward CAR.
         */
		void set_assumption_M(State *s, const int frame_level);
		// same as ↑, with `prefers` flattened and placed in the front.
        void set_assumption_M(State *s, const int frame_level, const std::vector<Cube>& prefers);

        // set assumption, with latches primed, rest untouched.
        void set_assumption_primed(Cube& uc_or_flags);
		

        // here we could like to reuse the data-structure for phase-saving within Minisat, in order to guide it when SAT.
        // FIXME: cause false SAFE now.
        void set_expectation(const std::vector<int>& expectations);


		// @note: can only be used if the unroll level is 1. For unroll level>1, use `get_states`.
		State* get_state();
        State* getInitState();

		// this version is used for bad check only
		Cube get_conflict_no_bad(const int bad);
		Cube get_shrunk_uc();
        Cube get_conflict_another(int option, int nth);

		inline void add_cube_negate(const Cube &cu)
		{
			CARSolver::add_cube_negate(cu);
		}

		void add_new_frame_M(const OFrame& frame, const int frame_level);
        void add_new_frame_P(const OFrame& frame, const int frame_level);

        // used in main solver
		void add_clause_from_cube_M(const Cube &cu, const int frame_level);
        // used in propagation solver.
        void add_clause_from_cube_P(const Cube &cu, const int frame_level);

        // without the frame level, without prime.
        // if flag > 0, a real flag : also add the flag
        // if not, just ignore the flag.
        // FIXME: this is ugly.
		void add_clause_from_cube(const Cube &cu, int flag, bool __placeholder);

		void shrink_model(Assignment &model);

	public:
        /**
         * @section flags.
         * Categories of flags:
         * 1. O flags used in main solver.      ==> MFlag
         * 2. O flags use din propagation solver.   ==> PFlag
         * 3. propagation flags used in propagation test(to add temporary clauses) ==> PTFlag
         */
        

        // 1) O flags used in main solver. Used to (de)activate particular frames of O frame.
        // flags[i] = flagOf(O[i])
		std::vector<int> MFlags;

		inline int MFlagOf(const int frame_level)
		{
			assert(frame_level >= 0);
			while(size_t(frame_level) >= MFlags.size())
			{
                int flag = max_flag++;
				MFlags.push_back(flag);
			}
			return MFlags.at(frame_level);
		}

        // backward mapping from Mflag to level.
        inline int MLevelOf(const int Mflag)
        {
            int abs_flag = abs(Mflag);
            for(size_t i = 0; i< MFlags.size(); ++i)
            {
                if(abs_flag == MFlags[i])
                    return i;
            }
            return NOT_M_FLAG; // does not exists.
        }

        // 2) Temporary flags used in inductive checking(propagation).
        std::vector<int> PTFlag;
        inline int getNewPTFlag()
        {
            int flag = max_flag++;
            PTFlag.push_back(flag);
            return flag;
        }
        
        // 3) Propagation flags 
        std::vector<int> PFlags;

		inline int PFlagOf(const int frame_level)
		{
			assert(frame_level >= 0);
			while(size_t(frame_level) >= PFlags.size())
			{
                int flag = max_flag++;
				PFlags.push_back(flag);
			}
			return PFlags.at(frame_level);
		}

        // backward mapping from Mflag to level.
        inline int PLevelOf(const int PFlag)
        {
            int abs_flag = abs(PFlag);
            for(size_t i = 0; i< PFlags.size(); ++i)
            {
                if(abs_flag == PFlags[i])
                    return i;
            }
            return NOT_P_FLAG; // does not exists.
        }
        

    public: 
		inline int lits_per_round() const
		{
			// flag for each level, with the first layer's flag remaining unused
			return _model->max_id() + 1;
		}
	};

}

#endif
