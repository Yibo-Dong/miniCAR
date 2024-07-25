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
 * @file carsolver.h
 * @author jianwenli(lijwen2748@gmail.com), yibodong (prodongf@gmail.com)
 * @brief A wrapper outside SAT solver.
 * @version 0.1.0
 * 
 * 
 */
 
#ifndef CAR_SOLVER_H
#define	CAR_SOLVER_H

#include "statistics.h"
#include <iostream>
#include "Solver.h"
#include <vector>

namespace CARSolverNS {
    #ifdef MINISAT
		using SolverType = Minisat::Solver;
    #else
		using SolverType = Glucose::Solver;
    #endif // MINISAT
}

namespace car
{
		extern Statistics CARStats;
		class CARSolver : public CARSolverNS::SolverType
		{
		public:
			CARSolver() {}

#ifdef MINISAT
			Minisat::Lit SAT_lit(int id);	        // create the Lit used in SAT solver for the id.
			int lit_id(Minisat::Lit) const;         // return the id of SAT lit
#else
			Glucose::Lit SAT_lit(int id);			// create the Lit used in SAT solver for the id.
			int lit_id(Glucose::Lit) const;			// return the id of SAT lit
#endif 
            /**
             * @brief Solve with the assumptions in _assumption. 
             * @note before this, make sure all the assumption lits are put into assumptions.
             */
			inline bool solve_assumption() { auto res = solve_(); return res == l_True; }								// Solve with the assumptions in _assumption.
			inline void clear_assumption() { assumptions.clear(); } // clear the assumptions
			inline int size() { return clauses.size(); }			// Solve with the assumptions in _assumption.
			std::vector<int> get_assumption() const;				// get the assumption
			std::vector<int> get_model() const;						// get the model from SAT solver
			std::vector<int> get_uc() const;						// get UC from SAT solver
            std::vector<int> get_uc_rand();						    // get UC from SAT solver
            std::vector<int> get_uc_fp_rev();                       // get UC from SAT solver
            std::vector<int> get_uc_fp();                           // get UC from SAT solver
			std::vector<int> get_uc_no_bad(int bad) const;			// get UC from SAT solver
            std::vector<int> get_uc_another(int option, int nth);	// get another UC from SAT solver
			void add_cube(const std::vector<int> &);				// add each element in uc as a clause
			void add_cube_negate(const std::vector<int> &cu);	    // add the negate of the cube
			void add_clause_internal(const std::vector<int> &);		// add the or clause

			template <typename... Args>
			void add_clause(Args... args)
			{
				std::vector<int> v = {args...};
				add_clause_internal(v);
			}

			// printers
			void print_clauses(std::ostream & out_stream);
			void print_assumption(std::ostream & out_stream);
		};
}

#endif
