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

/*
	Author: Jianwen Li
	Update Date: October 23, 2017
	Start Solver in CAR
*/

#ifndef START_SOLVER_H
#define START_SOLVER_H

#include "carsolver.h"
#include "model.h"
#include "data_structure.h"
#include "statistics.h"

namespace car {
    extern Statistics CARStats;
    class StartSolver : public CARSolver {
    public:
        StartSolver (const Model* m, const int bad, const bool forward, const bool verbose = false)
        {
            if (!forward)
                add_cube (const_cast<Model*>(m)->init ());
            else
            {
                // FIXME: what does here mean?
                for (int i = 0; i < const_cast<Model*>(m)->latches_start (); i ++)
                    add_clause (const_cast<Model*>(m)->element (i));
                assumptions.push (SAT_lit (bad));
            }
            model_ = m;
            forward_ = forward;
            max_id_ = m->max_id () + 1;
            // NOTE: Because they are different solver objects, this flag will not conflict with flags in main solver.
            flag_ = max_id_;
            fresh = true;
        }
        ~StartSolver () {}
        
        inline bool solve_with_assumption ()
        {
            bool res = solve_assumption ();
		    return res;
        }
        
        // block last flag, and set a new flag.
        // NOTE: Together with "add_clause_with_flag()", this makes sure all the clauses inserted before are invalidated.
        inline void reset ()
        {
            if(fresh)
            {
                fresh = false;
                assumptions.push (SAT_lit (flag_));
                return;
            }
            assumptions.pop();
            assumptions.push(~SAT_lit (flag_));
            assumptions.push(SAT_lit (++flag_));
        }

        // This method is temporarily add, will be blocked later by reset().
        // To permanently add, use update_constraint() instead.
        inline void add_clause_with_flag (const Cube& cu)
        {
            std::vector<int> cl = {-flag_};
            for (int i = 0; i < cu.size (); i ++)
                cl.push_back (-cu[i]);
            add_clause (cl);
        }
        
        // This method is permanently add. 
        // To temporarily add, use add_clause_with_flag() instead.
        inline void update_constraint (Cube &cu)
        {
            CARSolver::add_cube_negate (cu);
        }
        inline bool get_direction() const
        {
            return forward_;
        }

        State* create_new_state()
        {
            int ninputs = model_->num_inputs();
            int nlatches = model_->num_latches();
            Assignment inputs(ninputs), latches(nlatches);
            for (int i = 0; i <ninputs ; i ++)
            {
                if (model[i] == l_True)
                    inputs[i] = i+1;
                else if (model[i] == l_False)
                    inputs[i] = -(i+1);
            }
            for (int i = 0; i < nlatches ; i ++)
            {
                if (model[i+ninputs] == l_True)
                    latches[i] = i+1+ninputs;
                else if (model[i+ninputs] == l_False)
                    latches[i] = -(i+1+ninputs);
            }
			State *res = new State(inputs,latches);
			return res;
        }
        
     private:
        int max_id_;
        int flag_;
        bool forward_;
        bool fresh; // whether it only contains initial assumptions.
        const Model* model_;
        
    };
}

#endif
