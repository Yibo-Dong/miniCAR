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

#include "mainsolver.h"
#include "definition.h"
#include <algorithm>
using namespace std;
namespace car
{

	// record the flag of flags of every frams in every O sequence(in bi-car, there will be multiple O sequences, each starts from a unique center state.)
	std::unordered_map<Osequence*,std::vector<int>> MainSolver::flag_of_O;

	/**
	 * @brief 把constraints、outputs和latches加入clause。
	 * 		
	 * @param m 
	 * @param verbose 
	 */
	MainSolver::MainSolver (Model* m,int rotate_is_on, const bool verbose, bool uc_no_sort) : rotate_is_on(rotate_is_on), uc_no_sort(uc_no_sort)
	{
		model_ = m;
		max_flag = m->max_id() + 1;
		unroll_level = 1;
		// BASIC STEP:
		// (1) create clauses for constraints encoding
		for (int i = 0; i < m->common_next_start(); ++i)
		{
			add_clause (m->element (i));
		}
		// (2) same next have same previous, initialized to 0.
		for( int i = m->common_next_start(); i< m->outputs_start(); ++i)
		{
			add_clause(m->element(i));
		}
		// (3) clause for encoding outputs.
		for (int i = m->outputs_start (); i < m->latches_start (); ++i)
		{
			add_clause (m->element (i));
		}
		// (4) clause for encoding latches's mapping relation
		// (5) create clauses for true and false
		for (int i = m->latches_start (); i < m->size (); ++i)
		{
			add_clause (m->element (i));
		}
	}

	/**
	 * @brief Idea of unroll is to encode T^k into the solver with k > 1
	 * 
	 * Inputs are independant, output is just a mark. 
	 * Relation between different levels is built up upon Latches:
	 * 	The next value of prior cycle === The Previous value of next cycle
	 * 
	 * E.g.
	 * Input	: 2
	 * Latch	: (4, 2)
	 * Output	: 4
	 * AndGate	: {}
	 * 
	 * ---> 
	 * 
	 * Input	: 2, 1 * 4 + 2
	 * Latch	: (4, 2), (1 * 4 + 4, 1 * 4 + 2)
	 * Output	: 1 * 4 + 4
	 * AndGate	: {}
	 * 
	 * We need to record: 8 == 2.
	 * 
	 * To unroll one level:
	 * First, copy the literals. 
	 * 
	 * @param m 
	 * @param verbose 
	 * @param unroll_level 
	 */
	MainSolver::MainSolver (Model* m, const bool verbose, int unroll_level, bool __placeholder) 
	{
		assert(unroll_level >=1);
		// no need to unroll if level == 1.
		model_ = m;
		this -> unroll_level = 1;
		// note: add flags for each round
		int lits_each_round = lits_per_round();

		max_flag = unroll_level * lits_each_round;
			
		// BASIC STEP:
		// (1) create clauses for constraints encoding
		for (int i = 0; i < m->common_next_start(); ++i)
		{
			add_clause (m->element (i));
		}
		// (2) same next have same previous, initialized to 0.
		for( int i = m->common_next_start(); i< m->outputs_start(); ++i)
		{
			add_clause(m->element(i));
		}
		// (3) clause for encoding outputs.
		for (int i = m->outputs_start (); i < m->latches_start (); ++i)
		{
			add_clause (m->element (i));
		}
		// (4) clause for encoding latches's mapping relation
		// (5) create clauses for true and false
		for (int i = m->latches_start (); i < m->size (); ++i)
		{
			add_clause (m->element (i));
		}

		// flag for the first level
		int flag_for_first_level = 1 * lits_each_round;
		add_clause(flag_for_first_level);
		
		if(unroll_level == 1)
		{
			return;
		}

		// UNROLL SECTION:

		for(int level = 2; level <= unroll_level; ++level)
		{
			unroll();
		}
	}
	
	/**
	 * @brief assumptions = {assum, bad_id}
	 * 
	 * @param assum 
	 * @param bad_id 
	 */
	void MainSolver::set_assumption (const Assignment& assum, const int bad_id)
	{
        assumptions.clear ();
		assumptions.push (SAT_lit (bad_id));
		for (const int &var :assum)
		{
			assumptions.push (SAT_lit (var));
		}
	}

	/**
	 * @brief set assumption = { s->s() , flag_of(Os[frame_level]) }
	 * 
	 * @param s 
	 * @param frame_level 
	 * @param forward 
	 */
    void MainSolver::set_assumption(Osequence* O, State*s,  const int frame_level, const bool forward)
    {
        assumptions.clear();
        if (frame_level > -1)
			assumptions.push (SAT_lit (flag_of(O,frame_level)));
		for (const int &id :s->s())
		{
			int target = forward ? model_->prime (id) : id;
			assumptions.push (SAT_lit (target));
		}
    }

	void MainSolver::set_assumption(Osequence* O, State*s,  const int frame_level, const bool forward, const std::vector<Cube> & prefers)
    {
        assumptions.clear();
        if (frame_level > -1) // should always satisfy
			assumptions.push (SAT_lit (flag_of(O,frame_level)));
        for(int i = 0; i < prefers.size(); ++i)
		{
			auto &vec = prefers[i];

			for(auto &id : vec)
			{
				int target = forward?model_->prime(id) : id;
				assumptions.push(SAT_lit(target));
			}
		}
        if(!rotate_is_on)
        {
            // with rotate on, it will already be the whole state.
            // therefore, no need to add the rest part.
            for (const int &id :s->s())
            {
                int target = forward ? model_->prime (id) : id;
                assumptions.push (SAT_lit (target));
            }
        }
    }


	bool MainSolver::solve_with_assumption (const Assignment& st, const int p)
	{
		set_assumption(st,p);
		// FIXME: inline this. Anyway, speed it up!
		if(unroll_level > 1)
		{
			enable_level(unroll_level);
		}
        return solve_with_assumption();
	}


	//NOTE: this State* is not owned by this solver, but the checker. It should be immediately added into clear duty.	
	State* MainSolver::get_state (const bool forward)
	{
		assert(unroll_level == 1);
		Assignment model = get_model();
		if(!forward)
			shrink_model (model);
		Assignment inputs(model.begin(),model.begin() + model_->num_inputs());
		Assignment latches(model.begin() + model_->num_inputs(),model.begin() + model_->num_inputs()+model_->num_latches());
		State* s = new State(inputs,latches);
		return s;
	}

	// the states in ret is not owned by MainSolver
	void MainSolver::get_states(std::vector<State*>& ret, const bool forward)
	{
		Assignment model = get_model();
		assert(model.size() == unroll_level * lits_per_round());
		for(int level = 0 ; level < unroll_level; ++level)
		{
			int offset = level * lits_per_round();

			Assignment model_for_this_level(offset + model.begin(), offset + model.begin() + model_->num_inputs() + model_->num_latches());
			//TODO: this needs to be tested.
			if(forward)
				shrink_model(model_for_this_level);
			for(auto &i : model_for_this_level)
			{
				// the behavior of negative number's modulo arithmatic is not well defined. do not use. 
				i -= ( i > 0 ? 1 : -1) * offset;
			}
			Assignment inputs(model_for_this_level.begin(), model_for_this_level.begin() + model_->num_inputs());
			Assignment latches(model_for_this_level.begin() + model_->num_inputs(), model_for_this_level.begin() + model_->num_inputs() + model_->num_latches());
			State* s = new State(inputs,latches);
			ret.push_back(s);
		}
	}

	Assignment MainSolver::get_state_full_assignment (const bool forward)
	{
		Assignment model = get_model();
		if(!forward)
			shrink_model (model);
		return std::move(model);
	}
	
	/**
	 * @brief get unsat core, and remove bad from it if exists.
	 * 
	 * @note this version is used for bad check only
	 * @note the last lit is pushed again.
	 * @param bad 
	 * @return Cube 
	 */
	Cube MainSolver::get_conflict_no_bad (const int bad)
	{
		Cube res = get_uc_no_bad (bad);
		if(res.empty())
			return res;
        if(!uc_no_sort)
        {
    		sort(res.begin(),res.end(),car::absIncr);
        }
		return std::move(res);
	}
	
	/**
	 * @brief get unsat core, get latch, shrink to previous if forward. 
	 * @note the last lit is pushed again.
	 * @param forward 
	 * @return Cube 
	 */
	Cube MainSolver::get_conflict (const bool forward)
	{
		Cube conflict = get_uc ();
        if(model_->latch_var(conflict.back()))
        {
            cout<<"This is interesting"<<endl;
            cout<<"State is";
            for(auto l: get_assumption())
            {
                cout<<l<<", ";
            }
            cout<<"uc is ";
            for(auto l:conflict)
            {
                cout<<l<<", ";
            }
        }

		assert(!conflict.empty());
		if (forward)
		{
		    model_->shrink_to_previous_vars (conflict);
		}
		else
		{
		    model_->shrink_to_latch_vars (conflict);
		}
        if(conflict.empty())
        {
            // should not enter. For debug only.
            cout<<"error. conflict is empty."<<endl;
            Cube conf = get_uc();
            for(auto l:conf)
                cout<<l<<", ";
            cout<<endl;
        }
		int conflict_back = conflict.back();

        if(!uc_no_sort)
        {
            std::sort (conflict.begin (), conflict.end (), car::absIncr);
        }
		return std::move(conflict);
	}

    Cube MainSolver::get_conflict_another (const bool forward,int option, int nth)
	{
		Cube conflict = get_uc_another (option, nth);
		assert(!conflict.empty());
		if (forward)
		{
		    model_->shrink_to_previous_vars (conflict);
		}
		else
		{
		    model_->shrink_to_latch_vars (conflict);
		}
		int conflict_back = conflict.back();

        #ifndef UC_NO_SORT
		std::sort (conflict.begin (), conflict.end (), car::absIncr);
        #endif
		return std::move(conflict);
	}
	
	void MainSolver::add_new_frame(const Frame& frame, const int frame_level,Osequence *O, const bool forward)
	{
		for (int i = 0; i < frame.size (); i ++)
		{
			add_clause_from_cube (frame[i], frame_level, O, forward);
		}
	}

	// Actually, each center_state should only exist in O sequence in one direction.
	// Should we just record this?
	void MainSolver::add_clause_from_cube(const Cube &cu, const int frame_level, Osequence *O, const bool forward)
	{
		int flag = flag_of(O,frame_level);
		vector<int> cl = {-flag};
		for (int i = 0; i < cu.size (); i ++)
		{
			if (!forward)
				cl.push_back (-model_->prime (cu[i]));
			else
				cl.push_back (-cu[i]);
		}
		add_clause (cl);
	}

	void MainSolver::shrink_model(Assignment& model)
	{
		Assignment res = model;

		int num_inputs = model_->num_inputs();
		int num_latches = model_->num_latches();
		
		for (int i = num_inputs + 1; i <= num_inputs + num_latches; i++)
		{
			int p = model_->prime(i);
			assert(p != 0);
			
			int abs_p = abs(p);
			assert(model.size() >= abs_p);
			
			int val = model[abs_p - 1];
			
			if (p == val)
				res[i - 1] = i;
			else
				res[i - 1] = -i;
		}

		model = res;
	}

	void MainSolver::unroll()
	{
		int level = unroll_level+1;
		max_flag = level * lits_per_round();
		int lits_each_round = lits_per_round();
		int flag_for_this_level = level * lits_each_round;
		// copy clauses.
		// @note: last 2: true & false. The two lits reserved.
		for (int i = 0; i < model_->size(); ++i)
		{
			// copy it first
			vector<int> unrolled_clause = model_->element(i);
			
			// unroll it
			for(int &lit : unrolled_clause)
			{
				lit += (lit > 0 ? 1 : -1) *  (level-1) * lits_each_round;
			}
			// add flag
			unrolled_clause.push_back(-flag_for_this_level);

			add_clause(unrolled_clause);
		}

		// add equivalent relation between levels
		// prior's next is present!
		for (int i = model_->num_inputs() + 1; i < model_->num_inputs() + model_->num_latches() + 1; ++i)
		{
			int present = (level - 1) * lits_each_round + i;
			int prior_next = (level - 2) * lits_each_round * (model_->next_map_[i] > 0 ? 1 : -1) + model_->next_map_[i];
			// equiv
			add_clause(-present,prior_next,-flag_for_this_level);
			add_clause(present,-prior_next,-flag_for_this_level);
		}
		unroll_level = level;
	}

    // deprecated!
	void MainSolver::enable_level(int level)
	{
		assert(level <= unroll_level);
		for(int i = 2; i <= level; ++i)
		{
			int flag_for_this_level = i * lits_per_round();
			assumptions.push (SAT_lit (flag_for_this_level));
		}
	}

    // deprecated!
	void MainSolver::enable_max()
	{
		for(int i = 1; i <= unroll_level; ++i)
		{
			int flag_for_this_level = i * lits_per_round();
			assumptions.push (SAT_lit (flag_for_this_level));
		}
	}
	
}
