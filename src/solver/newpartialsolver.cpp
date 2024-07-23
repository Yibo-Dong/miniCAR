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
	Update Date: October 6, 2017
	Main Solver in CAR
*/

#include "newpartialsolver.h"
#include "definition.h"
#include <algorithm>
using namespace std;
namespace car
{
	PartialSolver::PartialSolver (Problem* m) 
	{
		model_ = m;
		max_partial_id = m->max_id() +1;

	    //constraints
		for (int i = 0; i < m->outputs_start (); i ++)
		{
			add_clause (m->element (i));
		}
		//outputs
		for (int i = m->outputs_start (); i < m->latches_start (); i ++)
		{
			add_clause (m->element (i));
		}
		//latches
		for (int i = m->latches_start (); i < m->size (); i ++)
		{
			add_clause (m->element (i));
		}

	}
	
	void PartialSolver::add_clause_with_flag(Assignment& cls)
	{
		vector<int> cl = cls;
		int flag = max_partial_id;
		cls.push_back(flag);
		add_clause(cl);
	}

	Assignment PartialSolver::get_conflict()
	{
		Assignment conflict = get_uc();
		model_ -> shrink_to_latch_vars(conflict);
		std::sort (conflict.begin (), conflict.end (), car::absIncr);
		return conflict;
	}

	Assignment PartialSolver::get_conflict_no_bad(int bad)
	{
		Assignment conflict = get_uc_no_bad(bad);
		model_ -> shrink_to_latch_vars(conflict);

		std::sort (conflict.begin (), conflict.end (), car::absIncr);
		return conflict;
	}
	

	bool PartialSolver::solve_with_assumption(const Assignment& assumption)
	{
		set_assumption(assumption);
		bool res = solve_assumption ();
		return res;
	}
	
	void PartialSolver::set_assumption (const Assignment& assum)
	{
		assumptions.clear ();
		for (const int &var :assum)
		{
			assumptions.push (SAT_lit (var));
		}
	}
}
