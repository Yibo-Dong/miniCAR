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
#include "SSLV.h"
#include <algorithm>
using namespace std;
using namespace Glucose;
namespace car
{
	/**
	 * @brief put all the clauses generated into SAT solver.
	 * 		
	 */
    MainSolver::MainSolver(Problem *m, bool forward, bool rotate_is_on, bool uc_no_sort, bool simp) : rotate_is_on(rotate_is_on), uc_no_sort(uc_no_sort), _model(m), unroll_level(1), reverseT(forward), bad(m->output(0))
    {
        max_flag = m->max_id() + 1;
        if(!simp)
        {            
            // (1) create clauses for constraints encoding
            for(int i = 0; i < m->common_next_start(); ++i)
            {
                add_clause(m->element(i));
            }
            // (2) same next have same previous, initialized to 0.
            for(int i = m->common_next_start(); i < m->outputs_start(); ++i)
            {
                add_clause(m->element(i));
            }
            // (3) clause for encoding outputs.
            for(int i = m->outputs_start(); i < m->latches_start(); ++i)
            {
                add_clause(m->element(i));
            }
            // (4) clause for encoding latches's mapping relation
            // (5) create clauses for true and false
            for(int i = m->latches_start(); i < m->size(); ++i)
            {
                add_clause(m->element(i));
            }
        }
        else{
            assert(_model->num_constraints() == 0);
            loadSimpCNF();
        }
    }

    void MainSolver::loadSimpCNF()
    {
        auto m = _model;
        int max_id = m->max_id();
        SSLV *sslv = new SSLV();
        for(size_t i = 0; i <= max_id; ++i)
        {
            sslv->newVar();
        }
        // minus 1 : we skip 0 in Aiger.
        for(size_t i = 1; i <= m->num_inputs(); ++i)
        {
            Var id = i - 1; 
            sslv->setFrozen(id, true);
        }
        for(size_t i = m->num_inputs() + 1; i <= m->num_inputs() + m->num_latches(); ++i)
        {
            Var id = i;
            sslv->setFrozen(id - 1, true);
            sslv->setFrozen(abs(m->prime(id)) - 1, true);
        }
        // true && false
        sslv->setFrozen(abs(m->output(0)) - 1, true);
        sslv->setFrozen(abs(m->true_id()) - 1, true);
        sslv->setFrozen(abs(m->false_id()) - 1, true);
        for (int i = 0; i < m->size(); ++i)
        {
            auto& cl = m->element(i);
            vec<Lit> lits(cl.size());
            int index = 0;
            for(int id : cl)
            {
                assert(id != 0);
                int var = abs(id) - 1;
                Lit l = (id > 0) ? mkLit(var) : ~mkLit(var);
                lits[index++] = l;
            }
            bool res = sslv->addClause(lits);
        }
        sslv->eliminate(true);

        while(nVars() < sslv->nVars())
            newVar();

        for (Glucose::ClauseIterator c = sslv->clausesBegin(); c != sslv->clausesEnd(); ++c)
        {
            const Glucose::Clause &cls = *c;
            Glucose::vec<Glucose::Lit> cls_;
            for (int i = 0; i < cls.size(); ++i)
            {
                cls_.push(cls[i]);
            }
            addClause(cls_);
        }

        for (auto c = sslv->trailBegin(); c != sslv->trailEnd(); ++c)
        {
            addClause(*c);
        }
        delete sslv;
    }

    bool MainSolver::badcheck(const Assignment &st, const int bad)
    {
        // assumption = {st, bad}
        assumptions.clear();
        assumptions.push (SAT_lit (bad));
        push_to_assumption(st);
        return solve_assumption();
    }

    void MainSolver::push_flags_M(const int level){
        assert(level > -1);
        // activate this frame
        auto flagOfThisFrame = MFlagOf(level);
        assumptions.push (SAT_lit (flagOfThisFrame));
        // deactivate others
        for(auto flg : MFlags)
        {
            if (flg != flagOfThisFrame)
                assumptions.push(SAT_lit(-flg));
        }
    }

    void MainSolver::push_to_assumption(const Cube &cu)
    {
        for(auto& lit :cu)
            assumptions.push(SAT_lit(lit));
    }

    void MainSolver::push_to_assumption_primed(const Cube &cu)
    {
        for(auto& lit :cu)
        {
            if(_model->latch_var(lit))
                assumptions.push(SAT_lit(_model->prime(lit)));
            else
                assumptions.push(SAT_lit(lit));
        }
    }

    void MainSolver::set_assumption_primed(Cube& uc_or_flags) {
        assumptions.clear();
        push_to_assumption_primed(uc_or_flags);
    }

    /**
	 * @brief set assumption = { s->getLatches() , MFlagOf(Os[frame_level]) }
	 * 
	 * @param s 
	 * @param frame_level 
	 */
    void MainSolver::set_assumption_M(State*s,  const int frame_level)
    {
        // clear the assumptions first
        assumptions.clear();
        
        // set flags, both activation flags and deactivation flags
        push_flags_M(frame_level);
        
        // push the latches of this state `s`.

        if(reverseT)
            push_to_assumption_primed(s->getLatches());
        else
            push_to_assumption(s->getLatches());
    }

	void MainSolver::set_assumption_M(State*s,  const int frame_level, const std::vector<Cube> &prefers)
    {
        // clear the assumptions first
        assumptions.clear();

        // if level == -1: badcheck.
        if(frame_level == -1)
        {
            assert(!reverseT);
            push_to_assumption({bad});
        }
        else
            // set flags, both activation flags and deactivation flags
            push_flags_M(frame_level);

        for(auto& cu:prefers)
        {
            if(reverseT)
                push_to_assumption_primed(cu);
            else
                push_to_assumption(cu);
        }

        {
            if(reverseT)
                push_to_assumption_primed(s->getLatches());
            else
                push_to_assumption(s->getLatches());
        }
    }

    void MainSolver::set_expectation(const std::vector<int>& expectations)
    {
        assert(!reverseT && "temporarily not avaiable for forward-car");
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
            if(reverseT)
                setPolarity(var, sgn);
            else
                setPolarity(id_prime, sgn);
        }
        
    }


	//NOTE: this State* is not owned by this solver, but the checker. It should be immediately added into clear duty.	
	State* MainSolver::get_state()
	{
		assert(unroll_level == 1);
		Assignment model = get_model();
		if(!reverseT)
			shrink_model (model);
		Assignment inputs(model.begin(),model.begin() + _model->num_inputs());
		Assignment latches(model.begin() + _model->num_inputs(),model.begin() + _model->num_inputs()+_model->num_latches());
		State* s = new State(inputs,latches);
		return s;
	}

    State* MainSolver::getInitState()
    {
        assert(unroll_level == 1);
		Assignment model = get_model();
        Assignment initLatches(model.begin() + _model->num_inputs(),model.begin() + _model->num_inputs()+_model->num_latches());
        State *newInit = new State(initLatches);
        return newInit;
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
		return res;
	}
	
	/**
	 * @brief get unsat core, get latch, shrink to previous if reverseT. 
	 * @note the last lit is pushed again.
	 * @return Cube 
	 */
	Cube MainSolver::get_shrunk_uc ()
	{
		Cube conflict = get_uc ();

		assert(!conflict.empty());
		if (reverseT)
		{
		    _model->shrink_to_previous_vars (conflict);
		}
		else
		{
		    _model->shrink_to_latch_vars (conflict);
		}

        if(!uc_no_sort)
        {
            std::sort (conflict.begin (), conflict.end (), car::absIncr);
        }
		return conflict;
	}

    Cube MainSolver::get_conflict_another (int option, int nth)
	{
		Cube conflict = get_uc_another (option, nth);
		assert(!conflict.empty());
		if (reverseT)
		{
		    _model->shrink_to_previous_vars (conflict);
		}
		else
		{
		    _model->shrink_to_latch_vars (conflict);
		}

        if(!uc_no_sort)
        {
            std::sort(conflict.begin(), conflict.end(), car::absIncr);
        }
        return conflict;
	}
	
	void MainSolver::add_new_frame_M(const OFrame& frame, const int frame_level)
	{
		for(size_t i = 0; i < frame.size (); i ++)
		{
			add_clause_from_cube_M (frame[i], frame_level);
		}
	}

	void MainSolver::add_new_frame_P(const OFrame& frame, const int frame_level)
	{
		for(size_t i = 0; i < frame.size (); i ++)
		{
			add_clause_from_cube_P (frame[i], frame_level);
		}
	}

	// Actually, each center_state should only exist in O sequence in one direction.
	// Should we just record this?
	void MainSolver::add_clause_from_cube_M(const Cube &cu, const int frame_level)
	{
		int flag = MFlagOf(frame_level);
		vector<int> cl = {-flag};
		for(size_t i = 0; i < cu.size (); i ++)
		{
			if (!reverseT)
				cl.push_back (-_model->prime (cu[i]));
			else
				cl.push_back (-cu[i]);
		}
		add_clause (cl);
	}

	void MainSolver::add_clause_from_cube_P(const Cube &cu, const int frame_level)
	{
		assert(!reverseT && "backward CAR only for propagation\n");
		int flag = PFlagOf(frame_level);
		vector<int> cl = {-flag};
		for(size_t i = 0; i < cu.size (); i ++)
		{
			cl.push_back (-cu[i]);
		}
		add_clause (cl);
	}

    // for propagation use only.
    void MainSolver::add_clause_from_cube(const Cube &cu, int flag, bool __placeholder)
    {
        vector<int> cl;
        cl.reserve(cu.size() +1);
        // flag
        if(flag >= 0)
            cl.push_back(flag);
		// negate of cu.
        for(size_t i = 0; i < cu.size (); i ++)
		{
            cl.push_back (-cu[i]);
		}
		add_clause (cl);
    }

    void MainSolver::shrink_model(Assignment& model)
    {
		Assignment res = model;

		int num_inputs = _model->num_inputs();
		int num_latches = _model->num_latches();
		
		for(int i = num_inputs + 1; i <= num_inputs + num_latches; i++)
		{
			int p = _model->prime(i);
			assert(p != 0);
			
			int abs_p = abs(p);
			assert(model.size() >= size_t(abs_p));
			
			int val = model[abs_p - 1];
			
			if (p == val)
				res[i - 1] = i;
			else
				res[i - 1] = -i;
		}
		model = res;
	}

}
