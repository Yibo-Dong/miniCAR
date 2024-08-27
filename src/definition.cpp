#include "definition.h"

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <map> 
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stack>
#include "definition.h"

using namespace std;

namespace car {

    // the special negp state
    // TODO: try unfold negp instead of this.
    const State* State::negp_state;

	int State::_maxid = 0;
    int State::_num_inputs = 0;
 	int State::_num_latches = 0;
 	
    /**
     * @brief Whether this state is already blocked by this cube.
     * @pre need state(this) to follow the literal order.
     * @param cu : The UC 
     * @return true : is blocked
     * @return false : not blocked
     */
    bool State::imply(const Cube &uc) const
    {
        if (isPartial())
        {
            std::unordered_set<int> lits;
            for (auto &lit : getLatches())
            {
                lits.insert(lit);
            }
            for (auto &l : uc)
            {
                if (lits.find(l) == lits.end())
                    return false;
            }
            return true;
        }
        else
        {
            for (int i = uc.size() - 1; i >= 0; --i)
            {
                // get the offset of this literal.
                int index = abs(uc[i]) - _num_inputs - 1;
                assert(index >= 0);
                if (_latches[index] ^ uc[i])
                {
                    return false;
                }
            }
            return true;
        }
    }

    /**
     * @brief calculate intersection of this->latches() and cu.
     * @pre need state(this) to follow the literal order.
     * @post follow the order in cu.
     * @param cu 
     * @return true 
     * @return false 
     */
    Cube State::intersect(const Cube &cu) const
    {
        Cube res;
        res.reserve(cu.size());
        if (isPartial())
        {
            std::unordered_set<int> lits;
            for (auto &lit : getLatches())
            {
                lits.insert(lit);
            }
            for (auto &l : cu)
            {
                if (lits.find(l) == lits.end())
                    res.push_back(l);
            }
        }
        else
        {
            for (size_t i = 0; i < cu.size(); i++)
            {
                int index = abs(cu[i]) - _num_inputs - 1;
                assert(index >= 0);
                if (_latches[index] == cu[i])
                    res.push_back(cu[i]);
            }
        }
        return res;
    }

    /**
     * @brief get input values as a string. Used in CEX printing.
     * 
     * @return string 
     */
	string State::getInputsStr() const
	{
		string res = "";
		for(size_t i = 0; i < _inputs.size(); i++)
			res += (_inputs.at(i) > 0) ? "1" : "0";
		return res;
	}

    /**
     * @brief get latch values as a string. Used in CEX printing.
     * 
     * @return string 
     */
	string State::getLatchesStr() const
	{
		string res = "";
		size_t j = 0;
		for(int i = 0; i < _num_latches; i++)
		{
			if (j == _latches.size())
				res += "x";
			else if (_num_inputs + i + 1 < abs(_latches.at(j)))
				res += "x";
			else
			{
				res += (_latches.at(j) > 0) ? "1" : "0";
				j++;
			}
		}
		return res;
	}

    Problem::Problem (aiger* aig, bool unroll): aig(aig)
	{
		// According to aiger format, inputs should be [1 ... _num_inputs]
		// and latches should be [num_inputs+1 ... num_latches+num_inputs]]
		_num_inputs         = aig->num_inputs;
		_num_latches        = aig->num_latches;
		_num_ands           = aig->num_ands;
		_num_constraints    = aig->num_constraints;
		_num_outputs        = aig->num_outputs;

		// preserve two more ids for TRUE (_max_id - 1) and FALSE (_max_id)
        _max_id = aig->maxvar + 2;
        _true = _max_id - 1;
		_false = _max_id;
		
		collectTrue (aig);
		
		set_constraints (aig);
		set_outputs (aig);
		set_init (aig);
		
		create_next_map (aig);
        if(unroll)
            unroll_prime(aig);
		create_clauses (aig);
	}
	
	// collect those that are trivially constant
	void Problem::collectTrue (const aiger* aig)
	{
		for(size_t i = 0; i < aig->num_ands; i ++)
		{
			aiger_and& aa = aig->ands[i];
			// lhs of an and gate is always even in aiger
			assert (aa.lhs % 2 == 0);
            // ? = T && T
			if (is_true (aa.rhs0) && is_true (aa.rhs1))
				_true_set.insert (aa.lhs);
            // ? = F && ? 
            // or
            // ? = ? && F
			else if (is_false (aa.rhs0) || is_false (aa.rhs1))
				_true_set.insert (aa.lhs + 1);
			// lhs = 8 && 7 = a4 && ~a3 is not constant
			// lhs = 9 && 8 = ~a4 && a4 is constant
			else if (aa.rhs0 == aa.rhs1+1 && aa.rhs0 % 2 == 1)
				_true_set.insert (aa.lhs + 1);
		}
	}

	void Problem::set_constraints (const aiger* aig)
	{
		for(size_t i = 0; i < aig->num_constraints; i++)
		{
			int lit = (int)aig->constraints[i].lit;
			_constraints.push_back(car_var(lit));
		}
	}
	
	void Problem::set_outputs (const aiger* aig)
	{
		for(size_t i = 0; i < aig->num_outputs; ++i)
		{
			int lit = (int)aig->outputs[i].lit;
			_outputs.push_back(car_var(lit));
		}
        for(size_t i = 0; i< aig->num_bad; ++i)
        {
            int lit = (int)aig->bad[i].lit;
			_outputs.push_back(car_var(lit));
        }
        _num_outputs = _outputs.size();
        assert(_outputs.size() == 1 && "Only 1 output is allowed!\n");
	}

	void Problem::set_init (const aiger* aig)
	{
		for(size_t i = 0; i < aig->num_latches; ++i)
		{
			if (aig->latches[i].reset == 0)
				_init.push_back (-(_num_inputs+1+i));
			else if (aig->latches[i].reset == 1)
				_init.push_back (_num_inputs+1+i);
			else
			{
                // do nothing. They are `don't care`.
			}
		}
	}
	
	void Problem::create_next_map (const aiger* aig)
	{
		for(size_t i = 0; i < aig->num_latches; i ++)
		{
			int val = (int)aig->latches[i].lit;
			//a latch should not be a negative number
			assert (val % 2 == 0);
			val = val / 2;
			//make sure our assumption about latches is correct
			assert (val == (_num_inputs + 1 + i));
			
			//pay attention to the special case when next_val = 0 or 1
			if (is_false (aig->latches[i].next))  //FALSE
			{
				_next_map.insert (std::pair<int, int> (val, _false));
				insert_to_reverse_next_map (_false, val);
			}
			else if (is_true (aig->latches[i].next)) //TRUE
			{
				_next_map.insert (std::pair<int, int> (val, _true));
				insert_to_reverse_next_map (_true, val);
			}
			else
			{
				int next_val = (int) aig->latches[i].next;
				next_val = (next_val % 2 == 0) ? (next_val/2) : -(next_val/2);
				_next_map.insert (std::pair<int, int> (val, next_val));
				insert_to_reverse_next_map (abs (next_val), (next_val > 0) ? val : -val);
			}
		}
	}
	
	void Problem::insert_to_reverse_next_map (const int index, const int val)
	{
	    auto it = _reverse_next_map.find (index);
	    if (it == _reverse_next_map.end ())
	    {
	        vector<int> v;
	        v.push_back (val);
	        _reverse_next_map.insert (std::pair<int, vector<int> > (index, v));
	    }
	    else
	        (it->second).push_back (val);
	}
	
	void Problem::create_clauses (const aiger* aig)
	{
		// contraints, outputs and latches gates are stored in order,
		// as the need for start solver construction
		std::set<unsigned> exist_gates, gates;

		// ============================================================================
		// (1) create clauses for constraints encoding
		// ============================================================================
		collect_necessary_gates(aig, aig->constraints, aig->num_constraints, exist_gates, gates);

		for (std::set<unsigned>::iterator it = gates.begin(); it != gates.end(); ++it)
		{
			aiger_and *aa = aiger_is_and(const_cast<aiger *>(aig), *it);
			assert(aa != NULL);
			add_clauses_from_equation(aa);
		}

        // ============================================================================
		// (1+) create clauses for constraints encoding
		// ============================================================================
        for(size_t i = 0; i < aig->num_constraints; ++i)
        {
            auto lit = aig->constraints[i].lit;
            auto var = car_var(lit);
            if(lit == 1)
            {
                // true, do nothing
            }
            else if(lit == 0)
            {
                cls_.push_back({false_id()});
            }
            else // not a constant
                cls_.push_back({var});
        }
		// ============================================================================
		// (2) same next have same previous, initialized to 0.
		// ============================================================================
		// if several latches' next value points to a same one, the values of these latches also have the constraint
		// But NOTE that all states in the aiger model must satisfy this constraint, EXCEPT the initial state
		// Also it is needed to add new contraints later
		
		set_common_next_start();
		// @note: even if no gates share a common next, the helper lit will be inserted.
		create_constraints_for_latches();
		// ============================================================================
		// (3) clause for encoding outputs.
		// ============================================================================
		// create clauses for outputs

		set_outputs_start();

		gates.clear();
		//  Use Outputs as the start point, recursively add all the rhs lits.
        
        // Assume: only one property to be checked.
        if(aig->num_outputs > 0)
    		collect_necessary_gates(aig, aig->outputs, aig->num_outputs, exist_gates, gates);
        else
            collect_necessary_gates(aig, aig->bad, aig->num_bad, exist_gates, gates);
    	for (std::set<unsigned>::iterator it = gates.begin(); it != gates.end(); ++it)
		{
			aiger_and *aa = aiger_is_and(const_cast<aiger *>(aig), *it);
			assert(aa != NULL);
			add_clauses_from_equation(aa);
		}
		// ============================================================================
		// (4) clause for encoding latches's mapping relation
		// ============================================================================
		// create clauses for latches

		set_latches_start();

		//  Use Next value of Latches as the start point, recursively add all the rhs lits.
		gates.clear();
		collect_necessary_gates(aig, aig->latches, aig->num_latches, exist_gates, gates, true);
		for (std::set<unsigned>::iterator it = gates.begin(); it != gates.end(); ++it)
		{
			aiger_and *aa = aiger_is_and(const_cast<aiger *>(aig), *it);
			assert(aa != NULL);
			add_clauses_from_equation(aa);
		}
		// ============================================================================
		//// (5) create clauses for true and false
		// ============================================================================
		cls_.push_back(clause(_true));
		cls_.push_back(clause(-_false));
	}

    // when dealing with cosntraints, we need to add primed constraints to solver.
    void Problem::unroll_prime(const aiger* aig)
    {
        // for each andgate, we would like to get its primed version.

        for(int i = 0; i < aig->num_ands; ++i)
		{
			aiger_and& aa = aig->ands[i];
            int primedLHS = ++_max_id;
            auto toInsert = std::pair<int,int>{car_var(aa.lhs), primedLHS};
            _next_map.insert(toInsert);
            insert_to_reverse_next_map(primedLHS, car_var(aa.lhs));
            int primedRHS0 = (input_var(car_var(aa.rhs0))) ? car_var(aa.rhs0) : prime(car_var(aa.rhs0));
            int primedRHS1 = (input_var(car_var(aa.rhs1))) ? car_var(aa.rhs1) :prime(car_var(aa.rhs1));
            cls_.push_back (clause ( (primedLHS), -(primedRHS0), -(primedRHS1)));
			cls_.push_back (clause (- (primedLHS),  (primedRHS0)));
			cls_.push_back (clause (- (primedLHS),  (primedRHS1)));
        }
    }

    /**
	 * @brief 
	 * flag1: All latches that share a common next have same value 
	 * flag2: All latches are initialized to 0
	 * flags are enabled when assigned to 1.
	 */
	void Problem::create_constraints_for_latches()
	{
		int flag1 = ++_max_id;

		bool exist = false;
		for (auto it = _reverse_next_map.begin(); it != _reverse_next_map.end(); ++it)
		{
			vector<int> &v = it->second;
			if (v.size() <= 1)
				continue;
			/**
			 * @brief As to latches { (l1,n), (l2,n), (l3, n) }
			 * we assume val(l1) == val(l2) == val(l3)
			 * which means, l1 <-> l2 <-> l3
			 * <=> (l1, ~l2), (~l1, l2), (l2, ~l3), (~l2, l3)
			 */
			// TODO: Try to record this equivalence relationship using a map and simplify the model.
			exist = true;
			for(size_t i = 0; i < v.size() - 1; i++)
			{
				cls_.push_back({v[i], -v[i + 1], -flag1});
				cls_.push_back({-v[i], v[i + 1], -flag1});
			}
		}
		if (!exist)
		{
			// at last add one clause for flag1 so that the sat solver will not treat flag1 as one literal
			cls_.push_back({++_max_id, -flag1});
		}

		// add initial state
		int flag2 = ++_max_id;
		for(size_t i = 0; i < _init.size(); i++)
		{
			cls_.push_back({_init[i], -flag2});
		}
		// either : all latches share a common next
		// or : this is the initial states, all latches are initialized.
		cls_.push_back({flag1, flag2});
	}

	void Problem::collect_necessary_gates (const aiger* aig, const aiger_symbol* as, const int as_size, 
	                                        std::set<unsigned>& exist_gates, std::set<unsigned>& gates, bool next)
	{
		for(int i = 0; i < as_size; i ++)
		{
			aiger_and* aa;
			if (next) 
			    aa = necessary_gate (as[i].next, aig);
			else
			{
			    aa = necessary_gate (as[i].lit, aig);
			    if (aa == NULL)
			    {
			    	if (is_true (as[i].lit))
			    		_outputs[i] = _true;
			    	else if (is_false (as[i].lit))
			    		_outputs[i] = _false;
			    }
			}
			iteratively_add (aa, aig, exist_gates, gates);	
		}
		
	}
	
	aiger_and* Problem::necessary_gate (const unsigned id, const aiger* aig)
	{
		if (!is_true (id) && !is_false (id))
			return aiger_is_and (const_cast<aiger*> (aig), (id % 2 == 0) ? id : (id-1));
			
		return NULL;
	}
	
	// deprecated. Recursion causes stack overflow on large cases.
	void Problem::recursively_add (const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates)
	{
		if (aa == NULL)
			return;
		if (exist_gates.find (aa->lhs) != exist_gates.end ())
			return;
		
		gates.insert (aa->lhs);
		exist_gates.insert (aa->lhs);
		aiger_and* aa0 = necessary_gate (aa->rhs0, aig);
		recursively_add (aa0, aig, exist_gates, gates);
		
		aiger_and* aa1 = necessary_gate (aa->rhs1, aig);
		recursively_add (aa1, aig, exist_gates, gates);
	}

	void Problem::iteratively_add(const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates)
	{
		if(aa == NULL)
			return;
		std::stack<const aiger_and*> stk;
		stk.push(aa);
		while (!stk.empty())
		{
			const aiger_and* current = stk.top();
			stk.pop();

			if (exist_gates.find(current->lhs) != exist_gates.end())
				continue;

			gates.insert(current->lhs);
			exist_gates.insert(current->lhs);

			aiger_and* aa0 = necessary_gate(current->rhs0, aig);
			if (aa0 != NULL)
				stk.push(aa0);

			aiger_and* aa1 = necessary_gate(current->rhs1, aig);
			if (aa1 != NULL)
				stk.push(aa1);
		}
	}
	
	/**
	 * @brief 
	 * 
	 * @param aa 
	 */
	void Problem::add_clauses_from_equation (const aiger_and* aa)
	{
		assert (aa != NULL);
		assert (!is_true (aa->lhs) && !is_false (aa->lhs));
		
		if (is_true (aa->rhs0))
		{
			/**
			 * @brief lhs = True /\ rhs1
			 * <=>	  lhs <-> rhs1
			 * <=>	  lhs -> rhs1 /\ rhs1 -> lhs
			 * <=>	  (~lhs \/ rhs1), (~rhs1, lhs)
			 */
			cls_.push_back (clause (car_var (aa->lhs), -car_var (aa->rhs1)));
			cls_.push_back (clause (-car_var (aa->lhs), car_var (aa->rhs1)));
		}
		else if (is_true (aa->rhs1))
		{
			cls_.push_back (clause (car_var (aa->lhs), -car_var (aa->rhs0)));
			cls_.push_back (clause (-car_var (aa->lhs), car_var (aa->rhs0)));
		}
		else
		{
			cls_.push_back (clause (car_var (aa->lhs), -car_var (aa->rhs0), -car_var (aa->rhs1)));
			cls_.push_back (clause (-car_var (aa->lhs), car_var (aa->rhs0)));
			cls_.push_back (clause (-car_var (aa->lhs), car_var (aa->rhs1)));
		}
			
	}
	
	/**
	 * @brief get the next id of given id
	 * 
	 * @param id 
	 * @return int 
	 */
	int Problem::prime (const int id) const
	{
		auto it = _next_map.find (abs (id));
		if (it == _next_map.end ())
		{
			cerr<<"[Fatal Error] Cannot find its prime:" << id<< endl;
			assert(false && "cannot find its prime");
		    return 0; //not found
		}
		return (id > 0 ? it->second : -(it->second));
	}
	
	/**
	 * @brief get the previous ids of given id
	 * @note there may be more than one.
	 * @param id 
	 * @return std::vector<int> 
	 */
	std::vector<int> Problem::previous (const int id) const
	{
	    vector<int> res;
	    auto it = _reverse_next_map.find (abs (id));
		if (it == _reverse_next_map.end ())
		    return res; //not found
		res = it->second;
		if (id < 0)
		{
		    for(size_t i = 0; i < res.size (); i ++)
		        res[i] = -res[i];
		}
		return res;
	}
	
	/**
	 * @brief collect the previous vars from uc.  
	 * @note multiple latch may share a common next, therefore one lit may have several previous.
	 * @param uc 
     * TODO: check it.
	 */
	void Problem::shrink_to_previous_vars (Cube& uc) const
	{
		Cube tmp;
		for(size_t i = 0; i < uc.size (); i ++)
		{
		    vector<int> ids = previous (abs (uc[i]));
			if (ids.empty ())
			    continue;
			// 每一个pre都放进去。
			for(size_t j = 0; j < ids.size (); j ++)
				tmp.push_back ((uc[i] > 0) ? ids[j] : (-ids[j]));	
		}
		uc = tmp;
	}

    void Problem::transform_uc_to_previous (Cube& uc) const
	{
    	for(size_t i = 0; i < uc.size(); ++i)
		{
		    vector<int> ids = previous (abs (uc[i]));
			if (ids.empty ())
            {
                assert(false && "should not happen");
            }
			else
                // take the first one.
				uc[i] = ((uc[i] > 0) ? ids[0] : (-ids[0]));	
		}
	}
	
	/**
	 * @brief collect the latch vars in uc.  
	 * @pre since we now set all the flags t obe at the beginning of the assumption, they shall be at the end.
	 * @param uc 
	 */
	void Problem::shrink_to_latch_vars (Cube& uc) const
	{
        if(!latch_var(abs(uc.back())))
            uc.pop_back();
        // for(auto i :uc)
        // {
        //     assert(latch_var(abs(i)));
        // }
	}


}
