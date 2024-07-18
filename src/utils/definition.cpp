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
#include "definition.h"

using namespace std;

namespace car {


const State* State::negp_state;
 	
    /**
     * @brief Whether this state is already blocked by this cube.
     * @pre need state(this) to follow the literal order.
     * @param cu 
     * @return true 
     * @return false 
     */
 	bool State::imply (const Cube& cu) const
	{
		for (int i = cu.size() - 1 ; i >= 0; --i)
		{
			int index = abs(cu[i]) - num_inputs_ - 1;
			assert (index >= 0);
			if (s_[index] ^ cu[i])
			{
				return false;
			}
		}
		return true;
	}
	
    /**
     * @brief calculate intersection of two cubes.
     * @pre need state(this) to follow the literal order.
     * @post follow the order in cu.
     * @param cu 
     * @return true 
     * @return false 
     */
	Cube State::intersect (const Cube& cu) 
	{
		Cube res;
		for (int i = 0; i < cu.size (); i ++)
		{
			int index = abs(cu[i]) - num_inputs_ - 1;
			assert (index >= 0);
			if (s_[index] == cu[i])
				res.push_back (cu[i]);
		}
		return res;
	}

	Cube State::intersect_no_dupl (const Cube& cu) 
	{
		Cube res;
		std::unordered_set<int> mark;
		for (int i = 0; i < cu.size (); i ++)
		{
			int index = abs(cu[i]) - num_inputs_ - 1;
			assert (index >= 0);
			if (s_[index] == cu[i] && mark.find(cu[i]) == mark.end())
			{
				res.push_back (cu[i]);
				mark.insert(cu[i]);
			}
		}
		mark.clear();
		return res;
	}

	Cube State::intersect (const State *s) 
	{
		Cube res;
		Cube latches = s->s();

		for (int i = 0; i < latches.size (); i ++)
		{
			if (s_[i] == latches[i])
				res.push_back( s_[i] > 0 ? (i +  num_inputs_ + 1) : -(i +  num_inputs_ + 1));
		}
		return res;
	}
 	

	string State::inputs() const
	{
		string res = "";
		for (int i = 0; i < inputs_.size(); i++)
			res += (inputs_.at(i) > 0) ? "1" : "0";
		return res;
	}

	string State::latches() const
	{
		string res = "";
		int j = 0;
		for (int i = 0; i < num_latches_; i++)
		{
			if (j == s_.size())
				res += "x";
			else if (num_inputs_ + i + 1 < abs(s_.at(j)))
				res += "x";
			else
			{
				res += (s_.at(j) > 0) ? "1" : "0";
				j++;
			}
		}
		return res;
	}

	int State::num_inputs_ = 0;
 	int State::num_latches_ = 0;

	unsigned State::next_id_ = 0;

	Cube State::next_latches()
	{
		// if it's negp, nothing is already known.
		if(is_negp)
			return {};

		Cube known;
		std::unordered_map<int,int> value;
		// value for 0/1 is constant.
		value[0] = 0;
		
		// value for inputs.
		// if(inputs_.empty())
		// {
			for (int i = 0; i < aig_->num_inputs; ++i)
			{
				value[i+1] = -1;
			}
		// }else{
		// 	for(int i = 0; i < aig_->num_inputs;++i)
		// 	{
		// 		value[i+1] = (inputs_[i] > 0) ? 1 : 0;
		// 	}
		// }

		// value for latches
		for(int i = 0; i < s_.size();++i)
		{
			int offset = abs(s_[i] -1 - int(aig_->num_inputs));
			value[1+aig_->num_inputs +i ] = (s_[i]>0) ? 1 : 0;
		}
		
		// auto getLitValue = [&value](int x){
		// 	assert(value.find(x) != value.end() || value.find(-x) != value.end());
		// 	if(x<0){return value[x] == -1 ? -1 : 1 - value[x];}
		// 	else{return value[x];}
		// };

		auto getAigValue = [&value](int x){
			// assert(value.find(x>>1) != value.end() || value.find(1 + (x>>1)) != value.end());
			if(x & 1){ assert(value.find(x>>1) != value.end());  return value[x>>1] == -1 ? -1 : 1 - value[x>>1];}
			else{assert(value.find(x>>1) != value.end()); return value[x>>1];}
		};
		
		// value for and gates
		for (int i = 0; i < aig_->num_ands; i ++)
		{
			aiger_and& aa = aig_->ands[i];
			int rhs0_value = getAigValue(aa.rhs0);
			int rhs1_value = getAigValue(aa.rhs1);
			int lhs_value;

			if(rhs0_value == 0 || rhs1_value == 0)
				lhs_value = 0;
			// taint will spread. mark as not determined
			else if(rhs0_value == -1 || rhs1_value ==-1)
			{
				lhs_value = -1;
				if(aa.rhs0 == aa.rhs1)
					lhs_value = 1;
				else if(aa.rhs0 == aa.rhs1 + 1 && aa.rhs0&1 )
					lhs_value = 0;
				// if(lhs_value != -1)
				// 	cerr<<"?= "<<aa.rhs0<<endl;
			}
			else if (rhs0_value == 1 && rhs1_value == 1)
				lhs_value = 1;
			// else
			// {
			// 	cerr<<"strange value,"<<aa.rhs0 <<":" <<rhs0_value<<", "<<aa.rhs1<<": "<<rhs1_value<<endl;
			// }
			value[aa.lhs >> 1] = lhs_value;

			// cerr<<aa.lhs<<" = "<<aa.rhs0 <<" + "<<aa.rhs1<<endl;
			// cerr<<getAigValue(aa.lhs)<<" = "<<getAigValue(aa.rhs0)<< " + "<<getAigValue(aa.rhs1) <<endl;
		}

		// calculate necessaries.
		for(int i = 0; i< aig_->num_latches; ++i)
		{
			int offset = 1+aig_->num_inputs;
			int next_value = getAigValue(aig_->latches[i].next);
			if(next_value!=-1)
			{
				// cerr<<(next_value>0);
				known.push_back(next_value==0? -(i + offset) : i+ offset);
			}
			// else{
			// 	cerr<<"x";
			// }
		}
		// if(id == 5)
		// {
		// 	std::vector<int> values(aig_->maxvar,0);
			
		// 	cout<<"values: ";
		// 	for(auto &p:value)
		// 	{
		// 		values[p.first] = p.second;
		// 	}	
		// 	for(int i = 0; i< values.size();++i)
		// 		cout<<i<<" : "<<values[i]<<", ";
		// 	cout<<endl;
		// 	cout<<"known :";
		// 	for(int i:known)
		// 	cout<<i<<", ";
		// 	cout<<endl;
		// }
		// cout<<"known:";
		// for(int i:known)
		// 	cout<<i<<", ";
		// cout<<endl;


		
		return known;
	}

	std::map<Assignment, int> State::ids;

Model::Model (aiger* aig, const bool verbose)
	{
		verbose_ = verbose;
		// According to aiger format, inputs should be [1 ... num_inputs_]
		// and latches should be [num_inputs+1 ... num_latches+num_inputs]]
		num_inputs_ = aig->num_inputs;
		num_latches_ = aig->num_latches;
		num_ands_ = aig->num_ands;
		num_constraints_ = aig->num_constraints;
		num_outputs_ = aig->num_outputs;

		//preserve two more ids for TRUE (max_id_ - 1) and FALSE (max_id_)
		max_id_ = aig->maxvar+2;
		max_model_id_ =  aig->maxvar+2;
		true_ = max_id_ - 1;
		false_ = max_id_;
		
		collect_trues (aig);
		
		set_constraints (aig);
		set_outputs (aig);
		set_init (aig);
		
		create_next_map (aig);
		create_clauses (aig);
	}
	
	// collect those that is trivially constant
	void Model::collect_trues (const aiger* aig)
	{
		for (int i = 0; i < aig->num_ands; i ++)
		{
			aiger_and& aa = aig->ands[i];
			// lhs of an and gate is always even in aiger
			assert (aa.lhs % 2 == 0);
			if (is_true (aa.rhs0) && is_true (aa.rhs1))
				trues_.insert (aa.lhs);
			else if (is_false (aa.rhs0) || is_false (aa.rhs1))
				trues_.insert (aa.lhs + 1);
			// lhs = 8 * 7 = a4 * ~a3 is not constant
			// lhs = 9 * 8 = ~a4 * a4 is constant
			else if (aa.rhs0 == aa.rhs1+1 && aa.rhs0 % 2 == 1)
				trues_.insert (aa.lhs + 1);
		}
	}

	void Model::set_constraints (const aiger* aig)
	{
		for (int i = 0; i < aig->num_constraints; i++)
		{
			int lit = (int)aig->constraints[i].lit;
			constraints_.push_back(car_var(lit));
		}
	}
	
	void Model::set_outputs (const aiger* aig)
	{
		for (int i = 0; i < aig->num_outputs; i++)
		{
			int lit = (int)aig->outputs[i].lit;
			outputs_.push_back(car_var(lit));
		}
	}

	void Model::set_init (const aiger* aig)
	{
		for (int i = 0; i < aig->num_latches; i ++)
		{
			if (aig->latches[i].reset == 0)
				init_.push_back (-(num_inputs_+1+i));
			else if (aig->latches[i].reset == 1)
				init_.push_back (num_inputs_+1+i);
			else
			{
				cout << "Error setting initial state!" << endl;
				exit (0);
			}
		}
	}
	
	void Model::create_next_map (const aiger* aig)
	{
		for (int i = 0; i < aig->num_latches; i ++)
		{
			int val = (int)aig->latches[i].lit;
			//a latch should not be a negative number
			assert (val % 2 == 0);
			val = val / 2;
			//make sure our assumption about latches is correct
			assert (val == (num_inputs_ + 1 + i));
			
			//pay attention to the special case when next_val = 0 or 1
			if (is_false (aig->latches[i].next))  //FALSE
			{
				next_map_.insert (std::pair<int, int> (val, false_));
				insert_to_reverse_next_map (false_, val);
			}
			else if (is_true (aig->latches[i].next)) //TRUE
			{
				next_map_.insert (std::pair<int, int> (val, true_));
				insert_to_reverse_next_map (true_, val);
			}
			else
			{
				int next_val = (int) aig->latches[i].next;
				next_val = (next_val % 2 == 0) ? (next_val/2) : -(next_val/2);
				next_map_.insert (std::pair<int, int> (val, next_val));
				insert_to_reverse_next_map (abs (next_val), (next_val > 0) ? val : -val);
			}
			
		}
	}
	
	void Model::insert_to_reverse_next_map (const int index, const int val)
	{
	    reverseNextMap::iterator it = reverse_next_map_.find (index);
	    if (it == reverse_next_map_.end ())
	    {
	        vector<int> v;
	        v.push_back (val);
	        reverse_next_map_.insert (std::pair<int, vector<int> > (index, v));
	    }
	    else
	        (it->second).push_back (val);
	}
	
	void Model::create_clauses (const aiger* aig)
	{
		// contraints, outputs and latches gates are stored in order,
		// as the need for start solver construction
		std::set<unsigned> exist_gates, gates;


		// ============================================================================
		// (1) create clauses for constraints encoding
		// ============================================================================
		collect_necessary_gates(aig, aig->constraints, aig->num_constraints, exist_gates, gates);

		for (std::set<unsigned>::iterator it = gates.begin(); it != gates.end(); it++)
		{
			aiger_and *aa = aiger_is_and(const_cast<aiger *>(aig), *it);
			assert(aa != NULL);
			add_clauses_from_equation(aa);
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
		collect_necessary_gates(aig, aig->outputs, aig->num_outputs, exist_gates, gates);
		for (std::set<unsigned>::iterator it = gates.begin(); it != gates.end(); it++)
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
		for (std::set<unsigned>::iterator it = gates.begin(); it != gates.end(); it++)
		{
			aiger_and *aa = aiger_is_and(const_cast<aiger *>(aig), *it);
			assert(aa != NULL);
			add_clauses_from_equation(aa);
		}

		// ============================================================================
		//// (5) create clauses for true and false
		// ============================================================================
		cls_.push_back(clause(true_));
		cls_.push_back(clause(-false_));
	}

	/**
	 * @brief 
	 * flag1: All latches that share a common next have same value 
	 * flag2: All latches are initialized to 0
	 * flags are enabled when assigned to 1.
	 */
	void Model::create_constraints_for_latches()
	{
		int flag1 = ++max_id_;

		bool exist = false;
		for (reverseNextMap::iterator it = reverse_next_map_.begin(); it != reverse_next_map_.end(); it++)
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
			for (int i = 0; i < v.size() - 1; i++)
			{
				cls_.push_back({v[i], -v[i + 1], -flag1});
				cls_.push_back({-v[i], v[i + 1], -flag1});
			}
		}
		if (!exist)
		{
			// at last add one clause for flag1 so that the sat solver will not treat flag1 as one literal
			cls_.push_back({++max_id_, -flag1});
		}

		// add initial state
		int flag2 = ++max_id_;
		for (int i = 0; i < init_.size(); i++)
		{
			cls_.push_back({init_[i], -flag2});
		}
		// either : all latches share a common next
		// or : this is the initial states, all latches are initialized.
		cls_.push_back({flag1, flag2});
	}

	void Model::collect_necessary_gates (const aiger* aig, const aiger_symbol* as, const int as_size, 
	                                        std::set<unsigned>& exist_gates, std::set<unsigned>& gates, bool next)
	{
		for (int i = 0; i < as_size; i ++)
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
			    		outputs_[i] = true_;
			    	else if (is_false (as[i].lit))
			    		outputs_[i] = false_;
			    }
			}
			recursively_add (aa, aig, exist_gates, gates);	
		}
		
	}
	
	aiger_and* Model::necessary_gate (const unsigned id, const aiger* aig)
	{
		if (!is_true (id) && !is_false (id))
			return aiger_is_and (const_cast<aiger*> (aig), (id % 2 == 0) ? id : (id-1));
			
		return NULL;
	}
	
	void Model::recursively_add (const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates)
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
	
	/**
	 * @brief 
	 * 
	 * @param aa 
	 */
	void Model::add_clauses_from_equation (const aiger_and* aa)
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
	int Model::prime (const int id)
	{
		nextMap::iterator it = next_map_.find (abs (id));
		if (it == next_map_.end ())
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
	std::vector<int> Model::previous (const int id)
	{
	    vector<int> res;
	    reverseNextMap::iterator it = reverse_next_map_.find (abs (id));
		if (it == reverse_next_map_.end ())
		    return res; //not found
		res = it->second;
		if (id < 0)
		{
		    for (int i = 0; i < res.size (); i ++)
		        res[i] = -res[i];
		}
		return res;
	}
	
	/**
	 * @brief collect the previous vars from uc.  
	 * @note multiple latch may share a common next, therefore one lit may have several previous.
	 * @param uc 
	 */
	void Model::shrink_to_previous_vars (Cube& uc)
	{
		Cube tmp;
		for (int i = 0; i < uc.size (); i ++)
		{
		    vector<int> ids = previous (abs (uc[i]));
			if (ids.empty ())
			    continue;
			// 每一个pre都放进去。
			for (int j = 0; j < ids.size (); j ++)
				tmp.push_back ((uc[i] > 0) ? ids[j] : (-ids[j]));	
		}
		uc = tmp;
	}
	
	/**
	 * @brief collect the latch vars in uc.  
	 * @pre since we now set all the flags t obe at the beginning of the assumption, they shall be at the end.
	 * @param uc 
	 */
	void Model::shrink_to_latch_vars (Cube& uc)
	{
        if(!latch_var(abs(uc.back())))
            uc.pop_back();
        // for(auto i :uc)
        // {
        //     assert(latch_var(abs(i)));
        // }
	}

bool absIncr (int i, int j)
{
	return abs (i) < abs(j);
}

bool imply(const std::vector<int>& v1, const std::vector<int>& v2)
{
    std::unordered_set<int> marks;
    for(auto lit: v1)
        marks.insert(lit);
    for(auto lit: v2)
        if(marks.find(lit) == marks.end())
            return false;
    marks.clear();
    return true;
}

Cube negate(const Cube& cu)
{
    Cube res;
    res.reserve(cu.size());
    for(auto &i : cu)
        res.emplace_back(-i);
    return res;
}

}
