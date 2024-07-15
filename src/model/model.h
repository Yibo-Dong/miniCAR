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

#ifndef MODEL_H
#define MODEL_H

extern "C" {
#include "aiger.h"
}
#include "assert.h"
#include "basic_data.h"
#include <unordered_map>
#include <unordered_set>
#include <set>


namespace car {

class Model {
public:
	Model (aiger*, const bool verbose = false);
	~Model () {}
	
	int prime (const int);
	std::vector<int> previous (const int);
	
	bool state_var (const int id)  {return (id >= 1) && (id <= num_inputs_+num_latches_);}
	bool latch_var (const int id)  {return (id >= num_inputs_+1) && (id <= num_inputs_+num_latches_);}

	inline int num_inputs() const { return num_inputs_; }
	inline int num_latches() const { return num_latches_; }
	inline int num_ands() const { return num_ands_; }
	inline int num_constraints() const { return num_constraints_; }
	inline int num_outputs() const { return num_outputs_; }
	inline int max_id() const { return max_id_; }
	inline int outputs_start() const { return outputs_start_; }
	inline int common_next_start() const { return common_next_start_; }
	inline int latches_start() const { return latches_start_; }
	inline int size() const { return cls_.size(); }
	std::vector<std::vector<int>> output_clauses()
	{
		// FIXME: Is this right?
		assert(latches_start_ > outputs_start_);
		return std::vector<car::Cube>(cls_.begin() + outputs_start_, cls_.begin() + latches_start_);
	};
	inline std::vector<int>& element (const int id) {return cls_[id];}
	inline int output (const int id)const  {return outputs_[id];}
	
	// return initial values for latches
	inline std::vector<int>& init () {return init_;}
	
	void shrink_to_previous_vars (Cube& cu);
	void shrink_to_latch_vars (Cube& cu);
	
	inline int true_id () const {return true_;}
	inline int false_id () const {return false_;}
	
private:
	//members
	bool verbose_;
	
	
	int num_inputs_;
	int num_latches_;
	int num_ands_;
	int num_constraints_;
	int num_outputs_;

public:
	int max_model_id_;
private:
	int max_id_;  //maximum used id in the model. Also preserve two more ids for TRUE (max_id_ - 1) and FALSE (max_id_)
	
	int true_;  //id for true
	int false_;  //id for false
	
	typedef std::vector<int> vect;
	typedef std::vector<vect> Clauses;
	
	vect init_;   //initial values for latches
	vect outputs_; //output ids
	vect constraints_; //constraint ids
public:
	Clauses cls_;  //set of clauses, it contains three parts:
	                //(1) clauses for constraints, i.e. those before position outputs_start_;
	                //(2) same next have same previous, initialized to 0.
					//(3) clauses for outputs, i.e. those before position latches_start_;
	                //(4) clauses for latches 
					//(5) clause for encoding our FALSE / TRUE
	
	int common_next_start_;
	int outputs_start_; //the index of cls_ to point the start position of outputs
	int latches_start_; //the index of cls_ to point the start position of latches

	
	typedef std::unordered_map<int, int> nextMap;
	nextMap next_map_;  //map from latches to their next values
	typedef std::unordered_map<int, std::vector<int> > reverseNextMap;
	reverseNextMap reverse_next_map_;  //map from the next values of latches to latches
	                                   //BE careful the situation when next (a) = c and next (b) = c!!
	
	std::unordered_set<unsigned> trues_;  //vars evaluated to be true, and their negation is false
	
	
	//functions

	inline bool is_true (const unsigned id)
	{
		return (id == 1) || (trues_.find (id) != trues_.end ());
	}
	
	inline bool is_false (const unsigned id)
	{
		return (id == 0) || (trues_.find ((id % 2 == 0) ? (id+1) : (id-1)) != trues_.end ());
	}
	
	inline int car_var (const unsigned id)
	{
		assert (id != 0);
		return ((id % 2 == 0) ? (id/2) : -(id/2));
	}
	
	inline vect clause (int id)
	{
		return (vect){id};
	}
	
	inline vect clause (int id1, int id2)
	{
		return (vect){id1, id2};
	}
	
	inline vect clause (int id1, int id2, int id3)
	{
		return (vect){id1, id2, id3};
	}
	
	inline void set_outputs_start ()
	{
	    outputs_start_ = cls_.size ();
	}
	
	inline void set_latches_start ()
	{
	    latches_start_ = cls_.size ();
	}

	inline void set_common_next_start()
	{
		common_next_start_ = cls_.size();
	}
	
	void collect_trues (const aiger* aig);
	void create_next_map (const aiger* aig);
	void create_clauses (const aiger* aig);
	void collect_necessary_gates (const aiger* aig, const aiger_symbol* as, const int as_size, std::set<unsigned>& exist_gates, std::set<unsigned>& gates, bool next = false);
	aiger_and* necessary_gate (const unsigned id, const aiger* aig);
	void recursively_add (const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates);
	void add_clauses_from_equation (const aiger_and* aa);
	void set_init (const aiger* aig);
	void set_constraints (const aiger* aig);
	void set_outputs (const aiger* aig);
	void insert_to_reverse_next_map (const int index, const int val);
	void create_constraints_for_latches ();
	
	
};

}

#endif
