/**
 * @file definition.h
 * @author yibodong (prodongf@gmail.com)
 * @brief The basic data structures
 * @version 0.1.0
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef DEFINITION_H
#define DEFINITION_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map> 
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <assert.h>
#include "statistics.h"
extern "C" {
#include "aiger.h"
}
#include "assert.h"

namespace car
{
    // ##################################################
    // #####           Basic Data Strucutre        ######
    // ##################################################


    // Problem ::= {Model, Property, (Optional) Constraints }
    // where:
    //      Model ::= BTS(Boolean Transition System) === {V, I, T}, where
    //          V = All the Variables
    //          I = Initial States
    //          T = Transition Relationship.
    //      Property ::= CNF to be satisfied
    //      Constraints ::= CNF that are always satisfied
    class Problem;
    // state ::= {inputs, latches}
    class State;

    // full assignment of a state
    using Assignment = std::vector<int>;
    // a CNF of literals
    using Cube = std::vector<int>; // TODO: change to LIT?
    // a DNF of literals
    using Clause = std::vector<int>; 
    // Many Cubes, representing an Over-Approximation
	using OFrame = std::vector<Cube>;
    // Many OFrames, representing a trial of Over-Appriximations, where O_{i+1} is the Over-appximation of R{O_i} (or R^{-1}). 
    using OSequence = std::vector<OFrame>;
    // Many States, representing an Under-Approximation.
	using USequence = std::vector<State *>;
    // <state, depth, target level> to be verified.
	using Obligation = std::tuple<car::State *, int, std::size_t>;
	// how to compare.
    struct CompareItem {
		bool operator()(const Obligation &a, const Obligation &b) const {
            // use the third member to sort.
			return int(std::get<2>(a)) >= int(std::get<2>(b)); 
		}
	};

    // ##################################################
    // #####                State                  ######
    // ##################################################

	class State
	{
        // members of a state
        private:
            // latch values
            Assignment _latches;
            // input values
			Assignment _inputs;
        
        public:
			// whether it is the special state 'neg P'. 
            // NOTE: At present, we do not unfold 'neg P', but use SAT solver to retrieve it's approximation (Follow simpleCAR). So we use one special state to represent it.
            const bool _isnegp = false;
			static const State* negp_state;

            // the problem imported from AIGER
            static const Problem* model_;
			static const aiger* aig_;
            
            // size information about the problem
			static int num_inputs_;
			static int num_latches_;
            static void set_num_inputs_and_latches (const int n1, const int n2) 
			{
				num_inputs_ = n1;
				num_latches_ = n2;
			}

		public:
            // Construct a special state 'neg P', whose id is -1.
			State(bool _isnegp) : _isnegp(_isnegp) { assert(_isnegp); negp_state=this; }
            const State *get_negp() const { return negp_state; }

            // Construct a state with only latches, inputs to be set later.
			State(const Assignment &latches) : _latches(latches), _isnegp(false) {}
            inline void set_inputs(const Assignment &st) { _inputs = st; }

            // Construct a state with inputs and latches.
			State(const Assignment &inputs , const Assignment &latches): _latches(latches), _inputs(inputs), _isnegp(false) {}

            // do nothing.
			~State() {}

			inline const Assignment& get_latches() const { return _latches; }
			inline const Assignment& get_inputs() const { return _inputs; }
			inline int latch_size() const { return _latches.size(); }
			inline int latch_at(int i) const { return _latches.at(i); }
			std::string get_inputs_str() const;
			std::string get_latches_str() const;

            // whether current state is blocked by this UC.
            // this is a basic subsumption test: check whether each literal in this UC appears in this state.
            // NOTE: latch values are stored in order ==> we could fetch-by-index.
			bool imply(const Cube &cu) const;

			// get the intersection with `cu`.
			Cube intersect(const Cube &cu) const;

		public:
			// calculate what is already known about next latches according to present latch values.
            // TODO: implement it.
			Cube probe();
	};


//elements in v1, v2 are in order
//check whether v2 is contained in v1 
bool imply (const std::vector<int>& v1, const std::vector<int>& v2);

bool absIncr (int i, int j);

Cube negate(const Cube& cu);

template <typename T>
void shuffle(std::vector<T>& vec) {
    #ifdef RANDSEED
        int seed = RANDSEED;
    #else
        int seed = 1;
    #endif
    std::mt19937 g(seed);
    std::shuffle(vec.begin(), vec.end(), g);
}


class Problem {
public:
	Problem (aiger*, const bool verbose = false);
	~Problem () {}
	
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