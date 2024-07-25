/**
 * @file definition.h
 * @author yibodong (prodongf@gmail.com)
 * @brief The basic data structures
 * @version 0.1.0
 * @date 2024-07-23
 * 
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
#include <string>
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
    // TODO: change to std::array.
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
			explicit State(bool _isnegp) : _isnegp(_isnegp) { assert(_isnegp); negp_state=this; }
            const State *get_negp() const { return negp_state; }

            // Construct a state with only latches, inputs to be set later.
			explicit State(const Assignment &latches) : _latches(latches), _isnegp(false) {}
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


    // ##################################################
    // #####                 Problem               ######
    // ##################################################
class Problem {
private:
    int num_inputs_;
	int num_latches_;
	int num_ands_;
	int num_constraints_;
	int num_outputs_;

    int max_id_;    // maximum used id in the model. Also preserve two more ids for TRUE (max_id_ - 1) and FALSE (max_id_)
	int true_;      // id for true
	int false_;     // id for false
	
	typedef std::vector<int> vect;
	typedef std::vector<vect> Clauses;

	std::unordered_set<unsigned> trues_;  //vars evaluated to be true, and their negation is false
	vect init_;   //initial values for latches
	vect outputs_; //output ids
	vect constraints_; //constraint ids

public:
    typedef std::unordered_map<int, int> nextMap;
	typedef std::unordered_map<int, std::vector<int> > reverseNextMap;
	nextMap next_map_;  // map from latches to their next values
	reverseNextMap reverse_next_map_;  // map from the next values of latches to latches
	                                   
public:
    // construct from an AIGER.
	explicit Problem (aiger*);
	~Problem () {}
	
    // get the 'next' of this literal.
	int prime (const int) const;
    // get those literals whose next is this.
	std::vector<int> previous (const int) const;
	
	bool state_var (const int id) const {return (id >= 1) && (id <= num_inputs_+num_latches_);}
	bool latch_var (const int id) const {return (id >= num_inputs_+1) && (id <= num_inputs_+num_latches_);}

    inline int num_inputs()         const { return num_inputs_; }
    inline int num_latches()        const { return num_latches_; }
    inline int num_ands()           const { return num_ands_; }
    inline int num_constraints()    const { return num_constraints_; }
    inline int num_outputs()        const { return num_outputs_; }

    inline int max_id()             const { return max_id_; }
    inline int true_id()           const { return true_;}
	inline int false_id()          const { return false_;}
    inline int size()               const { return cls_.size(); }
    inline int output(const int id) const { return outputs_[id]; }

    inline int outputs_start()      const { return outputs_start_; }
    inline int common_next_start()  const { return common_next_start_; }
    inline int latches_start()      const { return latches_start_; }
    
    inline const std::vector<int> &element(const int id) const { return cls_.at(id); }
    // return initial values for latches
    inline const std::vector<int> &init() const { return init_; }

    void shrink_to_previous_vars    (Cube &cu) const;
    void shrink_to_latch_vars       (Cube &cu) const;
	
public:
	Clauses cls_;   //set of clauses, it contains several parts:
	                //(1) clauses for constraints, i.e. those before position outputs_start_;
	                //(2) same next have same previous, initialized to 0.
					//(3) clauses for outputs, i.e. those before position latches_start_;
	                //(4) clauses for latches 
					//(5) clause for encoding our FALSE / TRUE
	
	int common_next_start_;
	int outputs_start_; //the index of cls_ to point the start position of outputs
	int latches_start_; //the index of cls_ to point the start position of latches

	inline bool is_true (const unsigned id) const
	{
		return (id == 1) || (trues_.find (id) != trues_.end ());
	}
	
	inline bool is_false (const unsigned id) const
	{
		return (id == 0) || (trues_.find ((id % 2 == 0) ? (id+1) : (id-1)) != trues_.end ());
	}
	
	inline int car_var (const unsigned id) const
	{
		assert (id != 0);
		return ((id % 2 == 0) ? (id/2) : -(id/2));
	}
	
	inline vect clause (int id) const
	{
		return (vect){id};
	}
	
	inline vect clause (int id1, int id2) const
	{
		return (vect){id1, id2};
	}
	
	inline vect clause (int id1, int id2, int id3) const
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
	
	void set_init (const aiger* aig);
	void set_constraints (const aiger* aig);
	void set_outputs (const aiger* aig);

    // collect those that are trivially constant
	void collect_trues (const aiger* aig);
    // previous -> next
	void create_next_map (const aiger* aig);
    // next -> previous(s)
	void insert_to_reverse_next_map (const int index, const int val);

	void collect_necessary_gates (const aiger* aig, const aiger_symbol* as, const int as_size, std::set<unsigned>& exist_gates, std::set<unsigned>& gates, bool next = false);
	aiger_and* necessary_gate (const unsigned id, const aiger* aig);
	void recursively_add (const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates);

	void add_clauses_from_equation (const aiger_and* aa);
	void create_constraints_for_latches ();
	void create_clauses (const aiger* aig);
};


    // ##################################################
    // #####            Tool Functions             ######
    // ##################################################

    /**
     * @brief whether v2 is contained in v1
     * @pre v2 and v1 can be unordered. 
     */
    inline bool imply_heavy (const std::vector<int>& v1, const std::vector<int>& v2)
    {
        std::unordered_set<int> marks(v1.begin(), v1.end());
        for(auto lit: v2)
            if(marks.find(lit) == marks.end())
                return false;
        return true;
    }

    inline bool absIncr (int i, int j) {return abs (i) < abs(j);}

    inline std::vector<int> negate(const std::vector<int>& cu) {
        std::vector<int> res = cu;
        for (auto &i : res)
            i = -i;
        return res;
    }

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


    // ##################################################
    // #####               CLI OPTIONS             ######
    // ##################################################
    struct OPTIONS{
        // use BMC rather than CAR
        bool bmc = false;
        // forward CAR / backward CAR
        bool forward = false;
        // use propagation
        bool propagate = false;
        // partial state enabled during forward-CAR.
        bool partial = false;
        // enable rotate
        bool enable_rotate = false;
        // do not check INVARIANT ==> incomplete
        bool inv_incomplete = false;
        // keep UC not sorted
        bool raw_uc = false;
        // #(UC) for intersection
        int inter_cnt = 0;

        // mUC mode
        int convMode = -1;
        // mUC param
        int convParam = 0;
        // #(UC)
        int convAmount = 0;
        // time limit to restart
        int time_limit_to_restart = -1;
        // remember option during restart
        int rememOption = 0;
        // literal ordering strategy
        int LOStrategy = 0;
        // method for calculate imply
        int impMethod = 0;
        // option for subsumption. NOTE: this is heavy!
        bool subStat = false;

        std::string inputPath;
        std::string outputPath;
    };
}

#endif