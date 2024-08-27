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
#define NOT_M_FLAG -2
#define NOT_P_FLAG -3
namespace car
{
    /// ##################################################
    /// #####           Basic Data Strucutre        ######
    /// ##################################################


    /// Problem to be check
    class Problem;
    /// state ::= {inputs, latches}
    class State;
    
    /// full assignment of a state
    using Assignment = std::vector<int>;/// TODO: change to std::array.
    /// a CNF of literals
    using Cube = std::vector<int>; /// TODO: change to LIT?
    /// a DNF of literals
    using Clause = std::vector<int>; 
    /// Many Cubes, representing an Over-Approximation
	using OFrame = std::vector<Cube>;
    /// Many OFrames, representing a trial of Over-Appriximations, where O_{i+1} is the Over-appximation of R{O_i} (or R^{-1}). 
    using OSequence = std::vector<OFrame>;
    /// Many States, representing an Under-Approximation.
	using USequence = std::vector<State *>;
    /// <state, depth, target level> to be verified.
	using Obligation = std::tuple<car::State *, int, std::size_t>;
	/// how to compare.
    struct CompareItem {
		bool operator()(const Obligation &a, const Obligation &b) const {
            /// use the third member to sort.
			return int(std::get<2>(a)) >= int(std::get<2>(b)); 
		}
	};

    /// ##################################################
    /// #####                State                  ######
    /// ##################################################

	class State
	{
        /// members of a state
        private:
            /// latch values
            Assignment _latches;
            /// input values
			Assignment _inputs;
    
            // size of the problem, not the state.
            static int _num_inputs;
			static int _num_latches;
            static int _maxid;
            
            /// whether it is the special state 'neg P'. 
            /// NOTE: At present, we do not unfold 'neg P', but use SAT solver to retrieve it's approximation (Follow simpleCAR). So we use one special state to represent it.
            const bool _isnegp = false;
            static const State* negp_state;
    
        public:
            /// ID, used to mark during blockedIn().
            const int id;

            void setInputs(const Assignment &st) { _inputs = st; }
            int getProblemInputSize() const { return _num_inputs; }
            int getProblemLatchSize() const { return _num_latches; }
            const Assignment &getInputs() const { return _inputs; }
            const Assignment &getLatches() const { return _latches; }
            int getLatchAt(int i) const
            {
                assert(i >= 0 && i < _latches.size());
                return _latches.at(i);
            }
            std::string getInputsStr() const;
            std::string getLatchesStr() const;

            bool isNegP() const { return _isnegp; }
            const State *getNegP() const { return negp_state; }
            bool isPartial() const { return _latches.size() != _num_latches; }

            static void setProblemSize(const int ni, const int nl)
            {
				_num_inputs = ni;
				_num_latches = nl;
			}

		public:
            /// Construct a special state 'neg P', whose id is -1.
			explicit State(bool _isnegp) : _isnegp(_isnegp),id(-1) { assert(_isnegp); negp_state=this; }
            /// Construct a state with only latches, inputs to be set later.
			explicit State(const Assignment &latches) : _latches(latches),id(++_maxid), _isnegp(false) {}
            /// Construct a state with inputs and latches.
			State(const Assignment &inputs , const Assignment &latches): _latches(latches), _inputs(inputs),id(++_maxid), _isnegp(false) {}
            /// do nothing.
			~State() {}

            /// whether current state is blocked by this UC.
            /// this is a basic subsumption test: check whether each literal in this UC appears in this state.
			bool imply(const Cube &cu) const;

			/// get the intersection with `cu`.
			Cube intersect(const Cube &cu) const;

		public:
			/// calculate what is already known about next latches according to present latch values.
            /// TODO: implement it.
			Cube probe();
	};


    /// ##################################################
    /// #####                 Problem               ######
    /// ##################################################

/**
 * @brief     
 *  Problem ::= {Model, Property, (Optional) Constraints }
 *  where:
 *      Model ::= BTS(Boolean Transition System) === {V, I, T}, where
 *          V = All the Variables,
 *          I = Initial States,
 *          T = Transition Relationship;
 *      Property ::= CNF to be satisfied,
 *      Constraints ::= CNF that are always satisfied.
 * 
 */
class Problem {
    /**
     * In AIG, lits are represented by unsigned int, with 0 reserved for FALSE:
     * 0: False. 1: True.
     * 2: v0, 3: ~v0, 4: v1, 5: ~v1, ... 
     * 
     */
private:
    aiger *aig;
    int _num_inputs;
    int _num_latches;
    int _num_ands;
    int _num_constraints;
    int _num_outputs;


public:
    int _max_id; /// maximum used id
    int _true;   /// id for true
    int _false;  /// id for false

	std::unordered_set<unsigned> _true_set; ///vars evaluated to be true
    Cube _init;                    /// initial values for latches
    std::vector<int> _outputs;     /// output ids
    std::vector<int> _constraints; /// constraint ids

	std::unordered_map<int, int> _next_map;  /// map from latches to their next values
	std::unordered_map<int, std::vector<int> > _reverse_next_map;  /// map from the next values of latches to latches
	                                   
public:
    /// construct from an AIGER.
	Problem (aiger*, bool unroll);
	~Problem () {}
	
    /// get the 'next' of this literal.
	int prime (const int) const;
    /// get those literals whose next is this.
	std::vector<int> previous (const int) const;
	
    bool input_var (const int id) const {return abs(id) >= 1 && (abs(id) <= _num_inputs);}
	bool state_var (const int id) const {return (abs(id) >= 1) && (abs(id) <= _num_inputs+_num_latches);}
	bool latch_var (const int id) const {return (abs(id) >= _num_inputs+1) && (abs(id) <= _num_inputs+_num_latches);}

    inline int num_inputs()         const { return _num_inputs; }
    inline int num_latches()        const { return _num_latches; }
    inline int num_ands()           const { return _num_ands; }
    inline int num_constraints()    const { return _num_constraints; }
    inline int num_outputs()        const { return _num_outputs; }

    inline int max_id()             const { return _max_id; }
    inline int true_id()           const { return _true;}
	inline int false_id()          const { return _false;}
    inline int size()               const { return cls_.size(); }
    inline int output(const int id) const { return _outputs[id]; }

    inline int outputs_start()      const { return outputs_start_; }
    inline int common_next_start()  const { return common_next_start_; }
    inline int latches_start()      const { return latches_start_; }
    
    inline const std::vector<int> &element(const int id) const { return cls_.at(id); }
    /// return initial values for latches
    inline const std::vector<int> &init() const { return _init; }

    void shrink_to_previous_vars    (Cube &cu) const;
    // shrink the uc to its 'previous' version.
    // note: maybe not unique. just take the first one.
    void transform_uc_to_previous      (Cube& uc) const;
    void shrink_to_latch_vars       (Cube &cu) const;
	
public:
    /**
     * set of clauses, it contains several parts: 
     * (1) clauses for constraints, i.e. those before position outputs_start_;
     * (2) same next have same previous, initialized to 0.
     * (3) clauses for outputs, i.e. those before position latches_start_;
     * (4) clauses for latches 
     * (5) clause for encoding our FALSE / TRUE
     * 
     */
	std::vector<std::vector<int>> cls_;   
	
	int common_next_start_;
    ///the index of cls_ to point the start position of outputs
	int outputs_start_; 
    ///the index of cls_ to point the start position of latches
	int latches_start_; 

	inline bool is_true (const unsigned id) const
	{
		return (id == 1) || (_true_set.find (id) != _true_set.end ());
	}
	
	inline bool is_false (const unsigned id) const
	{
		return (id == 0) || (_true_set.find ((id % 2 == 0) ? (id+1) : (id-1)) != _true_set.end ());
	}
	
	inline int car_var (const unsigned id) const
	{
		assert (id != 0);
		return ((id % 2 == 0) ? (id/2) : -(id/2));
	}
	
	inline std::vector<int> clause (int id) const
	{
		return std::vector<int>{id};
	}
	
	inline std::vector<int> clause (int id1, int id2) const
	{
		return std::vector<int>{id1, id2};
	}
	
	inline std::vector<int> clause (int id1, int id2, int id3) const
	{
		return std::vector<int>{id1, id2, id3};
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

    /// collect those that are trivially constant
	void collectTrue (const aiger* aig);
    /// previous -> next
	void create_next_map (const aiger* aig);
    /// next -> previous(s)
	void insert_to_reverse_next_map (const int index, const int val);

	void collect_necessary_gates (const aiger* aig, const aiger_symbol* as, const int as_size, std::set<unsigned>& exist_gates, std::set<unsigned>& gates, bool next = false);
	aiger_and* necessary_gate (const unsigned id, const aiger* aig);
    
    // deprecated. Recursion causes stack overflow on large cases.
	void recursively_add (const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates);
    void iteratively_add (const aiger_and* aa, const aiger* aig, std::set<unsigned>& exist_gates, std::set<unsigned>& gates);

	void add_clauses_from_equation (const aiger_and* aa);
	void create_constraints_for_latches ();
	void create_clauses (const aiger* aig);

    void unroll_prime(const aiger* aig);
};


    /// ##################################################
    /// #####            Tool Functions             ######
    /// ##################################################

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

    /**
     * @brief Sort according to abs value.
     */
    inline bool absIncr (int i, int j) {return abs (i) < abs(j);}

    /**
     * @brief Get the negation of a CNF, which is a DNF.
     */
    inline std::vector<int> negateCube(const std::vector<int>& cu) {
        std::vector<int> res(cu.size());
        std::transform(cu.begin(), cu.end(), res.begin(), [](int i) { return -i; });
        return res;
    }

    /**
     * @brief Shuffle std::vector, using Mersenne Twister.
     */
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


    /// ##################################################
    /// #####               CLI OPTIONS             ######
    /// ##################################################
    /// Options of the checker.
    struct OPTIONS{
        /// use BMC rather than CAR
        bool bmc = false;
        /// forward CAR / backward CAR
        bool forward = false;

        /// partial state enabled during forward-CAR.
        bool partial = false;
        /// enable rotate
        bool enable_rotate = false;
        /// do not check INVARIANT ==> incomplete
        bool inv_incomplete = false;
        /// keep UC not sorted
        bool raw_uc = false;
        /// #(UC) for intersection
        int inter_cnt = 0;

        /// mUC mode
        int convMode = -1;
        /// mUC param
        int convParam = 0;
        /// #(UC)
        int convAmount = 0;
        /// time limit to restart
        int time_limit_to_restart = -1;
        /// remember option during restart
        int rememOption = 0;
        /// literal ordering strategy
        int LOStrategy = 0;
        /// method for calculate imply
        int impMethod = 0;
        /// option for subsumption. NOTE: this is heavy!
        bool subStat = false;
        /// use a seperate solver for each frame.
        bool multi_solver = false;

        /// use propagation
        int propMode = 0;
        int propParam = 0;

        // whether to simplify the CNF, using another simpSolver.
        bool simplifyCNF = false;

        // unroll primed version. Must ON for CONSTRAINTS. Think about it: Whether we should used bad' instead of O0?
        bool unrollPrime = false;

        std::string inputPath;
        std::string outputPath;
    };

    /// The status of checking.
    enum RESEnum{
        RES_UNSAFE  = 0,
        RES_SAFE    = 1,
        RES_UNKNOWN = 2,
        RES_FAIL    = 3,
        RES_RESTART = 4,
    };
}

#endif