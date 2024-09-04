#ifndef carChecker_H
#define carChecker_H

#include "definition.h"
#include "implysolver.h"
#include "invsolver.h"
#include "mainsolver.h"
#include <assert.h>
#include "statistics.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stack>
#include <queue>
#include <fstream>
#include <map>

/**
 * @file carChecker.h
 * @author yibodong (prodongf@gmail.com) , jianwenli (lijwen2748@gmail.com)
 * @brief A checker, using Complementary Approximate Reachability algorithm.
 * @version 0.1.0
 * 
 * 
 */

namespace car
{
    extern Statistics CARStats; /// defined in main.cpp
    class Checker;

    class Checker
    {
    private:
        /**
         * @section options
         *
         */
        int inter_cnt = 0;
        bool rotate_enabled = false;
        bool inv_incomplete =false;
        bool uc_no_sort = false;
        int imply_decision = -1;
        int time_limit_to_restart = -1;
        int rememOption = 0;
        int impMethod = 0;
        int LOStrategy = 0;
        int convAmount = 0;
        int convMode=-1;
        int convParam=0;
        bool subStat = false;
        bool partial = false;
        bool restart_enabled = false;
        bool multi_solver = false;
        bool simp = false;

    public:
        /**
         * @section auxiliary data structures.
         */

        ///////////////////////////////////
        /// @subsection Propagation     ///
        ///////////////////////////////////

        enum PropModeEnum
        {
            PropNone        = 0,    // not activate.
            PropAlways      = 1,    // Do propagation check always
            PropShort       = 2,    // Only Propagate when it's a short UC
            PropContinue    = 3,    // If propagation succeed, continue for the next frame.
            PropShortCont   = 4,    // Short && Continue
            PropFresh       = 5,    // Only Propagate when it's at highest level.
            PropShortFresh  = 6,    // Only Short and Fresh.
            PropEarlyStop   = 7,    // Only early at this level.
            PropEarlyShortFresh = 8,// Early+Short+Fresh
        };
        int propMode = 0;
        int propParam = 0;

        // should only call after an UNSAT main call.
        void propagate(Cube& uc, int level);

    public:
        /**
         * @section auxiliary data structures.
         */

        ///////////////////////////////////
        /// @subsection mUC             ///
        ///////////////////////////////////

        /// we need to record whether it used another UC in the history.
        /// level -> array[INT_MAX]
        /// we use the bit to record whether it pushed another UC.
        std::unordered_map<int, unsigned short> conv_record;
        enum convModeEnum{
            ConvModeAlways = 0,
            ConvModeLow = 1,
            ConvModeHigh = 2,
            ConvModeStuck = 3,
            ConvModeRand = 4
        };
        // whether to calculate more UC.
        bool convTriggered(int level);
        void mUC(const Cube &uc, int level);
        
        ///////////////////////////////////
        /// @subsection container       ///
        ///////////////////////////////////
        ContainerEnum container_option;
        
        std::stack<Obligation> _container_stk;
        std::priority_queue<Obligation, std::vector<Obligation>, CompareItem> _container_pq;
        
        void containerClear(){
            switch (container_option)
            {
            case CONT_PRIQ:
                while(!_container_pq.empty())
                    _container_pq.pop();
                break;
            case CONT_STACK:
            default:
                while(!_container_stk.empty())
                    _container_stk.pop();
                break;
            }
        };
        void containerPush(const Obligation &o){
            switch (container_option)
            {
            case CONT_PRIQ:
                _container_pq.push(o);
                break;
            case CONT_STACK:
            default:
                _container_stk.push(o);
                break;
            }
        }
        void containerPop(){
            switch (container_option)
            {
            case CONT_PRIQ:
                _container_pq.pop();
                break;
            case CONT_STACK:
            default:
                _container_stk.pop();
                break;
            }
        }
        const Obligation & containerTopOrFront(){
            switch (container_option)
            {
            case CONT_PRIQ:
                return _container_pq.top();
            case CONT_STACK:
            default:
                return _container_stk.top();
            }
        }
        bool containerEmpty(){
            switch (container_option)
            {
            case CONT_PRIQ:
                return _container_pq.empty();
            case CONT_STACK:
            default:
                return _container_stk.empty();
            }
        }


        ///////////////////////////////////
        /// @subsection restart         ///
        ///////////////////////////////////
        /// remember Option: whether we shall remember something during restart
        enum rememberEnum{
            remem_None = 0,
            remem_O0 = 1,
            remem_short = 2,
            remem_Ok = 3,
            remem_Uk = 4,
        };

        /// last checker, which we can get information from.
        /// TODO: rewrite using copy-constructors
        Checker* last_chker;

        /// timer to calculate time consumption.
        clock_high sat_timer;

        ///////////////////////////////////
        /// @subsection assum ordering  ///
        ///////////////////////////////////
        enum LiteralOrderingEnum{
            LO_Classic = 0,     /// Intersection + Rotation for (1), Reverse it for (2), Random for the subsequent ones.
            LO_Rand = 1,        /// Just do random test.
            // LO_Fixpoint = 2,    /// Rotation for (1), {not_used, used} for later
            LO_Classic_Pos = 3, /// don't do for level == -1
        };

        inline void set_inter_cnt(int cnt)      { inter_cnt = cnt; }
        inline int  get_inter_cnt()             { return inter_cnt; }
        inline void set_rotate(bool rtt)        { rotate_enabled = rtt; }
        inline bool get_rotate()                { return rotate_enabled; }
        void get_inter_res(State*s, int level, std::vector<Cube>&);
        void get_rotate_res(State*s, int level, Cube& rcube, Cube& rest);
        std::vector<Cube> reorderAssumptions(State*s, int level, const std::vector<Cube>& inter, const Cube &rres, const Cube &rtmp);

        /// rotation calculation:
        std::vector<Cube> rotates;  /// corresponding to O[i]

        ///////////////////////////////////
        /// @subsection Expectation     ///
        ///////////////////////////////////
        
        
        /// Reuse the phase-saving strategy to implement a biased-SAT instead.


        ///////////////////////////////////
        /// @subsection Implication     ///
        ///////////////////////////////////

        /// how to calculate imply.
        enum ImpHowEnum{
            Imp_Manual      = 0,
            Imp_Solver      = 1,
            Imp_Sample      = 2,
            /// Imp_Sort    = 3, /// proved not to be useful. depricated.
            Imp_Exp         = 4,
            Imp_Thresh      = 5,
            /// Imp_MOM     = 6, /// proved not to be useful. depricated.
            Imp_Fresh       = 7, /// record the checking process of state.
            Imp_Bit         = 8, /// Use Bit Masks to represent UCs.
            Imp_BitFresh    = 9, /// Use Bit Masks, also record fresh indexes.
        };

        inline bool needImpSolver()
        {
            static std::unordered_set<int> __needImpSolver{Imp_Solver, Imp_Sample, Imp_Exp, Imp_Thresh};
            return __needImpSolver.count(impMethod)!=0;
        }

        /// level -> id -> freshIndex
        std::unordered_map<int, std::unordered_map <int, int>> impFreshRecord;

        std::vector<std::vector<UCMask>> ucs_masks;

        ///////////////////////////////////
        /// @subsection Garbage Collect ///
        ///////////////////////////////////
        /// the states to be cleared
        std::unordered_set<State *> clear_duties;
        /// add into clear duty, waiting for clear in the end.
        inline void clear_defer(State *s) { clear_duties.insert(s); };
        
        ///////////////////////////////////
        /// @subsection Invariant Check ///
        ///////////////////////////////////
        /// the map from O sequence to its minimal_level
        int fresh_levels;

    public:     
        /// @section counter-example
        std::ostream &out;
        /// record the prior state in this trail. This will be used in counter example printing.
        std::map<State *, State *> prior_in_trail_f;
        std::map<State *, State *> prior_in_trail_b;
        inline std::map<State *, State *> &whichPrior() { return backwardCAR ? prior_in_trail_b : prior_in_trail_f; }
        State *counter_start_f = nullptr;
        State *counter_start_b = nullptr;
        inline State *&whichCEX() { return backwardCAR ? counter_start_b : counter_start_f; }
        inline State *&otherCEX() { return backwardCAR ? counter_start_f : counter_start_b; }

    public:
        /// init all the solvers.
        explicit Checker(Problem *model, const OPTIONS& opt, std::ostream &out, Checker* last_chker);
        /// delete all the solvers.
        ~Checker();

        Checker(const Checker &) = delete;
        const Checker &operator=(const Checker &) = delete;

        /**
         * @brief entrance for CAR checking
         * 
         * @return ResEnum
         */
        RESEnum check();

    public:
        const bool backwardCAR;
        const int bad_;
        Problem         *model_         = nullptr;
        InvSolver       *inv_solver     = nullptr;
        /// @brief if one solver shared among all the frames, use this
        MainSolver      *main_solver    = nullptr;
        /// else, use this.
        std::vector<MainSolver*> main_solvers = {};
        // Propagation Solver.
        MainSolver      *prop_solver    = nullptr;
        /**
         * @brief Get the Main Solver of `level`. 
         * @note if level < 0, reuse that of level 0.
         */
        MainSolver* getMainSolver(int level);

        USequence Uf, Ub; /// Uf[0] is not explicitly constructed
        OSequence Onp, OI;
        inline USequence &whichU() { return backwardCAR ? Ub : Uf;}
        inline USequence &otherU() { return backwardCAR ? Uf : Ub;}
        inline OSequence &whichO() { return backwardCAR ? Onp :OI;}
        
        inline OFrame& whichFrame(int index){ assert(index <= whichO().size()); return whichO().at(index);}
        // minus one: for Otmp
        inline int OSize() {return whichO().size() - 1;}
        /// used in picking state randomly
        std::vector<std::pair<State *, int>> Uset;

    private:
        /// entrance for CAR
        RESEnum car();

        /**
         * @brief Check for immediate safe or unsafe
         *
         */
        RESEnum trivialCheck();

        /**
         * @brief Searching procedure of bi-car. Enumerate states in U and use the OSequence as the guide.
         *
         * @return true : successfully reached O[0]. which means a cex is found.
         * @return false : all states in present U has been checked. No cex is found.
         */
        RESEnum trySAT();

        /// @brief print evidence.
        void print_evidence() const;

        /**
         * @brief iteratively pick state from the given sequence.
         * @return State*
         * 	then this state will be marked as visited at this round.
         * @return nullptr every state has been tried.
         *  then all the markers will be cleaned. Therefore, this method can be reused for another round.
         * @post Anytime it ends the iterative round earlier(before nullptr), it must has found the counter example in trySAT. Therefore, there will be no picking in the future.
         */
        State *pickState();
        int pickStateLastIndex = -1;

        /**
         * @brief Get the solution from SAT solver, update it to U sequence.
         *
         * @return State* : the state retrieved.
         */
        State *getModel(State *s, int level);

        /**
         * @brief SAT_Assume(assum, clauses)
         *
         * @return true
         * @return false
         */
        bool satAssume(State *, int, bool &safe_reported);

        /**
         * @brief Interface for cleaning.
         *
         */
        void clean();

        /**
         * @brief first create a clause with target uc and flag of the corresponding O[dst], then add this clause to main solver or start solver, depending on the level.
         *
         * @param uc
         * @param O
         * @param dst_level
         */
        void addUCtoSolver(Cube &uc, int dst_level);

        /**
         * @brief init special sequences: Uf, Ub, Oi, Onp
         * Will call ImmediateCheck.
         */
        RESEnum initSequence();

        /**
         * @brief Check SAT_ASSUME(I, ~P).
         *
         * to generate uc from init
         * ~p should be in ~uc
         * use uc to initialize O[0] is suitable.
         *
         */
        RESEnum immediateCheck(State *from, OFrame &UCs);


        /**
         * @brief Whether this state is blocked in this O frame, namely O[frame_level] (or Otmp, if frame_level == O.size)
         *
         * @param s
         * @param frame_level
         * @param O
         * @return true
         * @return false
         */
        bool blockedIn(State *s, const int frame_level);

        /**
         * @brief Use `blocked in` to iterate, from min to max, to find the minimal level where this state is not blocked.
         *
         * @param s
         * @param min
         * @param max
         * @param O
         * @return int
         */
        int minNOTBlocked(State *s, const int min, const int max);

        /**
         * @brief Check whether there exists an invariant in this sequence. Which means the center state of this OSequence is not reachable, and should be later blocked
         *
         * @param O the sequence to be checked
         * @return true : invariant is found
         * @return false : no invariant is found
         */
        bool InvFound();
        bool InvFoundAt(int check_level, int minimal_update_level);

    };
}

#endif
