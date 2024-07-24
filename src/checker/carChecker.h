#ifndef carChecker_H
#define carChecker_H

#include "definition.h"
#include "invsolver.h"
#include "startsolver.h"
#include "mainsolver.h"
#include "newpartialsolver.h"
#include <assert.h>
#include "statistics.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <fstream>
#include <map>

namespace car
{
    extern Statistics CARStats; // defined in main.cpp
    class Checker;

    class Checker
    {
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

        // we need to record whether it used another UC in the history.
        // level -> array[INT_MAX]
        // we use the bit to record whether it pushed another UC.
        std::unordered_map<int, unsigned short> conv_record;
        enum convModeEnum{
            ConvModeAlways = 0,
            ConvModeLow = 1,
            ConvModeHigh = 2,
            ConvModeStuck = 3,
            ConvModeRand = 4
        };

        enum LiteralOrderingEnum{
            LO_Classic = 0,     // Intersection + Rotation for (1), Reverse it for (2), Random for the subsequent ones.
            LO_Rand = 1,        // Just do random test.
            LO_Fixpoint = 2,    // Rotation for (1), {not_used, used} for later
            LO_Exp = 3,         // Reuse the phase-saving strategy to implement a biased-SAT instead.
        };

        // remember Option: whether we shall remember something during restart
        enum rememberEnum{
            remem_None = 0,
            remem_O0 = 1,
            remem_short = 2,
            remem_Ok = 3,
            remem_Uk = 4,
        };
        // last checker, which we can get information from.
        Checker* last_chker;


        // how to calculate imply.
        enum ImpHowEnum{
            Imp_Manual = 0,
            Imp_Solver = 1,
            Imp_Sample = 2,
            // Imp_Sort = 3, // proved not to be useful. depricated.
            Imp_Exp = 4,
            Imp_Thresh = 5,
            // Imp_MOM = 6, // proved not to be useful. depricated.
            Imp_Fresh = 7,
        };

        // level -> id -> freshIndex
        std::unordered_map<int, std::unordered_map <int, int>> impFreshRecord;

        /**
         * @section rotation technique
         *
         */
        // corresponding to O[i]
        std::vector<Cube> rotates;
        // corresponding to Otmp
        Cube rotate;



         
        clock_high sat_timer;

        

        inline void set_inter_cnt(int cnt) { inter_cnt = cnt; }
        inline int get_inter_cnt() { return inter_cnt; }
        inline void set_rotate(bool rtt) { rotate_enabled = rtt; }
        inline int get_rotate() { return rotate_enabled; }

        bool restart_enabled = false;

        bool importO = false;
        

    public:
        bool ppstoped = false;
        

    public:
        Checker(Problem *model, const OPTIONS& opt, std::ostream &out, Checker* last_chker);

        /**
         * @brief Destroy the Checker object
         * delete all the solvers.
         */
        ~Checker();

        // entrance for CAR checking
        bool check();

    private:
        // entrance for CAR
        bool car();

        /**
         * @brief Check for immediate safe or unsafe
         *
         * @param out where to print
         * @param res the check result
         * @return true : Bad == True / False
         * @return false
         */
        bool trivialCheck(bool &res);

        /**
         * @brief Searching procedure of bi-car. Enumerate states in U and use the OSequence as the guide.
         *
         * @param O the O sequence with a singleton (a state that is picked from the U sequence) making up its level 0.
         * @param forward whether the Transformation relationship should be used forward or backward
         * @return true : successfully reached O[0]. which means a cex is found.
         * @return false : all states in present U has been checked. No cex is found.
         */
        bool trySAT(OSequence *O, bool forward, bool &safe_reported);

    public:
        Problem *model_;
        bool backward_first;
        int bad_;
        StartSolver *start_solver;

        // mapping from state to its corresponding o sequence.
        std::unordered_map<const State *, OSequence *> SO_map;
        // the main solver shared.
        MainSolver *main_solver;
        // the partial solver shared.
        PartialSolver *partial_solver;

        // the map from O sequence to its minimal_level
        std::unordered_map<const OSequence *, int> fresh_levels;
        USequence Uf, Ub; // Uf[0] is not explicitly constructed
        OSequence Onp, OI;
        // used in picking state randomly
        std::vector<std::pair<State *, int>> Uset;

    public:
        /**
         * @section Counter Example
         *
         */

        inline USequence &whichU() 
        {
            return backward_first ? Ub : Uf;
        }

        inline USequence &otherU()
        {
            return backward_first ? Uf : Ub;
        }

        // record the prior state in this trail. This will be used in counter example printing.
        std::map<State *, State *> prior_in_trail_f;
        std::map<State *, State *> prior_in_trail_b;
        inline std::map<State *, State *> &whichPrior()
        {
            return backward_first ? prior_in_trail_b : prior_in_trail_f;
        }

        State *counter_start_f = nullptr;
        State *counter_start_b = nullptr;
        inline State *&whichCEX()
        {
            return backward_first ? counter_start_b : counter_start_f;
        }
        inline State *&otherCEX()
        {
            return backward_first ? counter_start_f : counter_start_b;
        }
        std::ostream &out;

        /**
         * @brief print evidence to out.
         *
         * @param out
         */
        void print_evidence() const;

    private:
        /**
         * @section Clear duties
         *
         */

        // the states to be cleared
        std::unordered_set<State *> clear_duties;
        // add into clear duty, waiting for clear in the end.
        inline void clear_defer(State *s) { clear_duties.insert(s); };
        // the main solver to be cleared
        std::unordered_set<MainSolver *> clear_duties_mainsolver;
        // add into clear duty, waiting for clear in the end.
        inline void clear_defer(MainSolver *s) { clear_duties_mainsolver.insert(s); };

    private:
        /**
         * @section Basic methods
         *
         */
        /**
         * @brief Create a O with given state object as its level 0.
         *
         * @param s
         * @return OSequence&
         */
        OSequence *createOWith(State *s);

        int pickStateLastIndex = -1;
        /**
         * @brief iteratively pick state from the given sequence.
         * @return State*
         * 	then this state will be marked as visited at this round.
         * @return nullptr every state has been tried.
         *  then all the markers will be cleaned. Therefore, this method can be reused for another round.
         * @post Anytime it ends the iterative round earlier(before nullptr), it must has found the counter example in trySAT. Therefore, there will be no picking in the future.
         */
        State *pickState();

        /**
         * @brief Get the partial state with assignment s. This is used in forward CAR.
         *
         * @param s
         * @param prior_state
         * @return State*
         */
        State *get_partial_state(Assignment &s, const State *prior_state);

        /**
         * @brief Get the solution from SAT solver.
         *
         * @return State* : the state representing the solution, which is to be added to the U sequence.
         */
        State *getModel();
        State *getModel(State *);

        /**
         * @brief Update U sequence, and push into cex vector
         *
         * @return whether it is already seen at a lower level.
         */
        bool updateU(State *, State *prior_in_trail);

        /**
         * @brief Update O sequence
         *
         * @param O
         * @param dst_level
         * @param Otmp
         * @param safe_reported
         */
        void updateO(OSequence *O, int dst_level, OFrame &Otmp, bool &safe_reported);

        /**
         * @brief SAT_Assume(assum, clauses)
         *
         * @return true
         * @return false
         */
        bool satAssume(OSequence *O, State *, int, OFrame &Otmp, bool &safe_reported);

        /**
         * @brief Interface for cleaning.
         *
         */
        void clean();

        /**
         * @brief When Os finds an invariant, clear Os and erase s from Os.
         *
         * @param s
         * @param Os
         */
        void clearO(State *s, OSequence *Os);

        /**
         * @brief first create a clause with target uc and flag of the corresponding O[dst], then add this clause to main solver or start solver, depending on the level.
         *
         * @param uc
         * @param O
         * @param dst_level
         */
        void addUCtoSolver(Cube &uc, OSequence *O, int dst_level, OFrame &Otmp);

        /**
         * @brief init special sequences: Uf, Ub, Oi, Onp
         */
        bool initSequence(bool &res);

        /**
         * @brief Use start solver to get a state. Usually in ~p.
         *
         * @param start_solver
         * @return State*
         */
        State *enumerateStartStates(StartSolver *start_solver);

        /**
         * @brief Check SAT_ASSUME(I, ~P).
         *
         * to generate uc from init
         * ~p should be in ~uc
         * use uc to initialize O[0] is suitable.
         *
         */
        bool immediateCheck(State *from, int target, bool &res, OFrame &O0);

        /**
         * @brief Check whether this in O[0] is actually a CEX.
         *
         * @param from
         * @param O
         * @return true
         * @return false
         */
        bool lastCheck(State *from, OSequence *O);

        /**
         * @brief Whether this state is blocked in this O frame, namely O[frame_level] (or Otmp, if frame_level == O.size)
         *
         * @param s
         * @param frame_level
         * @param O
         * @param Otmp
         * @return true
         * @return false
         */
        bool blockedIn(State *s, const int frame_level, OSequence *O, OFrame &Otmp);

        /**
         * @brief Use `blocked in` to iterate, from min to max, to find the minimal level where this state is not blocked.
         *
         * @param s
         * @param min
         * @param max
         * @param O
         * @param Otmp
         * @return int
         */
        int minNOTBlocked(State *s, const int min, const int max, OSequence *O, OFrame &Otmp);

    public:
        /**
         * @brief Helper function to print SAT query.
         *
         * @param solver
         * @param s
         * @param o
         * @param level
         * @param res
         * @param out
         */
        void print_sat_query(MainSolver *solver, State *s, OSequence *o, int level, bool res, std::ostream &out);
        // count of tried before

    private:
        /**
         * invariant section
         *
         */

        /**
         * @brief Check whether there exists an invariant in this sequence. Which means the center state of this OSequence is not reachable, and should be later blocked
         *
         * @param O the sequence to be checked
         * @return true : invariant is found
         * @return false : no invariant is found
         */
        bool InvFound(OSequence *O);

        bool InvFoundAt(OSequence &O, int check_level, int minimal_update_level, InvSolver *inv_solver);
    };

    namespace inv
    {
        void InvAddOR(OFrame &frame, int level, InvSolver *inv_solver_);

        void InvAddAND(OFrame &frame, int level, InvSolver *inv_solver_);

        void InvRemoveAND(InvSolver *inv_solver_, int level);

        bool InvFoundAt(OSequence &F_, const int frame_level, InvSolver *inv_solver_, int minimal_update_level_);
    }
}

#endif
