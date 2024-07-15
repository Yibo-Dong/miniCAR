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
    Author: Jianwen Li, Yibo Dong
    Init Date: September 8, 2017
    Update Date: Feb 12, 2023
    Interface for the checker class
*/

#ifndef carChecker_H
#define carChecker_H

#include "data_structure.h"
#include "invsolver.h"
#include "startsolver.h"
#include "mainsolver.h"
#include "newpartialsolver.h"
#include "model.h"
#include <assert.h>
#include "utility.h"
#include "statistics.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <fstream>
#include <map>

namespace car
{
    extern Statistics CARStats; // defined in main.cpp
    extern bool verbose_;    // defined in main.cpp
    extern Model *model;     // defined in main.cpp
    extern int storage_id;
    class Checker;
    extern Checker ch;

    class Checker
    {
    public:
        /**
         * @brief Construct a new Checker object
         *
         * @param model aiger model to be checked
         * @param out
         * @param forward the search direction
         * @param evidence whether to print the CEX if unsafe
         * @param index_to_check the index of property to check. At present, only one property is allowed

         */
        Checker(Model *model, std::ostream &out, bool forward = true, bool evidence = false, int index_to_check = 0, int convMode = -1, int convParam = 0, bool enable_rotate = false, int inter_cnt=0, bool inv_incomplete = false, bool uc_raw = false, int impMethod = 0, int LOStrategy = 0, int ConvAmount = 0, bool subStat = false);

        Checker(int time_limit, Model *model, std::ostream &out, bool forward = true, bool evidence = false, int index_to_check = 0, int convMode = -1, int convParam = 0, bool enable_rotate = false, int inter_cnt=0, bool inv_incomplete = false, bool uc_raw = false, int impMethod = 0, int LOStrategy = 0, int ConvAmount = 0, bool subStat = false);

        Checker(int time_limit, Checker* last_chker, int rememOption, Model *model, std::ostream &out, bool forward = true, bool evidence = false, int index_to_check = 0, int convMode = -1, int convParam = 0, bool enable_rotate = false, int inter_cnt=0, bool inv_incomplete = false, bool uc_raw = false, int impMethod = 0, int LOStrategy = 0, int ConvAmount = 0, bool subStat = false);

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
         * @brief Searching procedure of bi-car. Enumerate states in U and use the Osequence as the guide.
         *
         * @param U the U sequence to enumerate.
         * @param O the O sequence with a singleton (a state that is picked from the U sequence) making up its level 0.
         * @param forward whether the Transformation relationship should be used forward or backward
         * @return true : successfully reached O[0]. which means a cex is found.
         * @return false : all states in present U has been checked. No cex is found.
         */
        bool trySAT(Usequence &U, Osequence *O, bool forward, bool &safe_reported);

    public:
        /**
         * @section Preprocessing technique
         *
         */

        bool restart_enabled = false;
        bool ppstoped = false;
        int time_limit_to_restart = -1;
        bool importO = false;
        enum rememberEnum{
            remem_None = 0,
            remem_O0 = 1,
            remem_short = 2,
            remem_Ok = 3,
            remem_Uk = 4,
        };
        Checker* last_chker;
        int rememOption = 0;

        int inter_cnt = 0;
        bool rotate_enabled = false;
        bool inv_incomplete =false;
        bool uc_no_sort = false;
        int imply_decision = -1;

        enum ImpHowEnum{
            Imp_Manual = 0,
            Imp_Solver = 1,
            Imp_Sample = 2,
            Imp_Sort = 3,
            Imp_Exp = 4,
            Imp_Thresh = 5,
            Imp_MOM = 6,
            Imp_Fresh = 7,
        };

        int impMethod = 0;

        // level -> id -> freshIndex
        std::unordered_map<int, std::unordered_map <int, int>> impFreshRecord;
         
        clock_high sat_timer;

        enum LiteralOrderingEnum{
            LO_Classic = 0,     // Intersection + Rotation for (1), Reverse it for (2), Random for the subsequent ones.
            LO_Rand = 1,        // Just do random test.
            LO_Fixpoint = 2,    // Rotation for (1), {not_used, used} for later
            LO_Exp = 3,         // Reuse the phase-saving strategy to implement a biased-SAT instead.
        };
        int LOStrategy = 0;
        int ConvAmount = 0;
        bool subStat = false;

        inline void set_inter_cnt(int cnt) { inter_cnt = cnt; }
        inline int get_inter_cnt() { return inter_cnt; }
        inline void set_rotate(bool rtt) { rotate_enabled = rtt; }
        inline int get_rotate() { return rotate_enabled; }

        /**
         * @brief extract the knowledge from prior preprocessing phases heavily
         *
         */
        void cook();

        /**
         * @brief extract the knowledge from prior preprocessing phases heavily
         *
         */
        void cook_light();

    public:
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
        int convMode=-1;
        int convParam=0;

        // for better manual method
        // maintain the index to visit.
        std::vector<std::set<std::pair<int,int>>> uc_len_indexes;
        void insert_to_uc_index(Cube &uc,int index, int level);

    public:
        Model *model_;
        bool evidence_;
        bool backward_first = true;
        int bad_;
        StartSolver *bi_start_solver;

        // mapping from state to its corresponding o sequence.
        std::unordered_map<const State *, Osequence *> SO_map;
        // the main solver shared.
        MainSolver *bi_main_solver;
        // the partial solver shared.
        PartialSolver *bi_partial_solver;
        // count of blocked states.
        int blocked_count = 0;

        // the map from O sequence to its minimal_level
        std::unordered_map<const Osequence *, int> fresh_levels;
        Usequence Uf, Ub; // Uf[0] is not explicitly constructed
        Osequence Onp, OI;
        // used in picking state randomly
        std::vector<std::pair<State *, int>> Uset;

    public:
        /**
         * @section Counter Example
         *
         */

        inline Usequence &whichU()
        {
            return backward_first ? Ub : Uf;
        }

        inline Usequence &otherU()
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
         * @return Osequence&
         */
        Osequence *createOWith(State *s);

        int pickStateLastIndex = -1;
        /**
         * @brief iteratively pick state from the given sequence.
         * @return State*
         * 	then this state will be marked as visited at this round.
         * @return nullptr every state has been tried.
         *  then all the markers will be cleaned. Therefore, this method can be reused for another round.
         * @post Anytime it ends the iterative round earlier(before nullptr), it must has found the counter example in trySAT. Therefore, there will be no picking in the future.
         */
        State *pickState(Usequence &);

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
        State *getModel(MainSolver *);
        State *getModel(MainSolver *, State *);

        /**
         * @brief Update U sequence, and push into cex vector
         *
         * @return whether it is already seen at a lower level.
         */
        bool updateU(Usequence &, State *, State *prior_in_trail);

        /**
         * @brief Update O sequence
         *
         * @param O
         * @param dst_level
         * @param Otmp
         * @param safe_reported
         */
        void updateO(Osequence *O, int dst_level, Frame &Otmp, bool &safe_reported);

        /**
         * @brief SAT_Assume(assum, clauses)
         *
         * @return true
         * @return false
         */
        bool satAssume(MainSolver *, Osequence *O, State *, int, Frame &Otmp, bool &safe_reported);

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
        void clearO(State *s, Osequence *Os);

        /**
         * @brief first create a clause with target uc and flag of the corresponding O[dst], then add this clause to main solver or start solver, depending on the level.
         *
         * @param uc
         * @param O
         * @param dst_level
         */
        void addUCtoSolver(Cube &uc, Osequence *O, int dst_level, Frame &Otmp);

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
        bool immediateCheck(State *from, int target, bool &res, Frame &O0);

        /**
         * @brief Check whether this in O[0] is actually a CEX.
         *
         * @param from
         * @param O
         * @return true
         * @return false
         */
        bool lastCheck(State *from, Osequence *O);

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
        bool blockedIn(State *s, const int frame_level, Osequence *O, Frame &Otmp);

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
        int minNOTBlocked(State *s, const int min, const int max, Osequence *O, Frame &Otmp);

    public:
        void print_flags(std::ostream &out);
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
        void print_sat_query(MainSolver *solver, State *s, Osequence *o, int level, bool res, std::ostream &out);
        std::vector<int> spliter;
        std::set<int> blocked_ids;
        std::vector<int> blocked_counter_array;
        // count of tried before
        int direct_blocked_counter = 0;

    private:
        /**
         * @section rotation technique
         *
         */
        // corresponding to O[i]
        std::vector<Cube> rotates;
        // corresponding to Otmp
        Cube rotate;

    private:
        /**
         * @section use failed states to score
         *
         */
        std::vector<std::unordered_map<int, int>> score_dicts;
        std::unordered_map<int, int> score_dict;
        std::unordered_map<int, int> decayStep;
        std::unordered_map<int, int> decayCounter;

    private:
        /**
         * invariant section
         *
         */

        /**
         * @brief Check whether there exists an invariant in this sequence. Which means the center state of this Osequence is not reachable, and should be later blocked
         *
         * @param O the sequence to be checked
         * @return true : invariant is found
         * @return false : no invariant is found
         */
        bool InvFound(Osequence *O);

        bool InvFoundAt(Osequence &O, int check_level, int minimal_update_level, InvSolver *inv_solver);
    };

    namespace inv
    {
        void InvAddOR(Frame &frame, int level, InvSolver *inv_solver_);

        void InvAddAND(Frame &frame, int level, InvSolver *inv_solver_);

        void InvRemoveAND(InvSolver *inv_solver_, int level);

        bool InvFoundAt(Fsequence &F_, const int frame_level, InvSolver *inv_solver_, int minimal_update_level_);
    }
}

#endif
