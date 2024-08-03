#include "carChecker.h"
#include "implysolver.h"
#include "newpartialsolver.h"
#include "statistics.h"
#include "definition.h"
#include <vector>
#include <iostream>
#include <stack>
#include <functional>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <tuple>
#include <chrono>
using namespace std;
using namespace std::chrono;
using clock_high = time_point<steady_clock>;

namespace car
{
    Checker::Checker(Problem *model, const OPTIONS &opt, std::ostream &out, Checker *last_chker) : 
    out(out), model_(model), rotate_enabled(opt.enable_rotate), inter_cnt(opt.inter_cnt), inv_incomplete(opt.inv_incomplete), uc_no_sort(opt.raw_uc), impMethod(opt.impMethod), time_limit_to_restart(opt.time_limit_to_restart), rememOption(opt.rememOption), LOStrategy(opt.LOStrategy), convAmount(opt.convAmount), convParam(opt.convParam), convMode(opt.convMode), subStat(opt.subStat), partial(opt.partial), last_chker(last_chker), fresh_levels(0), backwardCAR(!opt.forward), bad_(model->output(0)), inv_solver(nullptr), multi_solver(opt.multi_solver),propMode(opt.propMode), propParam(opt.propParam)
    {
        if(!multi_solver)
        {
            main_solver = new MainSolver(model, opt.forward, get_rotate(), uc_no_sort);
            // also initialize the prop solver.
            prop_solver = new MainSolver(model, opt.forward, get_rotate(), uc_no_sort);
        }    

        if (!backwardCAR)
        {
            // only forward needs these
            partial_solver = new PartialSolver(model);
            // TODO: extend to multi-properties.
            start_solver = new StartSolver(model, bad_, true);
        }
        else
        {
            partial_solver = nullptr;
            start_solver = nullptr;
        }
        rotates.clear();
        rotate.clear();
        if(opt.time_limit_to_restart > 0)
            restart_enabled = true;
    }

    Checker::~Checker()
    {
        clean();
        if (main_solver)
        {
            delete main_solver;
            main_solver = nullptr;
        }
        if(!main_solvers.empty())
        {
            for(auto& slv: main_solvers)
            {
                delete slv;
                slv = nullptr;
            }
        }
        if (start_solver)
        {
            delete start_solver;
            start_solver = nullptr;
        }
        if (partial_solver)
        {
            delete partial_solver;
            partial_solver = nullptr;
        }
        if (prop_solver)
        {
            delete prop_solver;
            prop_solver = nullptr;
        }
    }

    bool Checker::check()
    {
        bool res;
        if (trivialCheck(res))
        {
            return res ? true : false;
        }

        res = car();
        if (res && ppstoped)
        {
            // should try the next strategy.
            // go out, build another checker.
            return true;
        }

        do
        {
            if (res)
            {
                out << "0" << endl;
                out << "b0" << endl;
                out << "." << endl;
                break;
            }

            print_evidence();
        } while (false);

        return res;
    }

    MainSolver *Checker::getMainSolver(int level)
    {
        if(!multi_solver)
            return main_solver;
        else
        {
            CARStats.count_1_begin();
            if(level >= main_solvers.size())
            {
                int needs = level - main_solvers.size() + 1;
                while(needs > 0)
                {
                    auto *new_slv = new MainSolver(model_, !backwardCAR, get_rotate(), uc_no_sort);
                    main_solvers.push_back(new_slv);                        
                    --needs;
                }
            }
            CARStats.count_1_end(1);
            if(level <=0)
                return main_solvers.at(0);
            else
                return main_solvers.at(level);
        }
    }

    bool Checker::car()
    {
        bool res = false;
        if (initSequence(res))
            return res;
        if (restart_enabled)
        {
            // record the time that we entered.
            sat_timer = steady_clock::now();
        }

        /**
         * As to backward search,
         * @post Ob, Ub are both initialized.
         * uc of UNSAT(I, ~P) is inserted into O[0].
         * O[0] is inserted into MainSolver.
         */

        while (true)
        {   
            if (trySAT(res))
            {
                // if this PP phase ends, but has not reached the whole end.
                if (ppstoped)
                {
                    return true;
                    // continue;
                }
                return res;
            }

            if (!inv_incomplete)
            {
                if (InvFound())
                    return true;
            }

        }

        // dead code. Should not reach
        return false;
    }

    bool Checker::trivialCheck(bool &res)
    {
        const Problem *model = model_;
        // FIXME: fix for multiple properties.
        if (bad_ == model->true_id())
        {
            out << "1" << endl;
            out << "b"
                << "0" << endl;
            {
                // print init state
                // FIXME: fix for latch inits.
                for(int i = 0; i < model->num_latches(); ++i)
                    out << "0";
                out << endl;
                // print an arbitary input vector
                for(int j = 0; j < model->num_inputs(); j++)
                    out << "0";
                out << endl;
            }
            out << "." << endl;
            res = true;
            return true;
        }
        else if (bad_ == model->false_id())
        {
            out << "0" << endl;
            out << "b" << endl;
            out << "." << endl;
            res = false;
            return true;
        }
        return false;
    }

    bool Checker::trySAT(bool &safe_reported)
    {
        // NOTE: can eliminate initialization.
        auto& O = whichO();
        safe_reported = false;
        refreshOtmp();
        CARStats.count_enter_new_ronud();
        /**
         * this procedure is like the old car procedure, but the OSequence is not bound to be OI or Onegp.
         * @param missionary the state to be checked
         */
        while (State *missionary = pickState())
        {
            /**
             * build a stack <state, depth, target_level>
             */
            stack<Obligation> stk;
            stk.push(Obligation(missionary, 0, OSize() - 1));
            while (!stk.empty())
            {
                CARStats.count_enter_new_try_by();
                State *s;
                int dst, depth;
                std::tie(s, depth, dst) = stk.top();
                if (blockedIn(s, dst + 1))
                {
                    // TODO: memorize state's blocking status. Since we do not remove UC, once blocked, forever blocked.
                    stk.pop();
                    CARStats.count_tried_before();

                    int new_level = minNOTBlocked(s, dst + 2, OSize() - 1);
                    if (new_level <= OSize())
                    {
                        stk.push(Obligation(s, depth, new_level - 1));
                    }
                    continue;
                }

                if (restart_enabled)
                {
                    auto now = steady_clock::now();
                    duration_high elapsed = now - sat_timer;
                    double time_delay = elapsed.count();
                    if (time_delay > time_limit_to_restart * 1000)
                    {
                        ppstoped = true;
                        return true;
                    }
                }

                if (satAssume(s, dst, safe_reported))
                {
                    if (dst == -1)
                    {
                        return true;
                    }
                    State *tprime = getModel(dst);
                    
                    updateU(tprime, s);

                    // NOTE: why here calculate minNOTBLOCKED, rather than next time when pop?
                    int new_level = minNOTBlocked(tprime, 0, dst - 1);               
                    if (new_level <= dst) // if even not one step further, should not try it
                    {
                        stk.push(Obligation(tprime, depth + 1, new_level - 1));
                    }
                }
                else
                {
                    stk.pop();
                    if (safe_reported)
                        return true;

                    int new_level = minNOTBlocked(s, dst + 2, OSize() - 1);
                    if (new_level <= OSize())
                    {
                        stk.push(Obligation(s, depth, new_level - 1));
                    }
                }
            }
        }

        // this is same with extend_F_sequence in old CAR.
        O.push_back(Otmp);

        if (rotate_enabled)
            rotates.push_back(rotate);
        getMainSolver(OSize() - 1)->add_new_frame_M(Otmp, OSize() - 1);
        if (propMode>0)
            prop_solver->add_new_frame_P(Otmp, OSize() - 1);
        return false;
    }

    /**
     * @brief
     *
     * @pre level less than check_level has been checked before
     * @param O
     * @param check_level
     * @return true : invariant found at target level
     * @return false : invariant does not exists at target level
     */
    bool Checker::InvFoundAt(int check_level, int minimal_update_level)
    {
        OSequence &O = whichO();
        // a portion of `InvFound()`
        if (check_level < minimal_update_level)
        {
            inv_solver->add_constraint_or(O[check_level],check_level);
            return false;
        }
        inv_solver->add_constraint_and(O[check_level], check_level);

        bool res = !inv_solver->solve_with_assumption();

        inv_solver->release_constraint_and(check_level);
        inv_solver->add_constraint_or(O[check_level],check_level);
        return res;
    }

    bool Checker::InvFound()
    {
        OSequence &O = whichO();
        bool res = false;
        // FIXME: Should we reuse inv_solver instead of recreating?
        inv_solver = new InvSolver(model_);
        // FIXED: shall we start from 0 or 1?
        // 0, to add O[0] into solver.

        /**
         * @brief About minimal update level:
         * 		It is related to specific O sequence.
         * 		Each time invariant check is done, minimal update level is reset to the highest level.
         * 		Each time a modification of low-level frame will possibly modify this minimal level.
         *
         */
        for(int i = 0; i < OSize(); ++i)
        {
            if (InvFoundAt(i, fresh_levels))
            {
                while (OSize() > i)
                    O.pop_back();
                res = true;
                // already found invariant.
                fresh_levels = -1;
                break;
            }
        }
        // NOTE: not OSize()-1. because that level is also checked.
        fresh_levels = OSize();
        delete (inv_solver);
        inv_solver = nullptr;
        return res;
    }

    // NOTE: when we try to pick a guide state, it is not used as part of the assumption but the target. Therefore, uc is not generated according to it, and therefore will not be updated.
    // Luckily, we do not use Onp or OI to guide. Therefore, here we have nothing to do with start_solver.
    State *Checker::pickState()
    {
        auto& U = whichU();
        if(pickStateLastIndex <= 0)
        {
            // one round is end.
            pickStateLastIndex = U.size();
            return nullptr;
        }
        --pickStateLastIndex;
        State* s = U[pickStateLastIndex];
        if(s->_isnegp)
        {
            // FIXME: check for forward CAR, where we need to enumerate by querying the start solver.
            // TODO: check for efficiency.
            State *new_s = enumerateStartStates();
            if(new_s)
            {
                updateU(new_s,nullptr);
                if(start_solver)
                    start_solver->reset();
                return new_s;
            }
            else
            {
                // negp should be at position 0.
                // so when we try to fetch a state, but get negp, it is also the end of one round.
                pickStateLastIndex = U.size();

                if(start_solver)
                    start_solver->reset();
                return nullptr;
            }
            return nullptr;
        }
        return s;
    }


    /**
     * @brief About the principle here, see newpartialsolver.h
     *
     * @param s
     * @param prior_state
     */
    State *Checker::get_partial_state(Assignment &s, const State *prior_state)
    {
        Assignment next_assumptions = s;
        Assignment old_inputs(s.begin(), s.begin() + model_->num_inputs());

        if (prior_state)
        // it is not start state
        {
            // negate the s' as the cls, put into the solver.
            Assignment cls = negateCube(prior_state->get_latches());
            partial_solver->new_flag();
            partial_solver->add_clause_with_flag(cls);

            // the assumption being t's input and latches.
            next_assumptions.push_back(partial_solver->get_flag());
            partial_solver->set_assumption(next_assumptions);

            // Therefore, it ought to be unsat. If so, we can get the partial state (from the uc)
            bool res = partial_solver->solve_assumption();
            // uc actually stores the partial state.
            s = partial_solver->get_conflict();

            if (res || s.empty())
            {
                // all states can reach this state! It means the counter example is found. Sure the initial state can reach this state.
                // set the next state to be the initial state.
                assert(Ub.size());
                State *init = Ub[0];
                s = init->get_latches();
            }

            // block this clause.
            int flag = partial_solver->get_flag();
            partial_solver->add_clause(flag);
        }
        else
        // for initial states, there is no such "prior state", only "bad"
        {
            next_assumptions.push_back(-bad_);
            partial_solver->set_assumption(next_assumptions);
            bool res = partial_solver->solve_assumption();
            s = partial_solver->get_conflict_no_bad(-bad_);
            // if one-step reachable, should already been found in immediate_satisfiable.
            assert(!s.empty());
        }

        // FIXME: is this `inputs` really useful?
        Assignment inputs, latches;
        for (auto &it : s)
        {
            if (abs(it) <= model_->num_inputs())
                inputs.push_back(it);
            else
                latches.push_back(it);
        }
        if (inputs.empty()) // use this inputs as the inputs.
        {
            inputs = old_inputs;
        }

        // this state is not a full state.
        State *pstate = new State(inputs, latches);
        clear_defer(pstate);
        return pstate;
    }

    /**
     * @brief print the evidence. reuse backwardCAR, which exactly reveals present searching direciton.
     * @pre counterexmple is met already.
     * @param 
     */
    void Checker::print_evidence() const
    {
        // Print Backward chain first.
        bool latch_printed = false;

        out << "1" << endl;
        out << "b" << 0 << endl;
        if (counter_start_b)
        // backward_chain is not empty
        {
            State *to_print = counter_start_b;
            std::stack<std::string> helper;
            // if this is just a portion of the chain, there may not be last_inputs.
            // if(!to_print->last_inputs().empty())
            // 	helper.push(to_print->last_inputs());
            while (to_print)
            {
                State *next = prior_in_trail_b.find(to_print) != prior_in_trail_b.end() ? prior_in_trail_b.at(to_print) : nullptr;
                if (next)
                    helper.push(to_print->get_inputs_str());
                else
                    helper.push(to_print->get_latches_str());
                to_print = next;
            }
            while (!helper.empty())
            {
                out << helper.top() << endl;
                helper.pop();
            }
            latch_printed = true;
        }
        // then print forward chain.
        if (counter_start_f)
        {
            State *to_print = counter_start_f;
            if (!latch_printed)
                out << to_print->get_latches_str() << endl;
            while (to_print)
            {
                State *next = prior_in_trail_f.find(to_print) != prior_in_trail_f.end() ? prior_in_trail_f.at(to_print) : nullptr;
                out << to_print->get_inputs_str() << endl;
                to_print = next;
            }
        }
        out << "." << endl;
    }

    bool Checker::updateU(State *s, State *prior_state_in_trail)
    {
        // Counter Example Issue.
        // Every time we insert a new state into U sequence, it should be updated.
        // FIXME: without considering bi-direction, one state will only have one prior. Shall we move this into "state"?
        {
            whichPrior()[s] = prior_state_in_trail;
        }
        whichU().push_back(s);

        return true;
    }

    vector<Cube> reorderAssum(const vector<Cube>& inter, const Cube &rres, const Cube &rtmp)
    {
        vector<Cube> pref;
        pref.reserve(inter.size() + 2);
        pref = inter;
        pref.push_back(std::move(rres));
        pref.push_back(std::move(rtmp));
        return pref;
    }

    bool Checker::satAssume(State *s, int level, bool& safe_reported)
    {
        vector<Cube> inter;
        Cube rres, rtmp;
        vector<Cube> ucs; // this is used for subsumption test, will cause low-efficiency, do not turn it on in practical scenarios.

        MainSolver* ms = getMainSolver(level);
        bool res = false;
        if (level == -1)
        {
            // NOTE: in backward CAR, here needs further check.
            CARStats.count_main_solver_original_time_start();
            res = finalCheck(s);
            if(res)
            {
                // if sat, no uc, we just update it here.
                CARStats.count_main_solver_original_time_end(res,0);
            }

        }
        else
        {
            do
            {
                // intersection:
                if (get_inter_cnt())
                {
                    // NOTE: the meaning of inter_cnt is not consistent as to ni > 1
                    const OFrame &frame = whichFrame(level + 1);
 
                    // this index is which iCube we want to create.
                    int index = 1;

                    int uc_index = frame.size();
                    // in Conv, not every UC is original UC.
                    int record_bits = conv_record[level + 1];

                    while (index <= get_inter_cnt() && frame.size() >= size_t(index))
                    {
                        Cube inter_next;

                        uc_index -= 1;
                        if (convMode >= 0)
                        {
                            // if we generated more than UC, we shall just use the original one.
                            uc_index -= (record_bits & (1 << (index - 1))) ? 1 : 0;
                            if (uc_index < 0)
                                break;
                            // cerr<<"uc index is "<< uc_index<<"/"<<frame.size()<<endl;
                            // cerr<<"index is "<<index<<", record bit is "<< ((record_bits & (1<<(index-1))) ? 1 : 0 )<<endl;
                            // cerr<<"picked uc at "<<uc_index<<endl;
                        }
                        const Cube &last_uc = frame[uc_index];

                        inter_next = s->intersect(last_uc);

                        // otherwise, do not do this!
                        inter.push_back(inter_next);
                        ++index;
                    }
                    if (inter.empty())
                        inter.push_back({});
                }
            } while (0);

            do
            {
                if(get_rotate())
                {
                // rotates[i]
                if (!rotate_enabled)
                    break;
                Cube &rcu = size_t(level + 1) < rotates.size() ? rotates[level + 1] : rotate;
                if (rcu.empty())
                {
                    rres = {};
                    rtmp = s->get_latches();
                    // TODO: try this to be inter?
                    break;
                }

                // a full state
                if (s->latch_size() == model_->num_latches())
                {
                    // calculate intersection and put the others behind.
                    for(size_t i = 0; i < rcu.size(); ++i)
                    {
                        if (s->latch_at(abs(rcu[i]) - model_->num_inputs() - 1) == rcu[i])
                            rres.push_back(rcu[i]);
                        else
                            rtmp.push_back(-rcu[i]);
                    }
                }
                else
                {
                    // a partial state
                    // TODO: test me!
                    int i = 0, j = 0;
                    while (i < s->latch_size() && size_t(j) < rcu.size())
                    {
                        if (rcu[j] == s->latch_at(i))
                        {
                            rres.push_back(rcu[j]);
                            i++;
                            j++;
                        }
                        else
                        {
                            rtmp.push_back(s->latch_at(i));
                            if (rcu[j] < s->latch_at(i))
                                j++;
                            else
                                i++;
                        }
                    }                    
                }
                }
            } while (0);


            vector<Cube> pref = reorderAssum(inter, rres, rtmp);
            

            switch (LOStrategy)
            {
            case 1: // random
            {
                std::random_device rd;
                std::mt19937 gen(rd());
                vector<int> shuff_helper = s->get_latches();
                for(int i = shuff_helper.size() - 1; i > 0; --i) {
                    // starting from 1, because 0 is the flag.
                    std::uniform_int_distribution<int> dis(0, i);
                    int j = dis(gen);
                    std::swap(shuff_helper[i], shuff_helper[j]);
                }
                ms->set_assumption(s, level,  {shuff_helper});
                break;
            }

            case 3: // expectation test
            {
                // FIXME: implement this part to be a real decision procedure.
                std::vector<int> exp = s->get_latches(); // tmp only.
                ms->set_expectation(exp);
                ms->set_assumption(s, level,  pref);
                break;
            }

            case 2:
            case 0: // traditional one, no need to shuffel
            default:
            {
                ms->set_assumption(s, level,  pref);
                // do nothing
                break;
            }
            }
            
            CARStats.count_main_solver_original_time_start();
            res = ms->solve_assumption();
            if(res)
            {
                // if sat, no uc, we just update it here.
                CARStats.count_main_solver_original_time_end(res,0);
            }
        }
        if(!res)
        {
            // update the UC.
            Cube uc = ms->get_conflict();
            CARStats.count_main_solver_original_time_end(res,uc.size());
            
            if (uc.empty())
            {
                cerr<<"uc is empty@1!"<<endl;
                safe_reported = true;
                cerr<<"error level is "<<level<<endl;
                cerr<<"corresponding O frame is "<<endl;
                for(auto &cls: whichFrame(level))
                {
                    for(auto &lit: cls)
                    {
                        cerr<<lit<<" ";
                    }
                    cerr<<endl;
                }
                cerr<<"end printing the O frame"<<endl;
            }
            if(subStat)
            {
                ucs.push_back(uc);
                CARStats.recordUC(false);
            }

            // note: should do this before propagation.
            addUCtoSolver(uc, level + 1);

            if (propMode > 0 && level > 0 && level + 1 < OSize())
                propagate(uc, level);
            
            if (uc.empty())
            {
                cerr<<"uc is empty@2!"<<endl;
                safe_reported = true;
            }


        }

        if (get_rotate() && !res)
        {
            // update rotate
            Cube &rcu = size_t(level + 1) < rotates.size() ? rotates[level + 1] : rotate;
            rcu = rres;
            rcu.insert(rcu.end(), rtmp.begin(), rtmp.end());
        }

        if(convMode >= 0)
        {
            // <others> <skipped> <previousUC>

            // whether to calculate another UC.
            bool trigger = false;
            switch(convMode)
            {
                case ConvModeAlways:
                {    
                    trigger = true;
                    // cerr<<"conv mode = Always"<<endl;
                    break;
                }
                case ConvModeHigh:
                {    
                    if(convParam <=0)
                    {
                        trigger = false;
                        break;
                    }
                    if (convParam == 1)
                    {
                        trigger = true;
                        break;
                    }
                    if(level > 0 && level > OSize() - (OSize() / convParam))
                    {
                        trigger = true;
                    }
                    else{
                        trigger = false;
                    }
                    // cerr<<"conv mode = High, res = "<<trigger<<" level="<<level<<"/"<<OSize()<<endl;
                    break;
                }
                case ConvModeLow:
                {
                    if(convParam <=0)
                    {
                        // static bool printed = false;
                        // if(!printed)
                        // {
                        //     printed = true;
                        //     cerr<<"always no"<<endl;
                        // }
                        trigger = false;
                        break;
                    }
                    if (convParam == 1)
                    {
                        // static bool printed = false;
                        // if(!printed)
                        // {
                        //     printed = true;
                        //     cerr<<"always yes"<<endl;
                        // }
                        trigger = true;
                        break;
                    }

                    if (level < 0 || level <  1 + (OSize() / convParam))
                    {
                        trigger = true;
                    }
                    else{
                        trigger = false;
                    }
                    // cerr<<"conv mode = Low, res = "<<trigger<<" level="<<level<<"/"<<OSize()<<endl;
                    break;
                }
                case ConvModeRand:
                {
                    assert(convParam != 0);
                    static mt19937 mt_rand(1);
                                
                    trigger = ((mt_rand() % convParam) == 0) ? true : false;
                    
                    // cerr<<"conv mode = Random, res = "<<trigger<<endl;
                    break;
                }
                case ConvModeStuck:
                {
                    // TODO: fill in here.
                    break;
                }

                default:
                    trigger = true;
            }

            if (!res && trigger)
            {
                int amount = 0;
                while(++amount <= convAmount)
                {
                // retrieve another bit.
                conv_record[level + 1] <<= 1;

                // mark as inserted
                conv_record[level + 1] += 1;


                CARStats.count_main_solver_convergence_time_start();
                // get another conflict!
                auto nextuc = ms->get_conflict_another(LOStrategy, amount);
                
                if(subStat)
                {
                    // note, this is slow..
                    bool sub = false;
                    for(auto& uc: ucs)
                    {
                        if(imply_heavy(nextuc, uc))
                        {
                            CARStats.recordUC(true);
                            sub = true;
                            break;
                        }
                    }
                    if(!sub)
                    {
                        CARStats.recordUC(false);
                        ucs.push_back(nextuc);
                    }
                }

                CARStats.count_main_solver_convergence_time_end(nextuc.size());
                //TODO: analyse, whether imply or implied.

                addUCtoSolver(nextuc, level + 1);
                }
                // TODO: subsumption test is not cost-effective as we tested before, but we may use it to illustrate the diminishing effects.
                
            }
        }    

        return res;
    }

    /**
     * @brief Check for UC propagation to lower frames(in Backward CAR)
     * 
     * basis: at level l+1: !uc
     *     |= !uc /\ T -> (!uc')
     * <=> UNSAT(!uc /\ T /\ (uc'))
     */
    void Checker::propagate(Cube& uc, int level)
    {
        if(level <= 0 ||  level + 1 >= OSize())
            return;
        /**
         * @todo Play with it!
         * 1. Only short UCs?
         * 2. Only to particular Frames? (close ones, low-level ones, small-size ones)
         * 3. Trigger time: Only if the 2nd UC == 1st UC? (It seems to be the unique UC?)
         */

        switch (propMode)
        {
        
            case PropAlways:
            {
                // do it.
                break;
            }
            case PropShortCont:
            case PropShort:
            {
                if(uc.size() > propParam)
                    return;
                // do it.
                break;
            }

            case PropContinue:
            {
                break;
            }
            
            default: // fall through.
            case PropNone:
            {
                return;
            }
        }

        // should be right after an UNSAT call.
        CARStats.count_prop_begin();
        prop_solver->clear_assumption();
        
        vector<int> assumption; 
        // assumption ::= <Oflag of l+1>, <deactivation flags>, <activation flag>,<uc'>
        assumption.reserve(1 + prop_solver->PTFlag.size() + 1 + uc.size());

        // two approches:
        // 1): insert all the fresh clauses of current state here.
        // 2): everytime we update the main solver, we also update the prop_solver.
        // Choose 2) at present to avoid maintaining the indexes of 'fresh'.

        int Oflag = prop_solver->MFlagOf(level+1);
        assumption.push_back(Oflag);

        // deactivate all the previous ones.
        for(int flg: prop_solver->PTFlag)
        {
            assumption.push_back(flg);
        }
        int new_flag = prop_solver->getNewPTFlag();
        // activate the current one:
        assumption.push_back(-new_flag);

        // the UC(will be primed later in `set_assumption_primed`):
        assumption.insert(assumption.end(), uc.begin(), uc.end());

        prop_solver->set_assumption_primed(assumption); 

        // clause ::= !uc , current flag.
        prop_solver->add_clause_from_cube(uc,new_flag, true);  // !uc
        
        // check if !uc is inductive
        bool prop_res = prop_solver->solve_assumption();

        // do minimization of UC:
        if(!prop_res)
        {
            auto upcoming_uc = prop_solver -> get_uc();
            /**
             * The upcoming UC ::= <subset of uc> [<activation flag>] [<Oflag>]
             * Reason:
             * 1) O flag is the first in the assumption vector;
             *  either it appears ==> should be in the back
             *  or it does not appear
             * 2) The deactivation literals are pure literals, should not appear.
             * 3) The activation flag is likely to appear -- Otherwise, T /\ uc' is tautology.
            */

            bool strongInductive = false;// regardless of the frame.
            // assert: flag should only appear in the back, or not.

            // again, remember:
            // assumption ::= <Oflag of l+1>, <deactivation flags>, <activation flag>,<uc'>
            if(prop_solver->MLevelOf(upcoming_uc.back()) == NOT_M_FLAG)
            {
                // O flag does not participate.
                strongInductive = true;
            }
            else
            {
                upcoming_uc.pop_back();
                strongInductive = false;
            }

            if(upcoming_uc.back() == -new_flag)
            {
                upcoming_uc.pop_back();
            }
            
            if(upcoming_uc.size() < uc.size())
            {
                // map prop_uc -> 'previous' version.
                // note: not unique, just take the first one.
                model_->transform_uc_to_previous(upcoming_uc);              
                //TODO:change me back!
                uc = upcoming_uc;

            }

            if(strongInductive)
            {
                // propagate to all previous levels                
                for (int f = level; f >= 0; --f)
                {
                    addUCtoSolver(uc, f);
                }                
            }
            else{
                // propagate to only level l.
                addUCtoSolver(uc, level);
                // we could check whether it is relative to previous frames.
                if(propMode == PropContinue || propMode == PropShortCont)
                    propagate(uc,level - 1);
            }
        }

        CARStats.count_prop_end(!prop_res);
    }

    State *Checker::getModel(int level)
    {
        State *s = getMainSolver(level)->get_state();
        // NOTE: if it is the last state, it will be set last-inputs later.
        clear_defer(s);
        return s;
    }


    void Checker::clean()
    {
        for (State *duty : clear_duties)
        {
            if (duty)
            {
                delete duty;
                duty = nullptr;
            }
        }
    }

    void Checker::addUCtoSolver(Cube &uc, int dst_level_plus_one)
    {
        if (dst_level_plus_one < fresh_levels)
            fresh_levels = dst_level_plus_one;

        OFrame &frame = whichFrame(dst_level_plus_one);

        frame.push_back(uc);

        ImplySolver::add_uc(uc,dst_level_plus_one);
        

        if (dst_level_plus_one < OSize())
        {
            getMainSolver(dst_level_plus_one)->add_clause_from_cube_M(uc, dst_level_plus_one);
            if(propMode > 0)
            {
                // also update the propagation solver
                prop_solver->add_clause_from_cube_P(uc, dst_level_plus_one);
            }
        }
        else if (dst_level_plus_one == OSize())
        {
            if (!backwardCAR)
            {
                // FIXME: test me later.
                // Not always. Only if the start state is ~p.
                start_solver->add_clause_with_flag(uc);
            }
        }
    }

    bool Checker::initSequence(bool &res)
    {
        State *init = new State(model_->init());
        State *negp = new State(true);
        clear_defer(init);
        clear_defer(negp);
        OFrame O0; // a frame with only one uc.

        switch (rememOption)
        {
            case(remem_None):
            {
                if (immediateCheck(init, res, O0))
                    return true;

                // forward inits.
                if (!backwardCAR)
                {
                    // Uf[0] = ~p;

                    // NOTE: we do not explicitly construct Uf[0] data strcutre. Because there may be too many states.  We actually need concrete states to get its prime. Therefore, we keep an StartSolver here, and provide a method to enumerate states in Uf[0].
                    // construct an abstract state
                    updateU(negp,nullptr);
                    pickStateLastIndex = Uf.size();
                    
                    // O_I
                    // NOTE: O stores the UC, while we actually use negations.
                    OFrame frame;
                    auto& latches = init->get_latches();
                    frame.resize(latches.size());
                    std::transform(latches.begin(), latches.end(), frame.begin(), [](auto lit) {
                        return std::vector<int>{-lit};
                    });
                    OI = {frame};
                }
                // backward inits.
                if (backwardCAR)
                {
                    // Ub[0] = I;
                    updateU(init, nullptr);
                    pickStateLastIndex = Ub.size();

                    // Ob[0] = uc0.
                    // clauses will be added by immediate_satisfible.
                    // SAT_assume(init, ~p)
                    // uc from init
                    // ~p in ~uc
                    // use uc to initialize O[0] is suitable.
                    for(size_t index = 0; index< O0.size(); ++index)
                    {
                        ImplySolver::add_uc(O0[index],0);
                    }
                    Onp = OSequence({O0});
                }
                break;
            }

            case(remem_O0):
            {
                // If we remember all, it will be too many.
                // NOTE: sonly implement backward now.
                assert(backwardCAR);
                O0 = last_chker->Onp[0];
                
                if (backwardCAR)
                {
                    // Ub[0] = I;
                    updateU(init, nullptr);
                    pickStateLastIndex = Ub.size();

                    // Ob[0] = uc0.
                    // clauses will be added by immediate_satisfible.
                    // SAT_assume(init, ~p)
                    // uc from init
                    // ~p in ~uc
                    // use uc to initialize O[0] is suitable.

                    for(size_t index = 0; index< O0.size(); ++index)
                    {
                        ImplySolver::add_uc(O0[index],0);
                    }
                    Onp = OSequence({O0});
                }
                break;
            }
            case(remem_short):
            {
                // If we remember all, it will be too many.
                // NOTE: sonly implement backward now.
                assert(backwardCAR);
                O0 = last_chker->Onp[0];
                
                sort(O0.begin(),O0.end(),[](const vector<int>& a, const vector<int>& b){return a.size() < b.size();});
                
                // int len = 0;
                // for(auto &uc:O0)
                //     len+=uc.size();
                // cerr<<"before: sz = "<<O0.size()<<", avg = "<<len/O0.size()<<endl;

                static float portion = 2;
                O0.resize(O0.size() * (1-(1/portion)));
                portion++;

                // int after_len = 0;
                // for(auto &uc:O0)
                //     after_len+=uc.size();

                // cerr<<"after: sz = "<<O0.size()<<", avg = "<<after_len/O0.size()<<endl;

                if (backwardCAR)
                {
                    // Ub[0] = I;
                    updateU(init, nullptr);
                    pickStateLastIndex = Ub.size();

                    // Ob[0] = uc0.
                    // clauses will be added by immediate_satisfible.
                    // SAT_assume(init, ~p)
                    // uc from init
                    // ~p in ~uc
                    // use uc to initialize O[0] is suitable.

                    for(size_t index = 0; index< O0.size(); ++index)
                    {
                        ImplySolver::add_uc(O0[index],0);
                    }
                    Onp = OSequence({O0});
                }
                break;
            }    

            case(remem_Ok):
            {
                const int boundK = 3;
                // If we remember all, it will be too many.
                // NOTE: sonly implement backward now.
                assert(backwardCAR);
                
                if (backwardCAR)
                {
                    // Ub[0] = I;
                    updateU(init, nullptr);
                    pickStateLastIndex = Ub.size();

                    Onp = last_chker->Onp;
                    Onp.resize(boundK);

                    for(size_t findex = 0; findex < Onp.size(); ++findex)
                    {
                        auto &frame = Onp[findex];
                        for(size_t index = 0; index < frame.size(); ++index)
                        {
                            ImplySolver::add_uc(frame[index],findex);
                        }
                        getMainSolver(0)->add_new_frame_M(Onp[0], Onp.size() - 1);
                    }

                    if(get_rotate())
                    {
                        rotates.push_back(init->get_latches());
                    }
                    return false;
                }
                break;
            }

            case(remem_Uk):
            {
                assert(backwardCAR);
                delete init;
                auto *&last_init = last_chker->Ub[0];
                init = new State(last_init->get_inputs(), last_init->get_latches());
                clear_defer(init);

                if (immediateCheck(init, res, O0))
                    return true;

                updateU(init, nullptr);
                // backward inits.
                // Ub[0] = I;   
                for(car::State *st: last_chker->Ub){
                    if(!st)
                        break;
                    // all its followings
                    if(last_chker->whichPrior()[st] == init)
                    {
                        State* s = new State(st->get_inputs(), st->get_latches());
                        clear_defer(s);
                        updateU(s, init);
                    }
                }

                pickStateLastIndex = Ub.size();

                // Ob[0] = uc0.
                // clauses will be added by immediate_satisfible.
                // SAT_assume(init, ~p)
                // uc from init
                // ~p in ~uc
                // use uc to initialize O[0] is suitable.
                for(size_t index = 0; index< O0.size(); ++index)
                {
                    ImplySolver::add_uc(O0[index],0);
                }
                Onp = OSequence({O0});
                break;
            }

        default:
            if (immediateCheck(init, res, O0))
                return true;
            if (!backwardCAR)
            {
                // Uf[0] = ~p;

                // NOTE: we do not explicitly construct Uf[0] data strcutre. Because there may be too many states.  We actually need concrete states to get its prime. Therefore, we keep an StartSolver here, and provide a method to enumerate states in Uf[0].
                // construct an abstract state
                updateU(negp,nullptr);

                // O_I
                // NOTE: O stores the UC, while we actually use negations.
                OFrame frame;
                auto& latches = init->get_latches();
                frame.resize(latches.size());
                std::transform(latches.begin(), latches.end(), frame.begin(), [](auto lit) {
                    return std::vector<int>{-lit};
                });
                OI = {frame};
            }
            // backward inits.
            if (backwardCAR)
            {
                // Ub[0] = I;
                updateU(init, nullptr);

                // Ob[0] = uc0.
                // clauses will be added by immediate_satisfible.
                // SAT_assume(init, ~p)
                // uc from init
                // ~p in ~uc
                // use uc to initialize O[0] is suitable.
                for(const auto& uc0:O0)
                {
                    ImplySolver::add_uc(uc0,0);
                }
                Onp = OSequence({O0});
            }
            break;
        }



        if(get_rotate())
        {
            rotates.push_back(init->get_latches());
        }

        if (backwardCAR)
            getMainSolver(0)->add_new_frame_M(Onp[0], Onp.size() - 1);
        if (!backwardCAR)
            getMainSolver(0)->add_new_frame_M(OI[0], OI.size() - 1);
        if (propMode > 0 && backwardCAR)
            prop_solver->add_new_frame_P(Onp[0], Onp.size() - 1);


        return false;
    }

    //////////////helper functions/////////////////////////////////////////////

    // NOTE: if not updated, it return the same state all the time?
    State *Checker::enumerateStartStates()
    {
        if(partial)
        {
            if (start_solver->solve_with_assumption())
            {
                Assignment ass = start_solver->get_model();
                ass.resize(model_->num_inputs() + model_->num_latches());
                State *partial_res = get_partial_state(ass, nullptr);
                clear_defer(partial_res);
                return partial_res;
            }
        }
        else
        {
            if (start_solver->solve_with_assumption())
            {
                State *res = start_solver->create_new_state();
                clear_defer(res);
                return res;
            }
        }
        return NULL;
    }

    // This is used in sequence initialization
    bool Checker::immediateCheck(State *from, bool &res, OFrame &O0)
    {
        MainSolver* ms = getMainSolver(0);
        // NOTE: init may not be already set.
        vector<int> latches = from->get_latches();

        int last_max = 0;
        do
        {
            CARStats.count_main_solver_original_time_start();
            if (ms->badcheck(latches, bad_))
            {
                CARStats.count_main_solver_original_time_end(true,0);
                // if sat. already find the cex.
                State *s = ms->get_state(); // no need to shrink here
                clear_defer(s);
                whichPrior()[s] = from;
                whichCEX() = s;
                res = false;
                return true;
            }
            // NOTE: the last bit in uc is added in.
            Cube cu = ms->get_conflict_no_bad(bad_); // filter 'bad'
            CARStats.count_main_solver_original_time_end(false,cu.size());
            
            if (cu.empty())
            {
                // this means, ~p itself is bound to be UNSAT. No need to check.
                res = true;
                return true;
            }

            if(convMode==0 || convMode == 1)
            {
                // this time's first lit.
                if (abs(cu[0]) <= last_max)
                    break;
                else
                {
                    last_max = abs(cu[0]);
                    O0.push_back(cu);
                }
                auto unfresh = [&cu, &last_max](int x)
                {
                    return abs(x) > last_max;
                };
                std::stable_partition(latches.begin(), latches.end(), unfresh);
            }
            else{
                O0.push_back(cu);
                break;
            }
        } while (true);
        return false;
    }

    bool Checker::finalCheck(State *from)
    {
        bool direction = !backwardCAR;
        MainSolver *ms = getMainSolver(0);
        if (direction)
        // in forward CAR, last target state is Init. Initial state is a concrete state, therefore it is bound to be reachable (within 0 steps). Here we only need to update the cex structure.
        {
            whichCEX() = from;
            return true;
        }
        else
        // in backward car, we use uc0 instead of ~P to initialize Ob[0]. This makes it possible to return false, because uc0 is necessary but not essential.
        {
            // check whether it's in it.
            bool res = ms->badcheck(from->get_latches(), bad_);
            if (res)
            {
                // OK. counter example is found.    
                State *s = ms->get_state();
                clear_defer(s);
                whichPrior()[s] = from;
                whichCEX() = s;
                return true;
            }
        }
        return false;
    }

    /**
     * @brief
     * @pre s->get_latches() is in abs-increasing order
     *
     * @param s
     * @param frame_level
     * @param O
     * @return true
     * @return false
     */
    bool Checker::blockedIn(State *s, const int frame_level)
    {
        bool res = false;
        switch (impMethod)
        {
            case(Imp_Manual):
            {
                OFrame &frame = whichFrame(frame_level);
                for (const auto &uc : frame)
                {
                    res = s->imply(uc);
                    if (res)
                    {
                        break;
                    }
                }
                break;
            }

            case(Imp_Solver):
            {
                res = ImplySolver::is_blocked(s,frame_level);
                break;   
            }

            case(Imp_Sample):
            {
                // -1: not decided
                // 0: use solver
                // 1: manually 

                if(imply_decision != 1)
                {
                    // use solver
                    if(imply_decision == -1)
                        CARStats.count_imply_dec_begin();
                    res = ImplySolver::is_blocked(s,frame_level);
                    if(imply_decision == -1)
                        CARStats.count_imply_dec_end(1);
                }

                // manually.
                if(imply_decision != 0)
                {
                    if(imply_decision == -1)
                        CARStats.count_imply_dec_begin();
                    OFrame &frame = whichFrame(frame_level);
                    for (const auto &uc : frame)
                    {
                        res = s->imply(uc);
                        if (res)
                        {
                            break;
                        }
                    }
                    if(imply_decision == -1)
                        CARStats.count_imply_dec_end(2);
                }
                // do dicision
                if(imply_decision == -1)
                    imply_decision = CARStats.imply_dec_decide();
                
                break;
            }


            case (Imp_Exp):
            {
                // count for sz and winning rate.

                // solver
                clock_high begin = steady_clock::now();
                res = ImplySolver::is_blocked(s,frame_level);
                clock_high end = steady_clock::now();
                duration_high elapsed = end - begin;
                double time_delay_for_solver = elapsed.count(); 

                // manual
                begin = steady_clock::now();
                
                OFrame &frame = whichFrame(frame_level);
                for (const auto &uc : frame)
                {
                    res = s->imply(uc);
                    if (res)
                    {
                        break;
                    }
                }
                end = steady_clock::now();
                elapsed = end - begin;
                bool solver_win = elapsed.count()> time_delay_for_solver ? true : false;
                CARStats.record_winner(solver_win,frame.size());

                break;
            }

            case (Imp_Thresh):
            {
                OFrame &frame = whichFrame(frame_level);
                if(frame.size() > 10000)
                {
                    res = ImplySolver::is_blocked(s,frame_level);
                }
                else{
                    for(int i = frame.size()-1; i>=0; --i)
                    {
                        const auto &uc = frame[i];
                        res = s->imply(uc);
                        if (res)
                        {
                            break;
                        }
                    }
                }
                break;
            }

            case (Imp_Fresh):
            {
                /**
                 * @brief Record what has been tested as to each state
                 * 
                 */
                // get its fresh index!
                if(impFreshRecord.find(frame_level) == impFreshRecord.end())
                {
                    // this level has not been tested before.
                    impFreshRecord[frame_level] = std::unordered_map<int,int>({});
                }
                auto &this_level_map = impFreshRecord[frame_level];

                // fresh index: the index we should start testing
                int freshIndex = 0;
                if(this_level_map.find(s->id) != this_level_map.end())
                {
                    // found the fresh index! We could start from the last time
                    freshIndex = this_level_map[s->id];
                }

                auto &frame = whichFrame(frame_level);
                int i = freshIndex;
                for (; i<frame.size(); ++i)
                {
                    auto& uc = frame[i];
                    res = s->imply(uc);
                    if (res)
                    {
                        break;
                    }
                }
                if(res)
                    this_level_map[s->id] = i;          // the one that blocks it
                else
                    this_level_map[s->id] = frame.size();  // end of the frame.
                // cerr<<s->id<<"@"<<frame_level<<":"<<this_level_map[s->id]<<" ";

                break;
            }

            default:
                break;
        }
        return res;
    }

    /**
     * @brief Find minimal level among [min,max] that s is not blocked in. 
     * 
     * @param s 
     * @param min 
     * @param max 
     * @param O 
     * @return int 
     */
    int Checker::minNOTBlocked(State *s, const int min, const int max)
    {
        CARStats.count_imply_begin();
        int start = min;
        while (start <= max)
        {
            if (!blockedIn(s, start))
            {
                break;
            }
            start++;
        }
        CARStats.count_imply_end();
        return start;
    }
}