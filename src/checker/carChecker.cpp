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
    model_(model), out(out), rotate_enabled(opt.enable_rotate), inter_cnt(opt.inter_cnt), inv_incomplete(opt.inv_incomplete), uc_no_sort(opt.raw_uc), impMethod(opt.impMethod), time_limit_to_restart(opt.time_limit_to_restart), rememOption(opt.rememOption), LOStrategy(opt.LOStrategy), convAmount(opt.convAmount), convParam(opt.convParam), convMode(opt.convMode), subStat(opt.subStat), partial(opt.partial), last_chker(last_chker), fresh_levels(0), backward_first(!opt.forward), bad_(model->output(0)), inv_solver(nullptr)
    {
        main_solver = new MainSolver(model, opt.forward, get_rotate(), uc_no_sort);
        if (!backward_first)
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
                for (int i = 0; i < model->num_latches(); ++i)
                    out << "0";
                out << endl;
                // print an arbitary input vector
                for (int j = 0; j < model->num_inputs(); j++)
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
                    State *tprime = getModel();
                    
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
        main_solver->add_new_frame(Otmp, OSize() - 1);
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
        for (int i = 0; i < OSize(); ++i)
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
     * @brief print the evidence. reuse backward_first, which exactly reveals present searching direciton.
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
        auto& O = whichO();
        bool forward = !backward_first;
        vector<Cube> inter;
        Cube rres, rtmp;
        vector<Cube> ucs; // this is used for subsumption test, will cause low-efficiency, do not turn it on in practical scenarios.

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

                    while (index <= get_inter_cnt() && frame.size() >= index)
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
                Cube &rcu = level + 1 < rotates.size() ? rotates[level + 1] : rotate;
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
                    for (int i = 0; i < rcu.size(); ++i)
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
                    while (i < s->latch_size() && j < rcu.size())
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
                for (int i = shuff_helper.size() - 1; i > 0; --i) {
                    // starting from 1, because 0 is the flag.
                    std::uniform_int_distribution<int> dis(0, i);
                    int j = dis(gen);
                    std::swap(shuff_helper[i], shuff_helper[j]);
                }
                main_solver->set_assumption(s, level,  {shuff_helper});
                break;
            }

            case 3: // expectation test
            {
                // FIXME: implement this part to be a real decision procedure.
                std::vector<int> exp = s->get_latches(); // tmp only.
                main_solver->set_expectation(exp);
                main_solver->set_assumption(s, level,  pref);
                break;
            }

            case 2:
            case 0: // traditional one, no need to shuffel
            default:
            {
                main_solver->set_assumption(s, level,  pref);
                // do nothing
                break;
            }
            }
            
            CARStats.count_main_solver_original_time_start();
            res = main_solver->solve_assumption();
            if(res)
            {
                // if sat, no uc, we just update it here.
                CARStats.count_main_solver_original_time_end(res,0);
            }
        }
        if(!res)
        {
            // update the UC.
            Cube uc = main_solver->get_conflict();

            if(subStat)
            {
                ucs.push_back(uc);
                CARStats.recordUC(false);
            }

            if (uc.empty())
            {
                cerr<<"uc is empty!"<<endl;
                safe_reported = true;
            }
            CARStats.count_main_solver_original_time_end(res,uc.size());

            addUCtoSolver(uc, level + 1);
        }

        if (get_rotate() && !res)
        {
            // update rotate
            Cube &rcu = level + 1 < rotates.size() ? rotates[level + 1] : rotate;
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
                auto nextuc = main_solver->get_conflict_another(LOStrategy, amount);
                
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

    State *Checker::getModel()
    {
        State *s = main_solver->get_state();
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
        OSequence &O = whichO();
        if (dst_level_plus_one < fresh_levels)
            fresh_levels = dst_level_plus_one;

        OFrame &frame = whichFrame(dst_level_plus_one);

        frame.push_back(uc);

        ImplySolver::add_uc(uc,dst_level_plus_one);
        

        if (dst_level_plus_one < OSize())
            main_solver->add_clause_from_cube(uc, dst_level_plus_one);
        else if (dst_level_plus_one == OSize())
        {
            if (!backward_first)
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
                if (!backward_first)
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
                if (backward_first)
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
                    for(int index = 0; index< O0.size(); ++index)
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
                assert(backward_first);
                O0 = last_chker->Onp[0];
                
                if (backward_first)
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

                    for(int index = 0; index< O0.size(); ++index)
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
                assert(backward_first);
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

                if (backward_first)
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

                    for(int index = 0; index< O0.size(); ++index)
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
                assert(backward_first);
                
                if (backward_first)
                {
                    // Ub[0] = I;
                    updateU(init, nullptr);
                    pickStateLastIndex = Ub.size();

                    Onp = last_chker->Onp;
                    Onp.resize(boundK);

                    for(int findex = 0; findex < Onp.size(); ++findex)
                    {
                        auto &frame = Onp[findex];
                        for(int index = 0; index < frame.size(); ++index)
                        {
                            ImplySolver::add_uc(frame[index],findex);
                        }
                        main_solver->add_new_frame(Onp[0], Onp.size() - 1);
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
                assert(backward_first);
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
                for(int index = 0; index< O0.size(); ++index)
                {
                    ImplySolver::add_uc(O0[index],0);
                }
                Onp = OSequence({O0});
                break;
            }

        default:
            if (immediateCheck(init, res, O0))
                return true;
            if (!backward_first)
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
            if (backward_first)
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

        if (backward_first)
            main_solver->add_new_frame(Onp[0], Onp.size() - 1);
        if (!backward_first)
            main_solver->add_new_frame(OI[0], OI.size() - 1);



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
        auto &solver = main_solver;
        // NOTE: init may not be already set.
        vector<int> latches = from->get_latches();

        int last_max = 0;
        do
        {
            CARStats.count_main_solver_original_time_start();
            if (solver->badcheck(latches, bad_))
            {
                CARStats.count_main_solver_original_time_end(true,0);
                // if sat. already find the cex.
                State *s = solver->get_state(); // no need to shrink here
                clear_defer(s);
                whichPrior()[s] = from;
                whichCEX() = s;
                res = false;
                return true;
            }
            // NOTE: the last bit in uc is added in.
            Cube cu = solver->get_conflict_no_bad(bad_); // filter 'bad'
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
        bool direction = !backward_first;
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
            bool res = main_solver->badcheck(from->get_latches(), bad_);
            if (res)
            {
                // OK. counter example is found.    
                State *s = main_solver->get_state();
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
        OSequence &O = whichO();
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
                    for (int i = frame.size()-1; i>=0; --i)
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