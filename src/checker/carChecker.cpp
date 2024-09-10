#include "carChecker.h"
#include "implysolver.h"
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
    out(out), model_(model), rotate_enabled(opt.enable_rotate), inter_cnt(opt.inter_cnt), inv_incomplete(opt.inv_incomplete), uc_no_sort(opt.raw_uc), impMethod(opt.impMethod), time_limit_to_restart(opt.time_limit_to_restart), rememOption(opt.rememOption), LOStrategy(opt.LOStrategy), convAmount(opt.convAmount), convParam(opt.convParam), convMode(opt.convMode), subStat(opt.subStat), partial(opt.partial), last_chker(last_chker), fresh_levels(0), backwardCAR(!opt.forward), bad_(model->output(0)), inv_solver(nullptr), multi_solver(opt.multi_solver),propMode(opt.propMode), propParam(opt.propParam), simp(opt.simplifyCNF), container_option(opt.container_option)
    {
        if(!multi_solver)
        {
            main_solver = new ModelSolver(model, simp);
            prop_solver = new ModelSolver(model, simp);
        }    

        rotates.clear();
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
        if (prop_solver)
        {
            delete prop_solver;
            prop_solver = nullptr;
        }
    }

    RESEnum Checker::check()
    {
        RESEnum res = trivialCheck();
        if (res == RES_UNKNOWN)
        {
            res = car();
        }

        switch (res)
        {
        case RES_SAFE:
            out << "0" << endl;
            out << "b0" << endl;
            out << "." << endl;
            break;

        case RES_UNSAFE:
            print_evidence();

        case RES_UNKNOWN:
            break;

        case RES_RESTART:
            break;

        default:
            assert(false);
        }

        return res;
    }

    ModelSolver *Checker::getMainSolver(int level)
    {
        if(!multi_solver)
            return main_solver;
        else
        {
            if(level >= main_solvers.size())
            {
                int needs = level - main_solvers.size() + 1;
                while(needs > 0)
                {
                    auto *new_slv = new ModelSolver(model_,simp);
                    main_solvers.push_back(new_slv);                        
                    --needs;
                }
            }
            if(level <=0)
                return main_solvers.at(0);
            else
                return main_solvers.at(level);
        }
    }

    RESEnum Checker::car()
    {
        RESEnum res = initSequence();
        if (RES_UNKNOWN != res)
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
            res = trySAT();
            if(RES_UNKNOWN != res)
                return res;

            if (!inv_incomplete)
            {
                if (InvFound())
                    return RES_SAFE;
            }
        }

        // dead code. Should not reach
        return RES_FAIL;
    }

    RESEnum Checker::trivialCheck()
    {
        const Problem *model = model_;
        // FIXME: fix for multiple properties.
        if (bad_ == model->true_id())
        {

            vector<int> all0(model_->num_inputs(),0);
            State *init = new State(model_->init());
            State *neg = new State(all0, model_->init());
            clear_defer(init);
            clear_defer(neg);
            whichPrior()[neg] = init;
            whichCEX() = neg;
            return RES_UNSAFE;
        }
        else if (bad_ == model->false_id())
        {
            out << "0" << endl;
            out << "b" << endl;
            out << "." << endl;
            return RES_SAFE;
        }
        return RES_UNKNOWN;
    }

    RESEnum Checker::trySAT()
    {
        // NOTE: can eliminate initialization.
        auto& O = whichO();
        bool safe_reported = false;

        // add a new frame to O, as Otmp
        O.push_back({});

        // if rotate enabled, push back a new vector to record.
        if (rotate_enabled)
            rotates.push_back({});

        CARStats.count_enter_new_ronud();
        /**
         * this procedure is like the old car procedure, but the OSequence is not bound to be OI or Onegp.
         * @param missionary the state to be checked
         */
        while (State *missionary = pickState())
        {
            /**
             * build a container <state, depth, target_level>
             */
            containerClear();
            containerPush(Obligation(missionary, 0, OSize() - 1));
            while (!containerEmpty())
            {
                CARStats.count_enter_new_try_by();
                State *s;
                int dst, depth;
                std::tie(s, depth, dst) = containerTopOrFront();
                if (blockedIn(s, dst + 1))
                {
                    containerPop();
                    CARStats.count_tried_before();

                    int new_level = minNOTBlocked(s, dst + 2, OSize() - 1);
                    if (new_level <= OSize())
                    {
                        containerPush(Obligation(s, depth, new_level - 1));
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
                        return RES_RESTART;
                    }
                }

                if (satAssume(s, dst, safe_reported))
                {
                    if (dst == -1)
                    {
                        return RES_UNSAFE;
                    }
                    State *tprime = getModel(s, dst);

                    // NOTE: why here calculate minNOTBLOCKED, rather than next time when pop?
                    int new_level = minNOTBlocked(tprime, 0, dst - 1);               
                    if (new_level <= dst) // if even not one step further, should not try it
                    {
                        containerPush(Obligation(tprime, depth + 1, new_level - 1));
                    }
                }
                else
                {
                    containerPop();
                    if (safe_reported)
                        return RES_FAIL;

                    int new_level = minNOTBlocked(s, dst + 2, OSize() - 1);
                    if (new_level <= OSize())
                    {
                        containerPush(Obligation(s, depth, new_level - 1));
                    }
                }
            }
        }

        return RES_UNKNOWN;
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
        return s;
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
                    helper.push(to_print->getInputsStr());
                else
                    helper.push(to_print->getLatchesStr());
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
                out << to_print->getLatchesStr() << endl;
            while (to_print)
            {
                State *next = prior_in_trail_f.find(to_print) != prior_in_trail_f.end() ? prior_in_trail_f.at(to_print) : nullptr;
                out << to_print->getInputsStr() << endl;
                to_print = next;
            }
        }
        out << "." << endl;
    }


    Cube Checker::reorderAssumptions(State*s, int level, const vector<Cube>& inter, const Cube &rcube, const Cube &rest)
    {
        switch (LOStrategy)
        {
            case LO_Rand: // random
            {
                std::random_device rd;
                std::mt19937 gen(rd());
                vector<int> shuff_helper = s->getLatches();
                for(int i = shuff_helper.size() - 1; i > 0; --i) {
                    // starting from 1, because 0 is the flag.
                    std::uniform_int_distribution<int> dis(0, i);
                    int j = dis(gen);
                    std::swap(shuff_helper[i], shuff_helper[j]);
                }
                return shuff_helper;
            }
            
            case LO_Classic_Pos:
            {
                if(level == -1)
                    return s->getLatches();
                Cube pref;
                for(const auto& cu : inter)
                    pref.insert(pref.end(), cu.begin(), cu.end());
                pref.insert(pref.end(), rcube.begin(), rcube.end());
                pref.insert(pref.end(), rest.begin(), rest.end());
                vector<int> sl = s->getLatches();
                pref.insert(pref.end(),sl.begin(),sl.end());
                return pref;
            }
            
            case LO_Classic: // traditional one, no need to shuffel
            default:
            {                
                Cube pref;
                for(const auto& cu : inter)
                    pref.insert(pref.end(), cu.begin(), cu.end());
                pref.insert(pref.end(), rcube.begin(), rcube.end());
                pref.insert(pref.end(), rest.begin(), rest.end());
                vector<int> sl = s->getLatches();
                pref.insert(pref.end(),sl.begin(),sl.end());
                return pref;
            }
        }
    }

    bool Checker::satAssume(State *s, int level, bool& safe_reported)
    {
        ModelSolver* ms = getMainSolver(level);
        bool res = false;

        vector<Cube> inter;
        get_inter_res(s,level, inter);
        Cube rcube, rest;
        get_rotate_res(s,level,rcube,rest);
        auto asm_vec = reorderAssumptions(s, level, inter, rcube, rest);

        CARStats.count_main_solver_original_time_start();
        if(level == -1)
            res = ms->zeroStepReachable(asm_vec, bad_);
        else
            res = ms->oneStepReachable(asm_vec, level, backwardCAR);

        if(res)
        {
            CARStats.count_main_solver_original_time_end(res,0); //uc size == 0
            if(level == -1)
            {
                // TODO: consider to use immediateCheck() here, to get more than 1 uc when level == -1.
                assert(backwardCAR);
                State *cex = getModel(s, level);
                whichCEX() = cex;
                return true;
            }
        }
        else
        {
            Cube uc = ms->getUCofLatch(!backwardCAR);
            CARStats.count_main_solver_original_time_end(res,uc.size());
            
            if (uc.empty())
            {
                safe_reported = true;
                cerr<<"uc is empty@1!"<<endl;
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

            // add uc to solver.
            // note: should do this before propagation.
            addUCtoSolver(uc, level + 1);

            if (propMode > 0 && level > 0 && level + 1 <= OSize())
                propagate(uc, level);
            
            if (uc.empty())
            {
                // watch dog for propagation.
                cerr<<"uc is empty@2!"<<endl;
                safe_reported = true;
            }

            if(rotate_enabled)
            {
                // update rotate
                Cube &rcu = rotates[level + 1];
                rcu = rcube;
                rcu.insert(rcu.end(), rest.begin(), rest.end());
            }

            if(convMode >= 0)
            {
                mUC(uc, level);
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
     * @todo record whether the UC to propagate is already within lower frame.
     * @pre  should be right after an UNSAT call.
     */
    void Checker::propagate(Cube& uc, int level)
    {
        if(level <= 0 ||  level + 1 > OSize())
            return;

        static std::map<int,int> prop_tick;
        ++prop_tick[level];

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
            case PropFresh:
            {
                if(level < OSize()-1)
                    return;
                break;
            }
            case PropShortFresh:
            {
                if(uc.size() > propParam)
                    return;
                if(level < OSize()-1)
                    return;
                break;
            }
            case PropEarlyStop:
            {
                if(prop_tick[level] > 2048) // a magic number.
                    return;
                break;
            }
            case PropEarlyShortFresh:
            {
                if(uc.size() > propParam)
                    return;
                if(level < OSize()-1)
                    return;
                if(prop_tick[level] > 2048)
                    return;
                break;
            }
            
            default: // fall through.
            case PropNone:
            {
                return;
            }
        }

        CARStats.count_prop_begin();
        prop_solver->SATslv->clear_assumption();
        
        vector<int> assumption; 
        // assumption ::= <Oflag of l+1>, <deactivation flags>, <activation flag>,<uc'>
        assumption.reserve(1 + FlagManager::PTFlags.size() + 1 + uc.size());

        int Oflag = FlagManager::FlagOf(level+1, SolverFlag_Main);
        assumption.push_back(Oflag);

        // deactivate all the previous ones.
        for(int flg: FlagManager::PTFlags)
        {
            assumption.push_back(flg);
        }
        int new_flag = FlagManager::getNewPTFlag();
        // activate the current one:
        assumption.push_back(-new_flag);

        // the UC(will be primed later in `set_assumption_primed`):
        assumption.insert(assumption.end(), uc.begin(), uc.end());

        // get its primed version.
        for(auto& lit: assumption)
        {
            if(model_ -> latch_var(lit))
                lit = model_->prime(lit);
        }
        
        prop_solver->SATslv->push_assumption(assumption); 

        // FIXME: CHECK ME.
        // clause ::= !uc , current flag.
        prop_solver->addNotCubeToLevel(uc,level,SolverFlag_PropTemp, false); // !uc
        
        // check if !uc is inductive
        bool prop_res = prop_solver->SATslv->solve_assumption();

        // do minimization of UC:
        bool strongInductive = false;

        int originalUCsize = uc.size();

        if(!prop_res)
        {
            auto upcoming_uc = prop_solver->SATslv->get_uc();
            /**
             * The upcoming UC ::= <subset of uc> [<activation flag>] [<Oflag>]
             * Reason:
             * 1) O flag is the first in the assumption vector;
             *  either it appears ==> should be in the back
             *  or it does not appear
             * 2) The deactivation literals are pure literals, should not appear.
             * 3) The activation flag is likely to appear -- Otherwise, T /\ uc' is tautology.
            */

            strongInductive = false;// regardless of the frame.

            // again, remember:
            // assumption ::= <Oflag of l+1>, <deactivation flags>, <activation flag>,<uc'>
            // assert: flag should only appear in the back, or not.
            if(FlagManager::LevelOf(upcoming_uc.back(),SolverFlag_Main) == NOT_M_FLAG)
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
                // we could check whether it could be propagate further.
                if(propMode == PropContinue || propMode == PropShortCont)
                    propagate(uc,level - 1);
            }
        }

        CARStats.count_prop_end(!prop_res,strongInductive, originalUCsize, prop_tick[level]);
    }

    bool Checker::convTriggered(int level)
    {
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
        return trigger;
    }

    Cube Checker::getNextAsm()
    {
        // FIXME: fill in here
        return Cube({});
    }

    void Checker::mUC(const Cube &uc, int level)
    {
        if(convMode < 0)
            return;

        // whether to calculate another UC.
        bool trigger = convTriggered(level);
        if (trigger)
        {
            vector<Cube> ucs; // this is used for subsumption test, will cause low-efficiency, do not turn it on in practical scenarios.
            if(subStat)
            {
                ucs.push_back(uc);
                CARStats.recordUC(false);
            }
            int amount = 0;
            while(++amount <= convAmount)
            {
                // retrieve another bit.
                conv_record[level + 1] <<= 1;

                // mark as inserted
                conv_record[level + 1] += 1;

                CARStats.count_main_solver_convergence_time_start();
                // get another conflict!

                ModelSolver* ms = getMainSolver(level);
                Cube next_asm = getNextAsm();
                if(level == -1)
                    ms->zeroStepReachable(next_asm, bad_);
                else
                    ms->oneStepReachable(next_asm,level, backwardCAR);
                auto nextuc = ms->getUCofLatch(backwardCAR);
                
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
        }
    }

    void Checker::get_inter_res(State *s, int level, vector<Cube>&inter)
    {        
        inter.clear();
        if (!get_inter_cnt()) return;

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
        return;
    }

    void Checker::get_rotate_res(State *s, int level, Cube &rcube, Cube &rest)
    {
        rcube.clear();
        rest.clear();
        if(!rotate_enabled) return;

        // rotates[i]
        Cube &rcu = rotates[level + 1];
        if (rcu.empty())
        {
            rcube = {};
            rest = s->getLatches();
            // TODO: try this to be inter?
            return;
        }

        // a full state
        if (!s->isPartial())
        {
            // calculate intersection and put the others behind.
            for(size_t i = 0; i < rcu.size(); ++i)
            {
                if (s->getLatchAt(abs(rcu[i]) - model_->num_inputs() - 1) == rcu[i])
                    rcube.push_back(rcu[i]);
                else
                    rest.push_back(-rcu[i]);
            }
        }
        else
        {
            // a partial state
            std::unordered_set<int> helper;
            for(auto lit: s->getLatches())
            {
                helper.insert(lit);
            }
            for(auto l: rcu)
            {
                if(helper.find(l) != helper.end())
                    rcube.push_back(l);
                else
                    rest.push_back(l);
            }                 
        }
    }

    State *Checker::getModel(State *s, int level)
    {
        State *t = getMainSolver(level)->getState(backwardCAR);
        clear_defer(t);
        whichU().push_back(t);
        if(s->isPartial())
        {
            // if predecessor is partial, we would like to initialize it.
            // why not s->set_latches() ? because s may be reused, but newInit will not be used in the search(but only in CEX printing).
            State* newInit = getMainSolver(level)->getState(false);
            clear_defer(newInit);
            whichPrior()[t] = newInit;
        }
        else{
            whichPrior()[t] = s;
        }
        // NOTE: if it is the last state, it will be set last-inputs later.
        return t;
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

        if (dst_level_plus_one <= OSize())
        {
            getMainSolver(dst_level_plus_one)->addNotCubeToLevel(uc,dst_level_plus_one,SolverFlag_Main, backwardCAR);
            if(propMode > 0)
            {
                // also update the propagation solver
                prop_solver->addNotCubeToLevel(uc,dst_level_plus_one,SolverFlag_Prop,!backwardCAR);
            }
        }
        
        if(impMethod == Imp_Bit || impMethod == Imp_BitFresh)
        {
            while(ucs_masks.size() <= dst_level_plus_one + 1)
            {
                ucs_masks.push_back({});
            }
            auto& masks = ucs_masks[dst_level_plus_one];
            const int nlatches = model_->num_latches();
            const int startPos = model_->num_inputs() + 1;
            masks.push_back(UCMask(uc,nlatches,startPos));
        }
    }

    // TODO: rewrite this part.
    RESEnum Checker::initSequence()
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
                RESEnum res = immediateCheck(init, O0);
                if (RES_UNKNOWN != res)
                    return res;
                // backward inits.
                assert(backwardCAR);
                {
                    // Ub[0] = I;
                    whichPrior()[init] = nullptr;
                    whichU().push_back(init);
                    pickStateLastIndex = Ub.size();

                    // Ob[0] = uc0.
                    // clauses will be added by immediate_satisfible.
                    // SAT_assume(init, ~p)
                    // uc from init
                    // ~p in ~uc
                    // use uc to initialize O[0] is suitable.
                    for(size_t index = 0; index< O0.size(); ++index)
                    {
                        if(impMethod == Imp_Bit || impMethod == Imp_BitFresh)
                        {
                            ucs_masks.push_back({});   
                            auto& masks = ucs_masks[0];
                            const int nlatches = model_->num_latches();
                            const int startPos = model_->num_inputs() + 1;
                            masks.push_back(UCMask(O0[index],nlatches,startPos));
                        }
                    }
                    Onp = OSequence({O0});
                }
                break;
            }

            case(remem_O0):
            {
                // If we remember all, it will be too many.
                // NOTE: only implement backward now.
                assert(backwardCAR);
                O0 = last_chker->Onp[0];
                
                if (backwardCAR)
                {
                    // Ub[0] = I;
                    whichPrior()[init] = nullptr;
                    whichU().push_back(init);
                    pickStateLastIndex = Ub.size();

                    // Ob[0] = uc0.
                    // clauses will be added by immediate_satisfible.
                    // SAT_assume(init, ~p)
                    // uc from init
                    // ~p in ~uc
                    // use uc to initialize O[0] is suitable.

                    Onp = OSequence({O0});
                }
                break;
            }

        default:
            assert(false && "should have missed something in remember\n");
        }

        if(get_rotate())
        {
            rotates.push_back(init->getLatches());
        }

        if (backwardCAR)
        {
            auto ms = getMainSolver(0);
            for(const auto& cu: Onp[0])
            {
                ms->addNotCubeToLevel(cu,Onp.size() - 1,SolverFlag_Main,backwardCAR);
            }
            if (propMode > 0 && backwardCAR)
            for(const auto& cu: Onp[0])
            {
                prop_solver->addNotCubeToLevel(cu,Onp.size() - 1,SolverFlag_Prop,backwardCAR);
            }
        }


        return RES_UNKNOWN;
    }

    //////////////helper functions/////////////////////////////////////////////


    // This is used in sequence initialization
    RESEnum Checker::immediateCheck(State *from, OFrame &O0)
    {
        ModelSolver* ms = getMainSolver(0);
        // NOTE: init may not be already set.
        vector<int> latches = from->getLatches();

        int last_max = 0;
        do
        {
            CARStats.count_main_solver_original_time_start();
            if (ms->zeroStepReachable(latches, bad_))
            {
                CARStats.count_main_solver_original_time_end(true,0);
                // if sat. already found the cex.
                State *cex = getModel(from, 0);
                whichCEX() = cex;
                return RES_UNSAFE;
            }
            // NOTE: the last bit in uc is added in.
            Cube cu = ms->getUCofLatch(false); // filter 'bad'
            CARStats.count_main_solver_original_time_end(false,cu.size());
            
            if (cu.empty())
            {
                // this means, ~p itself is bound to be UNSAT. No need to check. But it should not happen in most cases.
                return RES_SAFE;
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
        return RES_UNKNOWN;
    }

    /**
     * @brief
     * @pre s->getLatches() is in abs-increasing order
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

            case (Imp_Bit):
            {
                UCMask u_s(s->getLatches());
                // get frame_level.
                if(frame_level >= ucs_masks.size())
                {
                    res = false;
                    break;
                }
                const auto& masks_this_level = ucs_masks.at(frame_level);
                for(const auto& m: masks_this_level)
                    if(u_s.imply(m))
                    {
                        res = true;
                        break;
                    }
                break;
            }

            case (Imp_BitFresh):
            {
                /**
                 * @brief Record what has been tested as to each state
                 * 
                 */

                // no uc yet. Return False.
                if(frame_level >= ucs_masks.size())
                {
                    res = false;
                    break;
                }

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

                
                UCMask u_s(s->getLatches());  // remember this also?

                const auto& masks_this_level = ucs_masks.at(frame_level);
                size_t framesz = whichFrame(frame_level).size();

                int i = freshIndex;
                for (; i<framesz; ++i)
                {
                    auto& uc = masks_this_level[i];
                    res = u_s.imply(uc);
                    if (res)
                    {
                        break;
                    }
                }
                if(res)
                    this_level_map[s->id] = i;          // the one that blocks it
                else
                    this_level_map[s->id] = framesz;  // end of the frame.

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