/**
 * @file statistics.h
 * @author yibodong (prodongf@gmail.com)
 * @brief The statistics we would like to collect
 * @version 0.1.0
 * @date 2024-07-25
 * 
 * 
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <chrono>
using namespace std::chrono;
using std::__cxx11::to_string;
using clock_high = steady_clock::time_point;
using duration_high = std::chrono::duration<double,std::milli>;


namespace car {

/**
 * @brief The statistic collector
 * 
 */
class Statistics 
{
    public:
        Statistics () {}
        ~Statistics () {}
        // main solver's information.

        /// how many are original ones
        int num_main_solver_original_calls_total = 0;
        /// how many succeeded
        int num_main_solver_original_calls_success = 0;
        /// how many failed
        int num_main_solver_original_calls_failed = 0;
        /// time spent
        double time_main_solver_original_calls_ = 0.0;
        /// uc's length in total.
        int sum_main_solver_original_uc_length_ = 0;
        clock_high main_solver_original_begin_, main_solver_original_end_;
        /// uc's distribution
        int num_uc_0_3 = 0;
        int num_uc_3_5 = 0;
        int num_uc_5_10 = 0;
        int num_uc_10_50 = 0;
        int num_uc_50_100 = 0;
        int num_uc_100_inf = 0;
        

        bool is_counting_main_original;

        inline void count_main_solver_original_time_start ()
        {
            main_solver_original_begin_ = steady_clock::now();
            is_counting_main_original = true;
        }

        int last_original_uc_sz = 0;

        inline void count_main_solver_original_time_end (bool res, int uc_length)
        {
            main_solver_original_end_ = steady_clock::now();
            duration_high elapsed = main_solver_original_end_ - main_solver_original_begin_;
            double time_delay = elapsed.count();

	        time_main_solver_original_calls_ += time_delay;
            is_counting_main_original = false;
	        num_main_solver_original_calls_total ++;
	        if(res)
            {
                ++num_main_solver_original_calls_success;
            }
            else
            {
                ++num_main_solver_original_calls_failed;
                sum_main_solver_original_uc_length_ += uc_length;
                last_original_uc_sz = uc_length;
                if(uc_length < 3)
                    ++num_uc_0_3;
                else if(uc_length < 5)
                    ++num_uc_3_5;
                else if(uc_length < 10)
                    ++num_uc_5_10;
                else if(uc_length < 50)
                    ++num_uc_10_50;
                else if(uc_length < 100)
                    ++num_uc_50_100;
                else
                    ++num_uc_100_inf;                
            }
        }

        /// convergence information:
        int num_main_solver_convergnece_calls_ = 0;
        int num_main_solver_convergence_shorter = 0;
        double time_main_solver_convergence_calls_ = 0.0;
        clock_high main_solver_convergence_calls_begin_, main_solver_convergence_calls_end_;
        int sum_main_solver_convergence_uc_length_ = 0;
        bool is_counting_convergence;

        inline void count_main_solver_convergence_time_start()
        {
            main_solver_convergence_calls_begin_ = steady_clock::now();
            is_counting_convergence = true;
        }
        inline void count_main_solver_convergence_time_end(int uc_length)
        {
            main_solver_convergence_calls_end_ = steady_clock::now();
            duration_high elapsed = main_solver_convergence_calls_end_ - main_solver_convergence_calls_begin_;
            double time_delay = elapsed.count();
            time_main_solver_convergence_calls_ += time_delay;
            is_counting_convergence = false;
            num_main_solver_convergnece_calls_ ++;
            
            sum_main_solver_convergence_uc_length_ += uc_length;
            if(uc_length > 0 && uc_length < last_original_uc_sz)
                ++num_main_solver_convergence_shorter;
        }

        /// search progress information:
        int num_rounds = 0;
        int num_try_by = 0;
        clock_high time_enter_new_try_by;
        inline void count_enter_new_ronud(){
            ++num_rounds;
        }
        inline void count_enter_new_try_by(){
            ++num_try_by;
            time_enter_new_try_by = steady_clock::now();
        }

        /// global time
        clock_high global_begin_, global_end_;
        double time_global = 0.0;
        inline void count_whole_begin() { 
            global_begin_ = steady_clock::now(); 
        }
        inline void count_whole_end()
        {
            global_end_ = steady_clock ::now();
            duration_high elapsed = global_end_ - global_begin_;
            time_global = elapsed.count();
        }

        /// tried before
        int num_tried_before = 0;
        inline void count_tried_before() { 
            clock_high tried_before_end = steady_clock::now();
            duration_high elapsed = tried_before_end - time_enter_new_try_by;
            time_imply += elapsed.count();
            ++num_tried_before; 
        }

        /// time for implication calculation
        clock_high imply_begin_, imply_end_;
        int count_imply = 0;
        double time_imply = 0.0;
        bool is_counting_imply;
        inline void count_imply_begin()
        {
            imply_begin_ = steady_clock::now();
            is_counting_imply = true;
        }
        inline void count_imply_end()
        {
            imply_end_ = steady_clock::now();
            duration_high elapsed = imply_end_ - imply_begin_;
            double time_delay = elapsed.count();
            time_imply += time_delay;
            count_imply++;
            is_counting_imply = false;
        }

        clock_high imply_dec_begin, imply_dec_end;
        double time_imply_dec_sol, time_imply_dec_man;
        int imply_decision = -1;
        inline void count_imply_dec_begin()
        {
            imply_dec_begin = steady_clock::now();
        }
        inline void count_imply_dec_end(int num)
        {
            imply_dec_end = steady_clock::now();
            duration_high elapsed = imply_dec_end - imply_dec_begin;
            double time_delay = elapsed.count();
            switch (num)
            {
            case 1:
                time_imply_dec_sol+= time_delay;
                break;
            case 2:
                time_imply_dec_man+= time_delay;
                break;     
            default:
                break;
            }
        }

        int imply_dec_cnter = 1000;
        void reset_imply_cnter()
        {
            imply_dec_cnter = 1000;
        }

        inline int imply_dec_decide(){
            imply_dec_cnter--;
            /// not the decision time
            if(imply_dec_cnter!=0)
                return -1;
            /// sample too few            
            if(std::max(time_imply_dec_sol,time_imply_dec_man*10) < 100)
            {
                imply_dec_cnter = 100;
                return -1;
            }/// decide.
            
            imply_decision = time_imply_dec_sol > (time_imply_dec_man*10) ? 1 : 0;
            return imply_decision;
        }

        std::map<int,int> solverWin, manualWin;
        
        inline void record_winner(bool isSolver, int framesz)
        {
            auto& ref = isSolver ?solverWin : manualWin;
            ref[framesz]++;
        }

        inline void serialize_winner(){
            std::cout<<"{"<<std::endl;
                std::cout<<"\"Solver\": {"<<std::endl;
                    for(auto& pr: solverWin){
                        std::cout<<"\""<<pr.first<<"\":"<<pr.second<<","<<std::endl;
                    }
                std::cout<<"},"<<std::endl;

                std::cout<<"\"Manual\": {"<<std::endl;
                    for(auto& pr: manualWin){
                        std::cout<<"\""<<pr.first<<"\":"<<pr.second<<","<<std::endl;
                    }
                std::cout<<"}"<<std::endl;
            std::cout<<"},"<<std::endl;
        }

        /// statistics about subsumption.
        int nUC = 0, nSubUC = 0;
        inline void recordUC(bool sub)
        {
            nUC++;
            nSubUC += sub;
        }

        /// propagation time
        clock_high prop_begin_, prop_end_;
        double time_prop = 0.0;
        int ntime_prop = 0;
        int ntime_prop_succ = 0;
        int ntime_prop_pure = 0;
        int num_prop_uc_0_3 = 0;
        int num_prop_uc_3_5 = 0;
        int num_prop_uc_5_10 = 0;
        int num_prop_uc_10_50 = 0;
        int num_prop_uc_50_100 = 0;
        int num_prop_uc_100_inf = 0;
        inline void count_prop_begin() { 
            prop_begin_ = steady_clock::now(); 
        }
        inline void count_prop_end(bool succeed, bool pureInductive, int uc_length)
        {
            ++ntime_prop;
            if(succeed)
            {
                ++ntime_prop_succ;
                if(pureInductive)
                    ++ntime_prop_pure;
                if(uc_length < 3)
                    ++num_prop_uc_0_3;
                else if(uc_length < 5)
                    ++num_prop_uc_3_5;
                else if(uc_length < 10)
                    ++num_prop_uc_5_10;
                else if(uc_length < 50)
                    ++num_prop_uc_10_50;
                else if(uc_length < 100)
                    ++num_prop_uc_50_100;
                else
                    ++num_prop_uc_100_inf;        
            }
            
            prop_end_ = steady_clock ::now();
            duration_high elapsed = prop_end_ - prop_begin_;
            time_prop += elapsed.count();
        }

        /// helper for debug.
        clock_high begin_1, end_1;
        double time_1_1 = 0.0;
        double time_1_2 = 0.0;
        double time_1_3 = 0.0;
        
        inline void count_1_begin()
        {
            begin_1 = steady_clock::now();
        }
        inline void count_1_end(int num)
        {
            end_1 = steady_clock::now();
            duration_high elapsed = end_1 - begin_1;
            double time_delay = elapsed.count();
            switch (num)
            {
            case 1:
                time_1_1+= time_delay;
                break;
            case 2:
                time_1_2+= time_delay;
                break;
            case 3:
                time_1_3+= time_delay;
                break;
            
            default:
                break;
            }
        }


        clock_high begin_2, end_2;
        double time_2 = 0.0;

        inline void count_2_begin()
        {
            begin_2 = steady_clock::now();
        }
        inline void count_2_end()
        {
            end_2 = steady_clock::now();
            duration_high elapsed = end_2 - begin_2;
            double time_delay = elapsed.count();
            time_2 += time_delay;
        }

        /// status
        std::string status = "cex found";

        inline void stop_everything(){
            count_whole_end();
            auto now = steady_clock::now();
            main_solver_original_end_ = now;
            main_solver_convergence_calls_end_ = now;
            status = "timeout";
            if(is_counting_convergence)
            {
                count_main_solver_convergence_time_end(0);
                --num_main_solver_convergnece_calls_;
            }
            if(is_counting_main_original)
            {
                count_main_solver_original_time_end(true,0);
                --num_main_solver_original_calls_total;
                --num_main_solver_original_calls_success;
            }
            if(is_counting_imply)
            {
                --count_imply;
                count_imply_end();
            }
        }

        void print() 
        {
            std::string uc_len_original = num_main_solver_original_calls_failed ? to_string(float(sum_main_solver_original_uc_length_) / num_main_solver_original_calls_failed) : "0";
            std::string uc_len_conv = num_main_solver_convergnece_calls_? to_string(float(sum_main_solver_convergence_uc_length_) / num_main_solver_convergnece_calls_) : "0";
            std::string imply_decision_str;
            switch (imply_decision)
            {
            case -1:
                imply_decision_str = "\"not decided\"";
                break;
            case 0:
                imply_decision_str = "\"Solver\"";
                break;
            case 1:
                imply_decision_str = "\"Manual\"";
                break;
            
            default:
                break;
            }

            std::cout<<"{"<<std::endl;

            std::cout << "      \"Status\": \"" << status << "\"," << std::endl;
            std::cout << "      \"Original main solver SAT Calls\": {" << std::endl;
            std::cout << "      \t\"Total Time\": " << time_main_solver_original_calls_ / 1000.0 << "," << std::endl;
            std::cout << "      \t\"Total Count\": " << num_main_solver_original_calls_total << "," << std::endl;
            std::cout << "      \t\"Success\": " << num_main_solver_original_calls_success << "," << std::endl;
            std::cout << "      \t\"Failed\": " << num_main_solver_original_calls_failed << "," << std::endl;
            std::cout << "      \t\"Tried before\": " << num_tried_before << "," << std::endl;

            std::cout << "      \"UC length avergage\": " << uc_len_original << "" << std::endl;
            std::cout << "      }," << std::endl;

            std::cout << "      \"Convergence main solver SAT Calls\": {" << std::endl;
            std::cout << "      \t\"Total Time\": " << time_main_solver_convergence_calls_ / 1000.0 << "," << std::endl;
            std::cout << "      \t\"Total Count\": " << num_main_solver_convergnece_calls_ << "," << std::endl;
            std::cout << "      \t\"shorter UC\": " << num_main_solver_convergence_shorter << "," << std::endl;
            std::cout << "      \t\"UC length avergage\": " << uc_len_conv << "" << std::endl;
            std::cout << "      }," << std::endl;

            std::cout << "      \"Implication\": {" << std::endl;
            std::cout << "      \t\"Decision\": " << imply_decision_str << "," << std::endl;
            std::cout << "      \t\"Solver\": " << time_imply_dec_sol << "," << std::endl;
            std::cout << "      \t\"Manual\": " << time_imply_dec_man << "," << std::endl;

            std::cout << "      \t\"Total Time\": " << time_imply / 1000.0 << "," << std::endl;
            std::cout << "      \t\"Group Count\": " << count_imply << std::endl;
            std::cout << "      }," << std::endl;
            if (!solverWin.empty() || !manualWin.empty())
            {
                std::cout << "      \"Winning History\": " << std::endl;
                serialize_winner();
            }
            std::cout << "      \"nUC\": " << nUC << "," << std::endl;
            std::cout << "      \"nSubUC\": " << nSubUC << "," << std::endl;
            std::cout << "      \"SubUCRate\": \"" << float(nSubUC) / float(nUC) << "\"," << std::endl;
            std::cout << "      \"Rounds of iteration\": " << num_rounds << "," << std::endl;
            std::cout << "      \"Counts of try_by\": " << num_try_by << "," << std::endl;
            std::cout << "      \"Propagation\": {" << std::endl;
            std::cout << "      \t\"Amount\": " << ntime_prop << "," << std::endl;
            std::cout << "      \t\"Succeed\": " << ntime_prop_succ << "," << std::endl;
            std::cout << "      \t\"Total Time\": " << time_prop / 1000.0 << "," << std::endl;
            std::cout << "      \"UC distribution\": {" << std::endl;
            std::cout << "      \t\"[0,3)\": " << num_prop_uc_0_3 << "," << std::endl;
            std::cout << "      \t\"[3,5)\": " << num_prop_uc_3_5 << "," << std::endl;
            std::cout << "      \t\"[5,10)\": " << num_prop_uc_5_10 << "," << std::endl;
            std::cout << "      \t\"[10,50)\": " << num_prop_uc_10_50 << "," << std::endl;
            std::cout << "      \t\"[50,100)\": " << num_prop_uc_50_100 << "," << std::endl;
            std::cout << "      \t\"[100,inf)\": " << num_prop_uc_100_inf << "," << std::endl;
            std::cout << "      }," << std::endl;
            std::cout << "      }," << std::endl;
            std::cout << "      \"UC distribution\": {" << std::endl;
            std::cout << "      \t\"[0,3)\": " << num_uc_0_3 << "," << std::endl;
            std::cout << "      \t\"[3,5)\": " << num_uc_3_5 << "," << std::endl;
            std::cout << "      \t\"[5,10)\": " << num_uc_5_10 << "," << std::endl;
            std::cout << "      \t\"[10,50)\": " << num_uc_10_50 << "," << std::endl;
            std::cout << "      \t\"[50,100)\": " << num_uc_50_100 << "," << std::endl;
            std::cout << "      \t\"[100,inf)\": " << num_uc_100_inf << "," << std::endl;
            std::cout << "      }," << std::endl;

            /// for test uses.
            if(time_1_1> 0 )
                std::cout << "      \"T1\": "     << time_1_1/ 1000.0 <<","<<std::endl;
            if(time_1_2> 0 )
                std::cout << "      \"T2\": "     << time_1_2/ 1000.0 <<","<<std::endl;
            if(time_1_3> 0 )
                std::cout << "      \"T3\": "     << time_1_3/ 1000.0 <<","<<std::endl;
            std::cout << "      \"Global Time\": "     << time_global/ 1000.0 <<""<<std::endl;


            std::cout<<"}"<<std::endl;

        }


};



}

#endif
