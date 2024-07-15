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

class Statistics 
{
    public:
        Statistics () {}
        ~Statistics () {}


        // main solver's information:
        
        // how many are original ones
        int num_main_solver_original_calls_total = 0;
        // how many succeeded
        int num_main_solver_original_calls_success = 0;
        // how many failed
        int num_main_solver_original_calls_failed = 0;
        // time spent
        double time_main_solver_original_calls_ = 0.0;
        // uc's length in total.
        int sum_main_solver_original_uc_length_ = 0;
        clock_high main_solver_original_begin_, main_solver_original_end_;

        bool is_counting_main_original;

        inline void count_main_solver_original_time_start ()
        {
#ifdef STAT
            main_solver_original_begin_ = steady_clock::now();
#endif
        is_counting_main_original = true;
        }

        int last_original_uc_sz = 0;

        inline void count_main_solver_original_time_end (bool res, int uc_length)
        {
#ifdef STAT
            main_solver_original_end_ = steady_clock::now();
            duration_high elapsed = main_solver_original_end_ - main_solver_original_begin_;
            double time_delay = elapsed.count();

	        time_main_solver_original_calls_ += time_delay;
#endif
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
            }
            
        }

        // convergence information:
        int num_main_solver_convergnece_calls_ = 0;
        int num_main_solver_convergence_shorter = 0;
        double time_main_solver_convergence_calls_ = 0.0;
        clock_high main_solver_convergence_calls_begin_, main_solver_convergence_calls_end_;
        int sum_main_solver_convergence_uc_length_ = 0;
        bool is_counting_convergence;

        inline void count_main_solver_convergence_time_start()
        {
#ifdef STAT
            main_solver_convergence_calls_begin_ = steady_clock::now();
#endif
            is_counting_convergence = true;
        }
        inline void count_main_solver_convergence_time_end(int uc_length)
        {
#ifdef STAT
            main_solver_convergence_calls_end_ = steady_clock::now();
            duration_high elapsed = main_solver_convergence_calls_end_ - main_solver_convergence_calls_begin_;
            double time_delay = elapsed.count();
            time_main_solver_convergence_calls_ += time_delay;
#endif
            is_counting_convergence = false;
            num_main_solver_convergnece_calls_ ++;
            
            sum_main_solver_convergence_uc_length_ += uc_length;
            if(uc_length > 0 && uc_length < last_original_uc_sz)
                ++num_main_solver_convergence_shorter;
        }

        // search progress information:
        int num_rounds = 0;
        int num_try_by = 0;
        clock_high time_enter_new_try_by;
        inline void count_enter_new_ronud(){
            ++num_rounds;
        }
        inline void count_enter_new_try_by(){
            ++num_try_by;
#ifdef STAT
            time_enter_new_try_by = steady_clock::now();
#endif
        }

        // global time
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

        // tried before
        int num_tried_before = 0;
        inline void count_tried_before() { 
#ifdef STAT
            clock_high tried_before_end = steady_clock::now();
            duration_high elapsed = tried_before_end - time_enter_new_try_by;
            time_imply += elapsed.count();
#endif
            ++num_tried_before; 
        }

        // time for implication calculation
        clock_high imply_begin_, imply_end_;
        int count_imply = 0;
        double time_imply = 0.0;
        bool is_counting_imply;
        inline void count_imply_begin()
        {
#ifdef STAT
            imply_begin_ = steady_clock::now();
#endif
        is_counting_imply = true;
        }
        inline void count_imply_end()
        {
#ifdef STAT
            imply_end_ = steady_clock::now();
            duration_high elapsed = imply_end_ - imply_begin_;
            double time_delay = elapsed.count();
            time_imply += time_delay;
#endif
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
            // not the decision time
            if(imply_dec_cnter!=0)
                return -1;
            // sample too few            
            if(std::max(time_imply_dec_sol,time_imply_dec_man*10) < 100)
            {
                imply_dec_cnter = 100;
                return -1;
            }// decide.
            
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

        // statistics about subsumption.
        int nUC = 0, nSubUC = 0;
        inline void recordUC(bool sub)
        {
            nUC++;
            nSubUC += sub;
        }

        // helper for debug.
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

        // status
        std::string status = "cex found";

        inline void stop_everything(){
            count_whole_end();
#ifdef STAT
            auto now = steady_clock::now();
            main_solver_original_end_ = now;
            main_solver_convergence_calls_end_ = now;
#endif
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

            std::cout << "      \"Status\": \""<<status <<"\","<<std::endl;
            std::cout << "      \"Original main solver SAT Calls\": {" <<std::endl;
            std::cout << "      \"Total Time\": "    << time_main_solver_original_calls_ / 1000.0 <<","<<std::endl;
            std::cout << "      \"Total Count\": "    << num_main_solver_original_calls_total <<","<<std::endl;            
            std::cout << "      \"Success\": "  << num_main_solver_original_calls_success <<","<<std::endl;        
            std::cout << "      \"Failed\": "   << num_main_solver_original_calls_failed <<","<<std::endl;
            std::cout << "      \"Tried before\": "   << num_tried_before <<","<<std::endl;
            
            std::cout << "      \"UC length avergage\": "   << uc_len_original <<""<<std::endl; 
            std::cout << "      },"<<std::endl;

            std::cout << "      \"Convergence main solver SAT Calls\": {" <<std::endl;
            std::cout << "      \"Total Time\": "    << time_main_solver_convergence_calls_ / 1000.0 <<","<<std::endl;
            std::cout << "      \"Total Count\": "    << num_main_solver_convergnece_calls_ <<","<<std::endl; 
            std::cout << "      \"shorter UC\": "   << num_main_solver_convergence_shorter <<","<<std::endl; 
            std::cout << "      \"UC length avergage\": "   << uc_len_conv <<""<<std::endl; 
            std::cout << "      },"<<std::endl;

            std::cout << "      \"Implication\": {" <<std::endl;
            std::cout << "      \"Decision\": "    << imply_decision_str <<","<<std::endl;
            std::cout << "      \"Solver\": "    << time_imply_dec_sol <<","<<std::endl;
            std::cout << "      \"Manual\": "    << time_imply_dec_man <<","<<std::endl;
            
            std::cout << "      \"Total Time\": "    << time_imply / 1000.0 <<","<<std::endl;
            std::cout << "      \"Group Count\": "    << count_imply <<std::endl; 
            std::cout << "      },"<<std::endl;
            if(!solverWin.empty() || !manualWin.empty())
            {
                std::cout << "      \"Winning History\": "    <<std::endl; 
                serialize_winner();
            }
            std::cout << "      \"nUC\": "  << nUC <<","<<std::endl;
            std::cout << "      \"nSubUC\": "     << nSubUC <<","<<std::endl;
            std::cout << "      \"SubUCRate\": "     << float(nSubUC) / float(nUC) <<","<<std::endl;
            std::cout << "      \"Rounds of iteration\": "  << num_rounds <<","<<std::endl;
            std::cout << "      \"Counts of try_by\": "     << num_try_by <<","<<std::endl;
            std::cout << "      \"Global Time\": "     << time_global/ 1000.0 <<""<<std::endl;
            
            std::cout<<"}"<<std::endl;

        }


};



}

#endif
