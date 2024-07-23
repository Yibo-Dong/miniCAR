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
	Author: Jianwen Li
	Update Date: October 4, 2017
	Invariant Solver in CAR
*/

#ifndef INV_SOLVER_H
#define INV_SOLVER_H

#include "definition.h"
#include "carsolver.h"
#include "statistics.h"
#include <vector>

namespace car 
{
	extern Statistics CARStats;
	class InvSolver : public CARSolver
	{
		bool verbose_;
		public:
			InvSolver (const Problem* m, bool verbose=false) : verbose_(verbose),id_aiger_max_ (const_cast<Problem*>(m)->max_id ())
			{
				model_ = const_cast<Problem*> (m);
			    int end =  model_->outputs_start ();
			    for (int i = 0; i < end ; i ++)
                    add_clause (model_->element (i));
			}
			~InvSolver () {}
		
			/**
			 * @brief call solve_assumption().
			 * 
			 * @return true 
			 * @return false 
			 */
			inline bool solve_with_assumption ()
			{
				if (verbose_)
					std::cout << "InvSolver::";
				bool res = solve_assumption ();
				return res;
			}
			
			/**
			 * @brief 
			 * Add the negation of this frame into the solver
			 * @param frame 
			 */
			inline void add_constraint_or (const OFrame &frame, int level)
			{
				// there should be no prior flag for this frame. 
				assert(or_flag.find(level) == or_flag.end());
				
				// get a new frame flag
				or_flag[level] = new_var();
				// std::cout<<"add new or frame flag for level "<<level <<" : "<<or_flag[level]<<std::endl;
				
				update_assumption_for_constraint(or_flag[level]);

				// uc_flags: at least one uc should be true.
				std::vector<int>& uc_flags = or_ucflags[level];
 				for (int i = 0; i < frame.size (); i ++)
 				{
					// flag for this uc.
 					int clause_flag = new_var ();
 					uc_flags.push_back (clause_flag);
 					for (int j = 0; j < frame[i].size (); j ++)
 					{
 						int id = frame[i][j];
 						add_clause (-clause_flag, id);
 					}
 				}
				std::vector<int> frame_clause = uc_flags;
				frame_clause.push_back(-or_flag[level]);
 				add_clause (frame_clause);
			}
			
			/**
			 * @brief 
			 * Add the real states into this solver.
			 * @param frame 
			 */
			inline void add_constraint_and (const OFrame &frame, int level)
			{
				// flag = 1 : Clauses in this frame is enabled.
				// The real states this frame Represent
				//	~uc1, ~uc2, ...
				// flag = 0 : Not enabled
				int frame_flag = inv_and_flag_of(level);
 				for (int i = 0; i < frame.size (); i ++)
 				{
 					add_uc_and(frame[i],level);
 				}
				#ifdef PRINT_INV
				std::cout<<"[flag] new and-frame flag: "<<frame_flag<<std::endl;
 				#endif
				update_assumption_for_constraint (frame_flag);
			}

			inline void add_uc(const Cube&uc, int level, int state_level)
			{
				// cout<<"update inv: level = "<<level <<" uc = ";
				// for(int i:uc)
				// cout<<i<<", ";
				// cout<<endl;

				// Oj. 
				if(level == state_level+1)
					add_uc_and(uc,level);
				else
				// Ok, k < j.
					add_uc_or(uc,level);
			}

			inline void add_uc_and(const Cube &uc, int level){
				int frame_flag = inv_and_flag_of(level);
				// v := ~l1 \/ ~l2 \/ ... \/ ~flag
				std::vector<int> v;
				for (int j = 0; j < uc.size (); ++j)
				{
					int id = uc[j];
					v.push_back (-id);
				}
				v.push_back (-frame_flag);
				add_clause (v);
			}

			/**
			 * @brief Add this uc to target level.
			 * 1. turn off old or_flag.
			 * 2. add new uc clause, create new or_flag, turn it on.
			 * 
			 */
			inline void add_uc_or(const Cube &uc, int level){
				// prior assumption
				std::vector<int> old_flags = get_assumption();
				
				// invalidate old frame flag.
				int old_flag = or_flag[level];
				assert(old_flag!=0);
				for(int i = 0; i < old_flags.size(); ++i)
				{
					if(old_flags[i] == old_flag)
					{
						// negate it
						old_flags[i] = -old_flag;
						#ifdef PRINT_INV
						std::cout<<"[flag] block old or-frame flag: "<<-or_flag[level]<<std::endl;
						#endif
						break;
					}
				}
				
				// generate a new frame flag
				or_flag[level] = new_var();
				// generate a new clause flag for new uc.
				int clause_flag = new_var ();
				// add clause for this uc.
				for (int j = 0; j < uc.size (); j ++)
				{
					add_clause (-clause_flag, uc[j]);
				}

				// add this uc's flag.
				or_ucflags[level].push_back(clause_flag);

				// clause to be inserted this time
				std::vector<int> frame_clause = or_ucflags[level];
				frame_clause.push_back(-or_flag[level]);
 				add_clause (frame_clause);
				
				#ifdef PRINT_INV
				std::cout<<"[flag] new flag for or-frame: "<<or_flag[level]<<std::endl;
				#endif
				old_flags.push_back(or_flag[level]);
				set_assumptions(old_flags);				
			}
			
			/**
			 * @brief add this flag into assumption.
			 * 
			 * @param frame_flag 
			 */
			inline void update_assumption_for_constraint (const int frame_flag)
			{
				assumptions.push (SAT_lit (frame_flag));
			}
			
			inline void set_assumptions(const std::vector<int>& assumes)
			{
				assumptions.clear();
				for(auto &lit:assumes)
				{
					assumptions.push(SAT_lit(lit));
				}
			}
			
			/**
			 * @brief pop the last lit of assumption, and negate it.
			 * 
			 */
			inline void release_constraint_and (int level)
			{
				// invalidate old and_flag
				std::vector<int> old_assumptions = get_assumption();
				for(int i = 0; i < old_assumptions.size();++i)
				{
					if(old_assumptions[i] == and_flag[level])
					{
						old_assumptions[i] = -and_flag[level];
						break;
					}
				}
				set_assumptions(old_assumptions);
				
				add_clause(-and_flag[level]);
				simplify();
			}
			
			inline int new_var () {return ++id_aiger_max_;}


			// FIXME: merge them with MainSolver
			inline State* get_state(const bool forward)
			{
				Assignment model = get_model();
				Assignment inputs(model.begin(),model.begin() + model_->num_inputs());
				Assignment latches(model.begin() + model_->num_inputs(),model.begin() + model_->num_inputs()+model_->num_latches());
				State* s = new State(inputs,latches);
				return s;
			}

			void shrink_model (Assignment& model)
			{
				Assignment res=model;
				for (int i = model_->num_inputs ()+1; i <= model_->num_inputs () + model_->num_latches (); i ++)
				{
					int p = model_->prime (i);
					assert (p != 0);
					assert (model.size () > abs (p));
					
					int val = model[abs(p)-1];
					if (p == val)
						res[i-1]=(i);
					else
						res[i-1]=(-i);
				}
				model = res;
			}

			void inv_update_U(USequence &U, State *s, int level, State * prior_state_in_trail)
			{
                
				while(U.size() <= level)
					U.push_back({});
				
				U.push_back(s);
				prior[s] = prior_state_in_trail;
			}

			void addUCtoSolver(Cube&uc, OSequence* O ,int level)
			{
				OFrame&frame = (*O)[level];
				//FIXME: add clause from cube.
			}

		public:
			std::map<State *, State*> prior;

		private:
			// we only store the flags for and, because we should never update a prior frame.
			std::map<int,int> and_flag;
			inline int inv_and_flag_of(int level)
			{
				if(and_flag.find(level) == and_flag.end())
				{
					and_flag[level] = new_var();
				}
				return and_flag[level];
			}

			std::map<int, int> or_flag;
			inline int inv_or_flag_of(int level)
			{
				if(or_flag.find(level) == or_flag.end())
				{
					or_flag[level] = new_var();
				}
				return or_flag[level];
			}

			std::map<int, std::vector<int>> or_ucflags;

			std::vector<int> frame_flags;	

		protected:
			Problem* model_;
			int id_aiger_max_;  	//to store the maximum number used in aiger model

	};
}
#endif
