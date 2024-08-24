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

#include "definition.h"
#include "carsolver.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
using namespace std;
using namespace Glucose;

namespace car
{
 	
	/**
	 * @brief int -> SAT lit
	 * 
	 * @param id 
	 * @return Lit 
	 */
 	Lit CARSolver::SAT_lit (int id) 
 	{
 		assert (id != 0);
        int var = abs(id)-1;
        int varsNeeded = var - nVars() + 1; // how many needed
        while (varsNeeded > 0) {
            newVar();
            --varsNeeded;
        }
		return ( (id > 0) ? mkLit(var) : ~mkLit(var) );
 	}
 	
	/**
	 * @brief SAT lit -> int
	 * 
	 * @param l 
	 * @return int 
	 */
 	int CARSolver::lit_id (Lit l) const
	{
		if (sign(l)) 
            return -(var(l) + 1);
		else 
            return var(l) + 1;
	}
 	
	
	/**
	 * @brief get the assumption in _assumption
	 * 
	 * @return std::vector<int> 
	 */
	std::vector<int> CARSolver::get_assumption() const
	{
		std::vector<int> res;
		res.reserve(assumptions.size());
		for(int i = 0 ; i < assumptions.size();++i)
		{
			res.push_back(lit_id(assumptions[i]));
		}
		return res;
	}

	/**
	 * @brief return the model from SAT solver when it provides SAT
	 * 
	 * @return std::vector<int> : the result from Solver.
	 */
	std::vector<int> CARSolver::get_model () const
	{
		std::vector<int> res;
		res.resize (nVars (), 0);
   		for(int i = 0; i < nVars (); i ++)
   		{
     		if (model[i] == l_True)
       			res[i] = i+1;
     		else
       			res[i] = -(i+1);
   		}
   		return res;
	}

	/**
	 * 
	 * @brief return the UC from SAT solver when it provides UNSAT
	 * 
	 * @return std::vector<int> 
	 */
 	std::vector<int> CARSolver::get_uc () const
 	{
 		std::vector<int> reason;
		reason.resize(conflict.size(),0);
		// 
 		for(int k = 0; k < conflict.size(); k++) 
 		{
        	Lit l = conflict[k];
        	reason[k] = -lit_id (l);
    	}
    	return reason;
  	}
		
	/**
	 * @brief return the UC without bad from SAT solver(in particular, from `conflict`) when it provides UNSAT
	 * @todo TODO: use std::transform() may help improve the efficiency by parallel.
	 * @return std::vector<int> 
	 */
	std::vector<int> CARSolver::get_uc_no_bad (int bad) const
 	{
 		std::vector<int> reason;
 		for(int k = 0; k < conflict.size(); k++) 
 		{
        	Lit l = conflict[k];
			int id = -lit_id (l);
			if(id!=bad)
	        	reason.push_back (id);
    	}
    	return reason;
  	}

    std::vector<int> CARSolver::get_uc_rand()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        // if std::vector:
        // std::shuffle(vec.begin(), vec.end(), gen); 
        for(int i = assumptions.size() - 1; i > 0; --i) {
            // starting from 1, because 0 is the flag.
            std::uniform_int_distribution<int> dis(1, i);
            int j = dis(gen);
            std::swap(assumptions[i], assumptions[j]);
        }
        solve_();
        return get_uc();
    }

    std::vector<int> CARSolver::get_uc_fp_rev()
    {
        /**
         * @brief previous assumption : <flag> <taken> <rest>
         * while <taken> == <uc> \Union <considered>
         */

        // <flag>. This is the flag of the target frame. It's easy.
        // Lit flag = assumptions[0];
        int sz = assumptions.size();

        for(int i = 1; i <= sz/2; ++i)
        {
            // swap!
            int tmp = assumptions[i].x;
            assumptions[i].x = assumptions[sz-i].x;
            assumptions[sz-i].x = tmp;
        }
        
        solve_();
        return get_uc();
    }

    std::vector<int> CARSolver::get_uc_fp()
    {
        // {not_used, used} for later
        // FIXME: It's too heavy here, as we tested before.
        // problem is: the assumption may have duplicate literals, which makes complex to figure out {not_used}, because the rest literals may have already appeared in the front part.
        assert(false && "Don't use me");
        return {};
    }

    /**
     * @brief Get another UC.
     * @pre Prior UC has been retrieved already.
     * @return Cube 
     */
    std::vector<int> CARSolver::get_uc_another(int option, int nth)
    {
        switch (option)
        {
        case 2: // LO_Fixpoint, 
        {
            return get_uc_fp();
        }
        
        case 1: // LO_Random
        {
            return get_uc_rand();
        }

        case 0: // LO_Classic
        default:
        {
            if(nth == 1)
            {
                return get_uc_fp_rev();
            }
            else{
                return get_uc_rand();
            }
            break;
        }
        }
    }
	 	
	/**
	 * @brief for each literal in the cube, make it a seperate clause 
	 * 
	 * @param cu 
	 */
 	void CARSolver::add_cube (const std::vector<int>& cu)
 	{
 	    for(size_t i = 0; i < cu.size (); i ++)
 	        add_clause (cu[i]);
 	}

	/**
	 * @brief negate each literal in the cube and make it one clause, representing its negation.
	 * 
	 * @param cu 
	 */
	void CARSolver::add_cube_negate (const std::vector<int>&cu)
	{
		add_clause (negateCube(cu));
	}

	/**
	 * @brief put the vector into SAT Solver. (as a clause)
	 * 
	 * @param v 
	 */
	void CARSolver::add_clause_internal (const std::vector<int>& v)
 	{
 		vec<Lit> lits(v.size());
		int index = 0;
		for(int id : v)
			lits[index++] = SAT_lit(id);
 		bool res = addClause (lits);
		assert(res && "Warning: Adding clause does not success\n");
 	}

	/**
	 * @brief helper function, print all the clauses in the Solver.
	 * 
	 */
 	void CARSolver::print_clauses(ostream & out)
	{
		out << "clauses in SAT solver: \n";
		for(int i = 0; i < clauses.size (); i ++)
		{
			Glucose::Clause& c = ca[clauses[i]];
			for(int j = 0; j < c.size (); j ++)
				out << lit_id (c[j]) << " ";
			out << "0 " << endl;
		}
	}
	
	/**
	 * @brief helper function, print the assumption in the Solver.
	 * 
	 */
	void CARSolver::print_assumption (ostream & out)
	{
	    out << "assumptions in SAT solver: \n";
		if (!assumptions.size())
			out<<" Empty ";
	    for(int i = 0; i < assumptions.size (); i ++)
	        out << lit_id (assumptions[i]) << " ";
	    out << endl;
	}
	
}
