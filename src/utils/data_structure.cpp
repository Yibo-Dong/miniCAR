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
	Update Date: September 20, 2017
	Data structures in CAR
*/

 
 #include <vector>
 #include <stdlib.h>
 #include <string.h>
 #include <assert.h>
 #include <map> 
 #include <unordered_map>
 #include <unordered_set>
 #include <algorithm>
 #include "utility.h"
 #include "data_structure.h"

 using namespace std;
 
 namespace car
 {
 	const State* State::negp_state;
 	
    /**
     * @brief Whether this state is already blocked by this cube.
     * @pre need state(this) to follow the literal order.
     * @param cu 
     * @return true 
     * @return false 
     */
 	bool State::imply (const Cube& cu) const
	{
		for (int i = cu.size() - 1 ; i >= 0; --i)
		{
			int index = abs(cu[i]) - num_inputs_ - 1;
			assert (index >= 0);
			if (s_[index] ^ cu[i])
			{
				return false;
			}
		}
		return true;
	}
	
    /**
     * @brief calculate intersection of two cubes.
     * @pre need state(this) to follow the literal order.
     * @post follow the order in cu.
     * @param cu 
     * @return true 
     * @return false 
     */
	Cube State::intersect (const Cube& cu) 
	{
		Cube res;
		for (int i = 0; i < cu.size (); i ++)
		{
			int index = abs(cu[i]) - num_inputs_ - 1;
			assert (index >= 0);
			if (s_[index] == cu[i])
				res.push_back (cu[i]);
		}
		return res;
	}

	Cube State::intersect_no_dupl (const Cube& cu) 
	{
		Cube res;
		std::unordered_set<int> mark;
		for (int i = 0; i < cu.size (); i ++)
		{
			int index = abs(cu[i]) - num_inputs_ - 1;
			assert (index >= 0);
			if (s_[index] == cu[i] && mark.find(cu[i]) == mark.end())
			{
				res.push_back (cu[i]);
				mark.insert(cu[i]);
			}
		}
		mark.clear();
		return res;
	}

	Cube State::intersect (const State *s) 
	{
		Cube res;
		Cube latches = s->s();

		for (int i = 0; i < latches.size (); i ++)
		{
			if (s_[i] == latches[i])
				res.push_back( s_[i] > 0 ? (i +  num_inputs_ + 1) : -(i +  num_inputs_ + 1));
		}
		return res;
	}
 	

	string State::inputs() const
	{
		string res = "";
		for (int i = 0; i < inputs_.size(); i++)
			res += (inputs_.at(i) > 0) ? "1" : "0";
		return res;
	}

	string State::latches() const
	{
		string res = "";
		int j = 0;
		for (int i = 0; i < num_latches_; i++)
		{
			if (j == s_.size())
				res += "x";
			else if (num_inputs_ + i + 1 < abs(s_.at(j)))
				res += "x";
			else
			{
				res += (s_.at(j) > 0) ? "1" : "0";
				j++;
			}
		}
		return res;
	}

	int State::num_inputs_ = 0;
 	int State::num_latches_ = 0;

	unsigned State::next_id_ = 0;

	Cube State::next_latches()
	{
		// if it's negp, nothing is already known.
		if(is_negp)
			return {};

		Cube known;
		std::unordered_map<int,int> value;
		// value for 0/1 is constant.
		value[0] = 0;
		
		// value for inputs.
		// if(inputs_.empty())
		// {
			for (int i = 0; i < aig_->num_inputs; ++i)
			{
				value[i+1] = -1;
			}
		// }else{
		// 	for(int i = 0; i < aig_->num_inputs;++i)
		// 	{
		// 		value[i+1] = (inputs_[i] > 0) ? 1 : 0;
		// 	}
		// }

		// value for latches
		for(int i = 0; i < s_.size();++i)
		{
			int offset = abs(s_[i] -1 - int(aig_->num_inputs));
			value[1+aig_->num_inputs +i ] = (s_[i]>0) ? 1 : 0;
		}
		
		// auto getLitValue = [&value](int x){
		// 	assert(value.find(x) != value.end() || value.find(-x) != value.end());
		// 	if(x<0){return value[x] == -1 ? -1 : 1 - value[x];}
		// 	else{return value[x];}
		// };

		auto getAigValue = [&value](int x){
			// assert(value.find(x>>1) != value.end() || value.find(1 + (x>>1)) != value.end());
			if(x & 1){ assert(value.find(x>>1) != value.end());  return value[x>>1] == -1 ? -1 : 1 - value[x>>1];}
			else{assert(value.find(x>>1) != value.end()); return value[x>>1];}
		};
		
		// value for and gates
		for (int i = 0; i < aig_->num_ands; i ++)
		{
			aiger_and& aa = aig_->ands[i];
			int rhs0_value = getAigValue(aa.rhs0);
			int rhs1_value = getAigValue(aa.rhs1);
			int lhs_value;

			if(rhs0_value == 0 || rhs1_value == 0)
				lhs_value = 0;
			// taint will spread. mark as not determined
			else if(rhs0_value == -1 || rhs1_value ==-1)
			{
				lhs_value = -1;
				if(aa.rhs0 == aa.rhs1)
					lhs_value = 1;
				else if(aa.rhs0 == aa.rhs1 + 1 && aa.rhs0&1 )
					lhs_value = 0;
				// if(lhs_value != -1)
				// 	cerr<<"?= "<<aa.rhs0<<endl;
			}
			else if (rhs0_value == 1 && rhs1_value == 1)
				lhs_value = 1;
			// else
			// {
			// 	cerr<<"strange value,"<<aa.rhs0 <<":" <<rhs0_value<<", "<<aa.rhs1<<": "<<rhs1_value<<endl;
			// }
			value[aa.lhs >> 1] = lhs_value;

			// cerr<<aa.lhs<<" = "<<aa.rhs0 <<" + "<<aa.rhs1<<endl;
			// cerr<<getAigValue(aa.lhs)<<" = "<<getAigValue(aa.rhs0)<< " + "<<getAigValue(aa.rhs1) <<endl;
		}

		// calculate necessaries.
		for(int i = 0; i< aig_->num_latches; ++i)
		{
			int offset = 1+aig_->num_inputs;
			int next_value = getAigValue(aig_->latches[i].next);
			if(next_value!=-1)
			{
				// cerr<<(next_value>0);
				known.push_back(next_value==0? -(i + offset) : i+ offset);
			}
			// else{
			// 	cerr<<"x";
			// }
		}
		// if(id == 5)
		// {
		// 	std::vector<int> values(aig_->maxvar,0);
			
		// 	cout<<"values: ";
		// 	for(auto &p:value)
		// 	{
		// 		values[p.first] = p.second;
		// 	}	
		// 	for(int i = 0; i< values.size();++i)
		// 		cout<<i<<" : "<<values[i]<<", ";
		// 	cout<<endl;
		// 	cout<<"known :";
		// 	for(int i:known)
		// 	cout<<i<<", ";
		// 	cout<<endl;
		// }
		// cout<<"known:";
		// for(int i:known)
		// 	cout<<i<<", ";
		// cout<<endl;


		
		return known;
	}

	std::map<Assignment, int> State::ids;
}
 		
