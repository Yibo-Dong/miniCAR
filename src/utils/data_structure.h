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

#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <map> 
#include "statistics.h"
#include "basic_data.h"
#include "model.h"

namespace car
{
	/* 1-100
	1,1, -1, -1 ,-1 ...
	1,2 √
	*/
	extern Statistics CARStats;

	// state
	class State
	{
		public:
			// construct a special state called negp.
			State(bool isnegp) : is_negp(isnegp), id(-1) { assert(isnegp); negp_state=this; }

			State(const Assignment &latches) : s_(latches), id(get_id(latches)) {}
			State(const Assignment &inputs , const Assignment &latches): s_(latches), inputs_(inputs), id(get_id(latches)) {}

			~State() {}

			// s_ -> cu
			// 当前状态是否在cu代表的状态中
			bool imply(const Cube &cu) const;

			// 取交
			Cube intersect(const Cube &cu);
			Cube intersect_no_dupl(const Cube &cu);
			Cube intersect(const State*);


			inline void set_inputs(const Assignment &st) { inputs_ = st; }

			inline void print() { std::cout << latches() << std::endl; }

			inline const Assignment& s() const { return s_; }
			inline Assignment &inputs_vec() { return inputs_; }
			std::string inputs() const;
			std::string latches() const;

			inline int size() const { return s_.size(); }
			inline int element(int i) const { return s_.at(i); }

		public:
			const unsigned id;
			static const State* negp_state;
			bool is_negp = false;
		
			static void set_num_inputs_and_latches (const int n1, const int n2) 
			{
				num_inputs_ = n1;
				num_latches_ = n2;
			}
			
		public:
			static const Model* model_;
			static const aiger* aig_;
		private:
			Assignment s_;
			std::vector<int> inputs_;

			static std::map<Assignment, int> ids;
			inline unsigned get_id(const Assignment& latches){
				return ++next_id_;
			}
		public:
			static int num_inputs_;
			static int num_latches_;
			static unsigned next_id_;

		public:

			// calculate what is already known about next latches according to present latch values.
			Cube next_latches();
	};

	typedef std::vector<std::vector<State *>> Bsequence;
	using Osequence = std::vector<Frame>;
	using Usequence = std::vector<State *>;
	using item = std::tuple<car::State *, int, std::size_t>;
	struct CompareItem {
		bool operator()(const item &a, const item &b) const {
			return int(std::get<2>(a)) >= int(std::get<2>(b)); // 使用元组的第三个成员进行比较
		}
	};

}
#endif
