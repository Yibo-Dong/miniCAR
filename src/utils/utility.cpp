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
	Update Date: October 30, 2017
	Utility functions
*/

#include "utility.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
using namespace std;

namespace car {

bool absIncr (int i, int j)
{
	return abs (i) < abs(j);
}

bool imply(const std::vector<int>& v1, const std::vector<int>& v2)
{
    std::unordered_set<int> marks;
    for(auto lit: v1)
        marks.insert(lit);
    for(auto lit: v2)
        if(marks.find(lit) == marks.end())
            return false;
    marks.clear();
    return true;
}

Cube negate(const Cube& cu)
{
    Cube res;
    res.reserve(cu.size());
    for(auto &i : cu)
        res.emplace_back(-i);
    return res;
}

}
