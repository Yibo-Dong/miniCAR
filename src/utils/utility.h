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

#ifndef UTILITY_H
#define UTILITY_H

#include "basic_data.h"
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>


namespace car{


//elements in v1, v2 are in order
//check whether v2 is contained in v1 
bool imply (const std::vector<int>& v1, const std::vector<int>& v2);

bool absIncr (int i, int j);

// my section

Cube negate(const Cube& cu);

// Cube minus(const Cube& c1, const Cube& c2);

// Cube intersect(const Cube& c1, const Cube& c2);

template <typename T>
void shuffle(std::vector<T>& vec) {
    #ifdef RANDSEED
        int seed = RANDSEED;
    #else
        int seed = 1;
    #endif
    std::mt19937 g(seed);
    std::shuffle(vec.begin(), vec.end(), g);
}

}

#endif


