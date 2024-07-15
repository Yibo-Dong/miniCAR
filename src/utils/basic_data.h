#ifndef BASIC_DATA_H
#define BASIC_DATA_H
#include <vector>
namespace car
{
    /* 1-100
	1,1, -1, -1 ,-1 ...
	1,2 √
	*/
	typedef std::vector<int> Assignment;
	// a /\ b /\ c
	typedef std::vector<int> Cube;
	// a \/ b \/ c
	typedef std::vector<int> Clause;
	// ~(c1 \/ c2 \/ c3)  c_i都是uc，不用管初始化（全集？）的问题
	typedef std::vector<Cube> Frame;
	//
	typedef std::vector<Frame> Fsequence;
}
#endif