#ifndef DEBUG_PRINTING
#define DEBUG_PRINTING

// #define PRINT_SAT
// #define PRINT_PROOF

// print prior map in `print_evidence()`
// #define PRINT_PRIOR

#define LOG(msg)


#ifdef PQUEUE
#define CONTAINER priority_queue<item,std::vector<item>, CompareItem>
#else
#define CONTAINER stack<item>
#endif

#ifdef PRINT_SIMPLE_SAT
#define PRINTIF_SIMPLE_SAT()                                                            \
	static int sat_cnt = 0;                                                           \
	sat_cnt++;                                                                        \
	cout << "Q" << sat_cnt << ":\tS = {";                                             \
	for (auto i : s->s())                                                             \
		cout << i << ", ";                                                            \
	cout << "},\tid = " << s->id << ",\tl = " << level << ",\tres = " << res << endl; \
	if (res == false)                                                                 \
	{                                                                                 \
		Cube uc = solver->get_conflict(!backward_first);                              \
		cout << "uc:";                                                                \
		for (int i : uc)                                                              \
			cout << i << ", ";                                                        \
		cout << endl;                                                                 \
	}
#else
#define PRINTIF_SIMPLE_SAT() 
#endif

#ifdef PRINT_PROOF
#define PRINTIF_PROOF()                                                \
    cout << "[ Print O sequence ]" << endl;                          \
    cout << "O size : " << o->size() << endl;                        \
    int level = 0;                                                   \
    for (auto &frame : *o)                                           \
    {                                                                \
        cout << "Level " << level++ << endl;                         \
        for (auto &cube : frame)                                     \
        {                                                            \
            cout << "(";                                             \
            for (auto l : cube)                                      \
            {                                                        \
                cout << l << ", ";                                   \
            }                                                        \
            cout << "), ";                                           \
        }                                                            \
        cout << endl;                                                \
    }                                                                \
    cout << "[ End of Print O sequence ]" << endl                    \
         << endl;                                                    \
    main_solver->print_last_3_clauses();                             \
    cout << "------------[End of One round]----------------" << endl \
         << endl;
#else
#define PRINTIF_PROOF()

#endif // PRINT_PROOF

#ifdef PRINT_PRIOR
#define PRINTIF_PRIOR()
cout << "UF" << endl;
print_U_sequnece(Uf, cout);
cout << "UB" << endl;
print_U_sequnece(Ub, cout);
cout << "F" << endl;
for (auto pair : prior_in_trail_f)
{
    cout << pair.first->id << " <- " << (pair.second ? to_string(pair.second->id) : "-1") << endl;
    cout << pair.first->id << " : " << endl;
    cout << "L:" << pair.first->latches() << endl;
    cout << "I:" << pair.first->inputs() << endl;
    cout << "Last:" << pair.first->last_inputs() << endl;
}

cout << "B" << endl;
for (auto pair : prior_in_trail_b)
{
    cout << pair.first->id << " <- " << (pair.second ? to_string(pair.second->id) : "-1") << endl;
    cout << pair.first->id << " : " << endl;
    cout << "L:" << pair.first->latches() << endl;
    cout << "I:" << pair.first->inputs() << endl;
    cout << "Last:" << pair.first->last_inputs() << endl;
}
#else
#define PRINTIF_PRIOR()
#endif // PRINT_PRIOR

#ifdef PRINT_QUERY
PRINTIF_QUERY() print_sat_query(solver, s, O, level, res);
#else
#define  PRINTIF_QUERY()   
#endif // PRINT_QUERY

#ifdef PRINT_UNREACHABLE
#define PRINTIF_UNREACHABLE()
    cout<<"This unreachable state: ";
    for(int i:otherCEX()->s())
    cout<<i<<", ";
    cout<<endl;
    print_sat_query(bi_main_solver,s,o,dst,false,cout);
#else
#define PRINTIF_UNREACHABLE()
#endif // PRINT_UNREACHABLE

// #define COUNT_POSSIBLE_REDUNDANT
// #define PRINT_SIMPLE_SAT
// #define PRINT_UC1
// #define PRINT_PRINT_SAT_IMMEDIATE_SAT
// #define PRINT_BACKTRACK
// #define PRINT_MAINSOLVER_INIT
// #define PRINT_ASSUMPTION
// #define PRINT_UNREACHABLE

#endif // !DEBUG_PRINTING


