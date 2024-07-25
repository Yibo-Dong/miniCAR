#include "bmcChecker.h"
using namespace std;
using namespace bmc;

void BMCChecker::check()
{
    bool res = false;
    
    while (!res)
    {
        cout << "unroll level is " << solver->unroll_level;
        // try with level = lev
        int target = model_->output(0) + (model_->output(0) > 0 ? 1 : -1) * lev * solver->lits_per_round();
        cout << ", target is " << target;
        res = solver->badcheck(init->get_latches(), target);
        cout << ", res is " << res << endl;
        if (res)
            break;
        ++lev;
        // unroll once ~
        solver->unroll();
    }
    solver->get_states(cex);
}

void BMCChecker::printEvidence(std::ostream& res_file)
{
    // print counter example
    res_file << "1" << endl
             << "b0" << endl;
    for (int i = 0; i < model_->num_latches(); ++i)
        res_file << "0";
    res_file << endl;
    for (auto s : cex)
        res_file << s->get_inputs_str() << endl;
    res_file << "." << endl;
    // end printing counter example

}
