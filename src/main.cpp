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

#include "bmcChecker.h"
#include "carChecker.h"
#include "statistics.h"
#include "definition.h"
#include "implysolver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <assert.h>
using namespace std;
using namespace car;

namespace car
{
    // global variables section.
    Statistics CARStats;
    Checker *chk;
    Problem *model;
}

/**
 * @brief do printing even if timeout.
 *
 * @param sig_num : the signal number.
 */
void signal_handler(int sig_num)
{
    CARStats.stop_everything();
    CARStats.print();
    exit(0);
}

void print_usage()
{
    printf("Usage: simplecar <-f|-b|-p|-e|-v|-h> <aiger file> <output directory>\n");
    printf("       -f          forward checking (Default = backward checking)\n");
    printf("       -b          backward checking \n");
    printf("       -p          enable propagation (Default = off)\n");
    printf("       -e          print witness (Default = off)\n");
    printf("       -h          print help information\n");
}

// cut filename from path.
// e.g. xxxxxx/xyz123.aig ==> xyz123
string get_file_name(string &s)
{
    size_t start_pos = s.find_last_of("/");
    if (start_pos == std::string::npos)
        start_pos = 0;
    else
        start_pos += 1;
    size_t end_pos = s.find(".aig", start_pos);
    if(end_pos == std::string::npos)
        end_pos = s.find(".aag", start_pos);
    assert(end_pos != std::string::npos);

    return s.substr(start_pos, end_pos - start_pos);
}

OPTIONS parse_args(int argc, char **argv)
{
    OPTIONS opt;
    for(int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--vb") == 0)
        {
            opt.inv_incomplete = true;
            opt.enable_rotate = true;
            opt.inter_cnt = 1;
            opt.raw_uc = true;
            opt.forward = false;
            opt.impMethod = 9;
            opt.LOStrategy = 3;
        }
        else if (strcmp(argv[i], "-f") == 0)
            opt.forward = true;
        else if (strcmp(argv[i], "-b") == 0)
            opt.forward = false;
        else if (strcmp(argv[i], "--propMode") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.propMode = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--propParam") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.propParam = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--partial") == 0)
            opt.partial = true;
        else if (strcmp(argv[i], "-h") == 0)
            print_usage();
        else if (strcmp(argv[i], "--rotate") == 0)
            opt.enable_rotate = true;
        else if (strcmp(argv[i], "--inter") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.inter_cnt = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--restart") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.time_limit_to_restart = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--rem") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.rememOption = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--bmc") == 0)
        {
            opt.bmc = true;
        }
        else if (strcmp(argv[i], "--imp") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.impMethod = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--incomplete") == 0)
        {
            opt.inv_incomplete = true;
        }
        else if (strcmp(argv[i], "--raw") == 0)
        {
            opt.raw_uc = true;
        }
        else if (strcmp(argv[i], "--convMode") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.convMode = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--convParam") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.convParam = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--convAmount") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.convAmount = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--order") == 0)
        {
            assert(i + 1 < argc);
            ++i;
            opt.LOStrategy = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--sub") == 0)
        {
            opt.subStat = true;
        }
        else if (strcmp(argv[i], "--ms") == 0)
        {
            opt.multi_solver = true;
        }
        else if (strcmp(argv[i], "--simp") == 0)
        {
            opt.simplifyCNF = true;
        }
        else if (opt.inputPath.empty())
        {
            opt.inputPath = string(argv[i]);
        }
        else if (opt.outputPath.empty())
        {
            opt.outputPath = string(argv[i]);
            if (opt.outputPath.at(opt.outputPath.size() - 1) != '/')
                opt.outputPath += "/";
        }
        else
        {
            cerr << "unrecognized option: " << argv[i] << endl;
            exit(1);
        }
    }
    if (opt.inputPath.empty() || opt.outputPath.empty())
    {
        cerr << "missing input or output" << endl;
        exit(1);
    }
    return opt;
}

void check_aiger(int argc, char **argv)
{
    OPTIONS opt = parse_args(argc, argv);

    std::string filename = get_file_name(opt.inputPath);
    std::string stdout_filename = opt.outputPath + filename + ".log";
    std::string res_file_name = opt.outputPath + filename + ".res";
    // redirect to log file.
    auto __fs = freopen(stdout_filename.c_str(), "w", stdout);
    ofstream res_file;
    res_file.open(res_file_name.c_str());

    // get aiger object
    aiger *aig = aiger_init();
    aiger_open_and_read_from_file(aig, opt.inputPath.c_str());
    const char *err = aiger_error(aig);
    if (err)
    {
        printf("read aiger file error!\n");
        exit(2);
    }
    if (!aiger_is_reencoded(aig))
        aiger_reencode(aig);

    car::model = new Problem(aig, opt.unrollPrime);
    State::set_num_inputs_and_latches(model->num_inputs(), model->num_latches());

    // assume that there is only one output needs to be checked in each aiger model,
    // which is consistent with the HWMCC format
    assert(model->num_outputs() == 1);

    // if BMC, then use BMC checker to check.
    if (opt.bmc)
    {
        auto bchker = new bmc::BMCChecker(model);
        bchker->check();
        bchker->printEvidence(res_file);
        return;
    }

    CARStats.count_whole_begin();
    chk = new Checker(model, opt, res_file, nullptr);
    RESEnum res = chk->check();

    if (opt.time_limit_to_restart <= 0)
    {
        // if Restart NOT enabled, use one checker to check.
        CARStats.count_whole_end();
        delete chk;
    }
    else
    {
        // if RESTART enabled:
        // TODO: verify the result of restart
        assert(opt.convMode >= 0 && "restart with multiple UCs, mode should be given");
        auto clear_delay = chk;
        while (RES_RESTART == res)
        {
            ++opt.convParam;
            ImplySolver::reset_all();
            CARStats.reset_imply_cnter();
            // TODO: check whether we need to clean sth else?

            // check, with information of last check avaiable.
            chk = new Checker(model, opt, res_file, clear_delay);
            res = chk->check();
            delete clear_delay;
            clear_delay = chk;
        }
        CARStats.count_whole_end();
        delete chk;
    }

    // cleaning work
    aiger_reset(aig);
    delete model;
    res_file.close();

    CARStats.print();
    return;
}

int main(int argc, char **argv)
{
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
    check_aiger(argc, argv);
    return 0;
}
