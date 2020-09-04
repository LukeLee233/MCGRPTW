#include <iostream>
#include <sys/timeb.h>
#include "file.h"
#include "utils.h"
#include "MCGRP.h"
#include "RNG.h"
#include "NeighborSearch.h"
#include "Memetic.h"
#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>
#include "ConstructPolicy.h"
#include "SingleInsert.h"
#include "DoubleInsert.h"
#include "Slice.h"
#include "Extraction.h"
#include "Swap.h"
#include "Invert.h"
#include "TwoOpt.h"
#include <boost/program_options.hpp>
#include <algorithm>
#include <numeric>
#include "json.hpp"

using namespace std;
namespace bpo = boost::program_options;
using json = nlohmann::json;

struct timeb epoch_start_time;

int main(int argc, char *argv[])
{
    /*-------------------------parse command line----------------------------*/
    string instance_directory;
    string config_file;
    bool freeze_config;

    int pool_size;
    int evolve_steps;
    int phase_number;
    int random_seed;
    int neighbor_size;
    double QNDF_weights;

    bpo::options_description opts("all options");
    bpo::variables_map vm;

    opts.add_options()
        ("help", "MCGRP program")
        ("directory,dir", bpo::value<std::string>(&instance_directory), "the instance directory")
        ("config_file,config", bpo::value<string>(&config_file)->default_value(""), "parameter_file")
        ("freeze_config", bpo::value<bool>(&freeze_config)->default_value(false), "whether to record the parameters")
        ("pool_size,ps", bpo::value<int>(&pool_size)->default_value(5), "size of population")
        ("evolve_step,es", bpo::value<int>(&evolve_steps)->default_value(5), "number of evolving step")
        ("phase_number,np", bpo::value<int>(&phase_number)->default_value(1))
        ("random_seed,rs", bpo::value<int>(&random_seed)->default_value(0), "number of evolving step")
        ("neighbor_size,ns", bpo::value<int>(&neighbor_size)->default_value(25), "neighbor size of the local search")
        ("qndf_weights,qs", bpo::value<double>(&QNDF_weights)->default_value(0.6), "weight of solution distance");

    bpo::store(bpo::parse_command_line(argc, argv, opts), vm);
    bpo::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return 0;
    }

    if (instance_directory.empty()) {
        std::cerr << "you need to specify a directory!\n";
        return -1;
    }

    if (!config_file.empty()) {
        json j;
        ifstream fin(config_file);
        fin >> j;

        j.at("pool_size").get_to(pool_size);
        j.at("evolve_steps").get_to(evolve_steps);
        j.at("phase_number").get_to(phase_number);
        j.at("random_seed").get_to(random_seed);
        j.at("neighbor_size").get_to(neighbor_size);
        j.at("QNDF_weights").get_to(QNDF_weights);
    }

    if (random_seed < 0 || random_seed > 100) {
        std::cerr << "Incorrect random seed base address!\n";
        return -1;
    }

    if (freeze_config) {
        json j;
        j["pool_size"] = pool_size;
        j["evolve_steps"] = evolve_steps;
        j["phase_number"] = phase_number;
        j["random_seed"] = random_seed;
        j["neighbor_size"] = neighbor_size;
        j["QNDF_weights"] = QNDF_weights;

        config_file = config_file.empty() ? "demo.json" : config_file;
        ofstream fout(config_file);
        fout << setw(4) << j << endl;
    }

    /*----------------------------------------------------------------*/

    /*----------------------Create search info log----------------------*/
    //create log floder
    string log_dir = instance_directory + "/log";
    boost::filesystem::path log_dir_path(log_dir);
    if (!(boost::filesystem::exists(log_dir_path))) {
        cout << "Create log directory..." << endl;
        if (boost::filesystem::create_directory(log_dir_path)) {
            cout << "Log Directory created successfully!\n\n";
        }
    }
    else {
        cout << "Log directory already exists!\n\n";
    }

    // create result folder
    boost::posix_time::ptime time_stamp = boost::posix_time::second_clock::local_time();
    stringstream buffer;
    string date_;
    buffer << time_stamp;
    date_ = buffer.str();
    string date_folder = instance_directory + "/log/" + date_;
    boost::filesystem::path date_dir_path(date_folder);
    if (!(boost::filesystem::exists(date_dir_path))) {
        cout << "Create" << date_ << " directory..." << endl;
        if (boost::filesystem::create_directory(date_dir_path)) {
            cout << "Directory created successfully!\n\n";
        }
    }
    else {
        cout << "Result directory already exists!\n\n";
    }

    // create result folder
    string result_dir = date_folder + "/result";
    boost::filesystem::path result_dir_path(result_dir);
    if (!(boost::filesystem::exists(result_dir_path))) {
        cout << "Create result directory..." << endl;
        if (boost::filesystem::create_directory(result_dir_path)) {
            cout << "Directory created successfully!\n\n";
        }
    }
    else {
        cout << "Result directory already exists!\n\n";
    }

/*----------------------Create info recorder----------------------*/
//  record searching parameters
    log_out.open(date_folder + "/parameters.txt", ios::out);
    print(cout, log_out, "Parameters settings:");
    print(cout, log_out, "directory " + instance_directory);
    print(cout, log_out, "base address of seeds: " + to_string(random_seed));
    print(cout, log_out, "search times: " + to_string(phase_number));
    print(cout, log_out, "evolution times: " + to_string(evolve_steps));
    print(cout, log_out, "pool size: " + to_string(pool_size));
    print(cout, log_out, "neighbor_size: " + to_string(neighbor_size));
    print(cout, log_out, "QNDF weights: " + to_string(QNDF_weights));
    log_out.close();
/*----------------------------------------------------------------*/


    vector<string> file_set = read_directory(instance_directory);
    for (auto file_name : file_set) {
        cout << string(2,'\n') << string (24,'-')
             << "Start instance: "
             << file_name
             << string (24,'-') << string(2,'\n')
             << flush;

        /* global info */
        double bestobj = numeric_limits<decltype(bestobj)>::max();
        double best_solution_time = numeric_limits<decltype(best_solution_time)>::max();
        int best_epoch = -1;
        vector<int> best_buffer;
        vector<double> fitnessVec;
        vector<double> BesttimeVec;
        vector<double> SearchtimeVec;

        /************************************************************************************************/
        /* initialize the instance object */
        instance_num_information instance_info;
        GetTasksNum(instance_directory + '/' + file_name, instance_info);

        RNG rng = RNG();
        MCGRP Mixed_Instance(instance_info, rng);

        Mixed_Instance.load_file_info(instance_directory + '/' + file_name, instance_info);
        Mixed_Instance.create_neighbor_lists(neighbor_size);


        log_out.open(date_folder + '/' + file_name + ".log", ios::out);

        for (int start_seed = random_seed; start_seed < random_seed + phase_number; start_seed++) {
            cout << string (24,'-') << "Start "
                 << start_seed - random_seed + 1
                 << "th times search" << string (24,'-') << endl;

            Mixed_Instance._rng.change(seed[start_seed]);
            Mixed_Instance.best_total_route_length = numeric_limits<double>::max();
            Mixed_Instance.best_sol_time = numeric_limits<double>::max();

            HighSpeedNeighBorSearch NBS(Mixed_Instance);

            ftime(&epoch_start_time);        //针对每轮搜索时间进行测试


            /*----------------------------------------------------------*/
            HighSpeedMemetic MA(NBS, pool_size, evolve_steps, QNDF_weights);
            cout << "Begin Memetic search...\n";
            ftime(&epoch_start_time);
            struct timeb cur_time;
            MA.memetic_search(Mixed_Instance);
            ftime(&cur_time);
            cout << "Finish " << start_seed - random_seed << "th search, spent: "
                 << (cur_time.time - epoch_start_time.time)
                     + ((cur_time.millitm - epoch_start_time.millitm) * 1.0 / 1000) << 's' << endl;

            MA.reset();
            log_out.close();
            cout << "Finished!\n\n\n";
            /*----------------------------------------------------------*/

            /* solution record */
            ftime(&cur_time);
            SearchtimeVec.push_back((cur_time.time - epoch_start_time.time) +
                ((cur_time.millitm - epoch_start_time.millitm) * 1.0 / 1000));

            /* record the best info of each epoch */
            fitnessVec.push_back(Mixed_Instance.best_total_route_length);
            BesttimeVec.push_back(Mixed_Instance.best_sol_time);

            cout << "Finish " << start_seed - random_seed + 1 << "th times\n";

            /* record the best solution during the whole searching process*/
            if (Mixed_Instance.best_total_route_length < bestobj) {
                best_solution_time = Mixed_Instance.best_sol_time;
                bestobj = Mixed_Instance.best_total_route_length;
                best_buffer = Mixed_Instance.best_sol_buff;
                best_epoch = start_seed - random_seed + 1;
            }

        }

        /***********************************output part********************************************/
        if (bestobj == numeric_limits<decltype(bestobj)>::max()) {
            cerr << "ERROR! Can't find a solution\n";
            abort();
        }
        double total_cost = accumulate(fitnessVec.begin(), fitnessVec.end(), 0);
        double best_total_time = accumulate(BesttimeVec.begin(), BesttimeVec.end(), 0);
        double total_time = accumulate(SearchtimeVec.begin(), SearchtimeVec.end(), 0);
        double average_cost = total_cost / fitnessVec.size();
        double average_best_time = best_total_time / BesttimeVec.size();
        double average_search_time = total_time / SearchtimeVec.size();

        cout << string(2,'\n') << string(24,'*')
             << "Searching Finish!" << string(24,'*') << endl;
        cout << "Searching result:"<<string(2,'\n') << flush;


        result_out.open(result_dir + '/' + file_name + ".res", ios::out);
        result_out << setprecision(2) << fixed;

        print(cout, result_out, "Instance: " + file_name);
        print(cout, result_out, "Best Cost: " + to_string(bestobj));
        print(cout, result_out, "Best Cost Epoch: " + to_string(best_epoch));
        print(cout, result_out, "Best Solution Search Time: " + to_string(best_solution_time) + 's');
        print(cout, result_out, "Best Solution: ");

        string buff = "";
        for (int i = 0; i < best_buffer.size(); i++) {
            if (i == best_buffer.size() - 1)
                buff += to_string(best_buffer[i]);
            else
                buff += (to_string(best_buffer[i]) + "->");
        }
        print(cout, result_out, buff);

        print(cout, result_out, "\n\nSearch time: " + to_string(fitnessVec.size()));


        print(cout, result_out, "\n\nAverage cost: " + to_string(average_cost));
        print(cout, result_out, "Average best solution search time: " + to_string(average_best_time) + 's');
        print(cout, result_out, "Average search time: " + to_string(average_search_time) + 's');
        result_out.close();
    }

    return 0;
}