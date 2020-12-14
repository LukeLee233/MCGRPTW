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
#include "slice.h"
#include "extraction.h"
#include "swap.h"
#include "invert.h"
#include "TwoOpt.h"
#include <boost/program_options.hpp>
#include <algorithm>
#include <numeric>
#include "json.hpp"
#include "config.h"
#include "ConstructPolicy.h"

using namespace std;
namespace bpo = boost::program_options;
using json = nlohmann::json;

struct timeb phase_start_time;

struct timeb cur_time;

int main(int argc, char *argv[])
{
    /*-------------------------parse command line----------------------------*/
    bpo::options_description opts("all options");
    bpo::variables_map vm;

    opts.add_options()
        ("help", "MCGRP program")
        ("directory,dir", bpo::value<std::string>(&instance_directory), "the instance directory")
        ("config_file,config", bpo::value<string>(&config_file)->default_value(""), "parameter_file")
        ("freeze_config", bpo::value<bool>(&freeze_config)->default_value(false), "whether to record the parameters")
        ("search_time", bpo::value<int>(&search_time)->default_value(60), "total search time")
        ("neighbor_search_mode",
         bpo::value<std::string>(&neighbor_search_mode)->default_value("random"),
         "select mode in neighbor search")
        ("significant_search",
         bpo::value<int>(&significant_search)->default_value(5),
         "if valid move steps doesn't less than this, descent search will be terminated")
        ("local_ratio",
         bpo::value<double>(&local_ratio)->default_value(0.8),
         "if the percentage of equal move larger than this number, descent search will be terminated")
        ("local_minimum_threshold",
         bpo::value<int>(&local_minimum_threshold)->default_value(40),
         "determine whether we arrived at a local optimal")
        ("tabu_step", bpo::value<int>(&tabu_step)->default_value(500), "tabu steps in the ascent search")
        ("infeasible_distance",
         bpo::value<double>(&infeasible_distance)->default_value(0.3),
         "the parameter determine the distance between infeasible region and feasible region")
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
        j.at("search_time").get_to(search_time);
        j.at("QNDF_weights").get_to(QNDF_weights);
        j.at("neighbor_search_mode").get_to(neighbor_search_mode);
        j.at("significant_search").get_to(significant_search);
        j.at("local_ratio").get_to(local_ratio);
        j.at("local_minimum_threshold").get_to(local_minimum_threshold);
        j.at("tabu_step").get_to(tabu_step);
        j.at("infeasible_distance").get_to(infeasible_distance);
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

    phase_number = (phase_number == -1) ? 1e6 : phase_number;
    evolve_steps = (evolve_steps == -1) ? 1e6 : evolve_steps;
    /*----------------------------------------------------------------*/

    /*----------------------Create search info log----------------------*/
    //create log folder
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

/*----------------------Create parameter logger----------------------*/
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


    vector<string> file_set{"TWe10_t.dat"};
    for (auto file_name : file_set) {
        cout << string(2, '\n') << string(24, '-')
             << "Start instance: "
             << file_name
             << string(24, '-') << string(2, '\n')
             << flush;

        /* global info */
        double bestobj = numeric_limits<decltype(bestobj)>::max();
        double best_solution_time = numeric_limits<decltype(best_solution_time)>::max();
        int best_phase = -1;
        vector<int> best_buffer;
        vector<double> FitnessVec;
        vector<double> BestTimeVec;
        vector<double> SearchTimeVec;
        vector<vector<int>> SolutionVec;

        /************************************************************************************************/
        /* initialize the instance object */
        InstanceNumInfo instance_info;
        GetTasksNum(instance_directory + '/' + file_name, instance_info);

        RNG rng = RNG();
        MCGRP Mixed_Instance(instance_info, rng);

        Mixed_Instance.load_file_info(instance_directory + '/' + file_name, instance_info);
        Mixed_Instance.create_neighbor_lists(neighbor_size);
        Mixed_Instance.register_distance("cost", unique_ptr<Distance>(new CostDistance(Mixed_Instance)));
        Mixed_Instance.register_distance("hybrid", unique_ptr<Distance>(new HybridDistance(Mixed_Instance, 0.3)));

#ifdef DEBUG
        log_out.open(date_folder + '/' + file_name + ".log", ios::out);
#endif
        struct timeb search_start_time;
        ftime(&search_start_time);
        ftime(&cur_time);
        for (int start_seed = random_seed; start_seed < random_seed + phase_number
            && get_time_difference(search_start_time, cur_time) < search_time; start_seed++) {
            cout << string(24, '-') << "Start "
                 << start_seed - random_seed + 1
                 << "th times search" << string(24, '-') << endl;

            Mixed_Instance._rng.change(seed[start_seed % seed_size]);
            Mixed_Instance.best_total_route_length = numeric_limits<double>::max();
            Mixed_Instance.best_sol_time = numeric_limits<double>::max();

            HighSpeedNeighBorSearch NBS(Mixed_Instance, tabu_step);

            /*----------------------------------------------------------*/
            cout << "Begin Optimal check..." << endl;
            vector<vector<int>> node_seq{
                {1, 2, 4, 69, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 20, 76, 20, 19, 21, 22, 75, 23, 26, 25, 24, 53, 51, 49, 47,
                 48, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 21, 22, 24, 53, 24, 25, 29, 25, 26, 23, 31, 32, 34, 32, 35,
                 41, 42, 43, 44, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 41, 35, 32, 33, 36, 33, 37, 39, 40, 41, 42, 43, 44, 59, 11, 10, 9, 8, 7,
                 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 21, 22, 24, 25, 24, 53, 24, 25, 26, 25, 24, 53, 51, 49, 47,
                 48, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 60, 61, 60, 67, 56, 42, 43, 44, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 41, 40, 39, 40, 41, 35, 39, 37, 39, 40, 39, 37, 39, 40, 41, 42, 43, 44,
                 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 11, 12, 76, 20, 18, 72, 73, 72, 18, 20, 18, 72, 71, 73, 72, 18, 20, 76, 12, 11, 10, 9,
                 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 21, 22, 24, 53, 52, 50, 19, 20, 76, 12, 11, 10, 9, 8, 7, 5, 4,
                 2, 1},
                {1, 2, 4, 69, 59, 11, 12, 16, 13, 14, 15, 17, 15, 14, 13, 16, 12, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 41, 35, 32, 33, 37, 39, 40, 41, 42, 43, 44, 59, 11, 10, 9, 8, 7, 5, 4, 2,
                 1},
                {1, 2, 4, 5, 4, 2, 3, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 69, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 11, 12, 16, 13, 14, 15, 77, 13, 14, 15, 77, 76, 12, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 60, 62, 63, 65, 63, 64, 63, 62, 67, 56, 42, 43, 44, 59, 11, 10, 9,
                 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 11, 12, 76, 20, 18, 15, 77, 76, 12, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 60, 62, 67, 56, 42, 43, 42, 41, 42, 43, 44, 59, 11, 10, 9, 8, 7,
                 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 21, 22, 75, 23, 31, 30, 31, 32, 35, 41, 42, 43, 44, 59, 11,
                 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 60, 67, 56, 55, 56, 42, 43, 44, 59, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 59, 44, 46, 47, 49, 50, 19, 18, 72, 73, 72, 71, 70, 71, 73, 72,
                 18, 20, 76, 12, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 21, 22, 75, 23, 26, 28, 29, 25, 26, 27, 26, 28, 29, 25, 24,
                 53, 51, 49, 47, 48, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 60, 67, 56, 42, 43, 44, 46, 45, 44, 59, 11, 10, 9, 8, 7, 5, 4, 2,
                 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 19, 21, 51, 49, 47, 48, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 46, 47, 49, 50, 52, 50, 19, 21, 51, 49, 47, 46, 47, 49, 50, 52, 54, 52, 50, 49,
                 47, 48, 11, 10, 9, 8, 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 44, 43, 42, 57, 58, 60, 62, 66, 68, 66, 62, 60, 67, 56, 42, 43, 44, 59, 11, 10, 9, 8,
                 7, 5, 4, 2, 1},
                {1, 2, 4, 69, 59, 11, 12, 16, 13, 14, 15, 77, 76, 12, 11, 10, 9, 8, 7, 5, 4, 2, 1}
            };
            int total_length = 0;
            for (const vector<int> &route_seq: node_seq) {
                for (int i = 1; i < route_seq.size(); i++) {
                    total_length += Mixed_Instance.min_cost[route_seq[i - 1]][route_seq[i]];
                }
            }
            cout << total_length << endl;

            vector<vector<vector<int>>> tasks = {
                {{11,11},{8,7},{7,5},{2,2}},
                {{50,50},{76,76},{76,20},{22,22},{75,23}},
                {{53,53},{29,29},{23,31},{34,34}},
                {{32,32},{36,33},{41,41}},
                {{4,69},{44,44},{59,59}},
                {{22,24},{25,25},{25,24},{24,53},{25,26}},
                {{61,60},{60,67}},
                {{39,40},{35,39},{37,37},{40,40},{39,37}},
                {{73,73},{20,20},{18,72},{71,71}},
                {{21,22},{53,52},{50,19}},
                {{13,14},{17,17},{15,14}},
                {{35,32},{33,37}},
                {{5,5},{5,4},{3,2},{2,1}},
                {{57,57},{10,10}},
                {{57,58},{58,69}},
                {{12,12},{16,16},{77,13},{15,77}},
                {{63,65},{64,63},{63,62}},
                {{59,11},{20,18},{18,15}},
                {{62,67},{67,56},{43,42},{42,41},{43,44}},
                {{30,30},{31,32}},
                {{55,56},{44,59}},
                {{58,59},{19,18},{73,72},{72,72},{71,70}},
                {{23,26},{26,28},{26,27},{28,29}},
                {{58,60},{56,42},{46,45}},
                {{44,46},{19,21},{51,49},{47,48},{48,48}},
                {{52,52},{21,51},{49,49},{46,47},{54,54}},
                {{68,66},{62,62},{60,60}},
                {{12,16},{77,76},{9,8}}
            };

            vector<int> task_seq = Mixed_Instance.parse_task(tasks);
//            vector<int> task_seq = {-5,76,98,-86,122,59,87,121,-79,112,46,60,81,-66,65,64,-106,134,45,43,42,-127,130,-62,34,132,-54,84,117,119,17,118,-133,110,35,-12,80,13,47,-50,16,-78,88,89,-115,15,-74,73,37,-63,101,-91,94,-100,77,111,97,-124,123,44,69,131,68,-92,105,-67,129,128,-95,120,36,75,-114,83,-125,113,82,116,-102,96,90,85,-61,19,21,-108,99,-93,20,22,126,-104,39,2,103,-107,4,109};

            auto sorted_solution = sort_solution(task_seq);
            for (const auto task : sorted_solution) {
                if (task < 0) cout << endl;
                cout << task << "->";
            }

            NBS.unpack_seq(get_delimiter_coding(task_seq), Mixed_Instance);
            cout << "unpack solution: " << NBS.get_cur_cost() << endl;

            NBS.parse_solution(Mixed_Instance);


            cout << "Finish " << start_seed - random_seed << "th search, spent: "
                 << get_time_difference(phase_start_time, cur_time) << 's' << endl;

            /*----------------------------------------------------------*/

            /* solution record */
            ftime(&cur_time);
            SearchTimeVec.push_back(get_time_difference(phase_start_time, cur_time));

            /* record the best info of each epoch */
            FitnessVec.push_back(Mixed_Instance.best_total_route_length);
            BestTimeVec.push_back(Mixed_Instance.best_sol_time);
            SolutionVec.push_back(Mixed_Instance.best_sol_neg);

            cout << "Finish " << start_seed - random_seed + 1 << "th times\n";

            /* record the best solution during the whole searching process*/
            if (Mixed_Instance.best_total_route_length < bestobj) {
                best_solution_time = Mixed_Instance.best_sol_time;
                bestobj = Mixed_Instance.best_total_route_length;
                best_buffer = Mixed_Instance.best_sol_neg;
                best_phase = start_seed - random_seed + 1;
            }

        }

#ifdef DEBUG
        log_out.close();
#endif

        /***********************************output part********************************************/
        if (bestobj == numeric_limits<decltype(bestobj)>::max()) {
            cerr << "ERROR! Can't find a solution\n";
            abort();
        }
        double total_cost = accumulate(FitnessVec.begin(), FitnessVec.end(), 0);
        double best_total_time = accumulate(BestTimeVec.begin(), BestTimeVec.end(), 0);
        double total_time = accumulate(SearchTimeVec.begin(), SearchTimeVec.end(), 0);
        double average_cost = total_cost / FitnessVec.size();
        double average_best_time = best_total_time / BestTimeVec.size();
        double average_search_time = total_time / SearchTimeVec.size();

        cout << string(2, '\n') << string(24, '*')
             << "Searching Finish!" << string(24, '*') << endl;
        cout << "Searching result:" << string(2, '\n') << flush;


        result_out.open(result_dir + '/' + file_name + ".res", ios::out);
        result_out << setprecision(2) << fixed;

        print(cout, result_out, "Instance: " + file_name);
        print(cout, result_out, "Best Cost: " + to_string(bestobj));
        print(cout, result_out, "Best Cost Epoch: " + to_string(best_phase));
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

        print(cout, result_out, "\n\nSearch time: " + to_string(FitnessVec.size()));


        print(cout, result_out, "\n\nAverage cost: " + to_string(average_cost));
        print(cout, result_out, "Average best solution search time: " + to_string(average_best_time) + 's');
        print(cout, result_out, "Average search time: " + to_string(average_search_time) + 's');
        result_out.close();
    }

    return 0;
}