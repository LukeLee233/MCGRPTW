#include <iostream>
#include <sys/timeb.h>
#include "utils.h"
#include "instance.h"
#include "RNG.h"
#include "local_search.h"
#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>
#include "slice.h"
#include "extraction.h"
#include "swap.h"
#include "invert.h"
#include "two_opt.h"
#include <boost/program_options.hpp>
#include <algorithm>
#include "json.hpp"
#include "config.h"
#include "initilizer.h"
//#include "Memetic.h"

using namespace std;
namespace bpo = boost::program_options;
using json = nlohmann::json;

int main(int argc, char *argv[])
{
    /*-------------------------parse command line----------------------------*/
    bpo::options_description opts("all options");
    bpo::variables_map vm;

    opts.add_options()
        ("help", "MCGRPTW program")
        ("directory,dir", bpo::value<std::string>(&instance_directory), "the instance directory")
        ("config_file,config", bpo::value<string>(&config_file)->default_value(""), "parameter_file")
        ("search_time", bpo::value<int>(&search_time), "total search time")
        ("neighbor_search_mode",
         bpo::value<std::string>(&neighbor_search_mode),
         "select mode in neighbor search")
        ("local_minimum_threshold",
         bpo::value<int>(&local_minimum_threshold),
         "determine whether we arrived at a local optimal")
        ("tabu_step", bpo::value<int>(&tabu_step), "tabu steps in the ascent search")
        ("infeasible_distance_threshold",
         bpo::value<double>(&infeasible_distance_threshold),
         "the parameter determine the distance between infeasible region and feasible region")
        ("evolve_step,es", bpo::value<int>(&iterations), "number of evolving step")
        ("intensive_local_search_portion", bpo::value<double>(&intensive_local_search_portion))
        ("merge_split_portion", bpo::value<double>(&merge_split_portion))
        ("random_seed,rs", bpo::value<int>(&random_seed), "number of evolving step")
        ("neighbor_size,ns", bpo::value<int>(&neighbor_size), "neighbor size of the local search")
        ("qndf_weights,qs", bpo::value<double>(&QNDF_weights)->default_value(0), "weight of solution distance");

    bpo::store(bpo::parse_command_line(argc, argv, opts), vm);
    bpo::notify(vm);

    if (vm.count("help")) {
        cout << opts << endl;
        return 0;
    }

    if (instance_directory.empty()) {
        cerr << "you need to specify a directory!\n";
        return -1;
    }

    if (!config_file.empty()) {
        json j;
        ifstream fin(config_file);
        fin >> j;

        j.at("iterations").get_to(iterations);
        j.at("random_seed").get_to(random_seed);
        j.at("neighbor_size").get_to(neighbor_size);
        j.at("search_time").get_to(search_time);
        j.at("neighbor_search_mode").get_to(neighbor_search_mode);
        j.at("local_minimum_threshold").get_to(local_minimum_threshold);
        j.at("tabu_step").get_to(tabu_step);
        j.at("infeasible_distance_threshold").get_to(infeasible_distance_threshold);
        j.at("intensive_local_search_portion").get_to(intensive_local_search_portion);
        j.at("merge_split_portion").get_to(merge_split_portion);
    }

    iterations = (iterations == -1) ? 1e6 : iterations;
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
    json j;
    j["iterations"] = iterations;
    j["random_seed"] = random_seed;
    j["neighbor_size"] = neighbor_size;
    j["search_time"] = search_time;
    j["neighbor_search_mode"] = neighbor_search_mode;
    j["local_minimum_threshold"] = local_minimum_threshold;
    j["tabu_step"] = tabu_step;
    j["infeasible_distance_threshold"] = infeasible_distance_threshold;
    j["intensive_local_search_portion"] = intensive_local_search_portion;
    j["merge_split_portion"] = merge_split_portion;


    ofstream fout(date_folder + + "/parameters.json");
    fout << setw(4) << j << endl;
    fout.close();
/*----------------------------------------------------------------*/

    vector<string> file_set = read_directory(instance_directory);
    for (auto file_name : file_set) {
        cout << string(2, '\n') << string(24, '-')
             << "Start instance: "
             << file_name
             << string(24, '-') << string(2, '\n')
             << flush;


        /************************************************************************************************/
        /* initialize the instance object */
        RNG rng = RNG();
        MCGRPTW instance(instance_directory + '/' + file_name, rng);

        instance.load_file_info();
        instance.create_neighbor_lists(neighbor_size);
        instance.register_distance("cost", unique_ptr<Distance>(new CostDistance(instance)));
        instance.register_distance("hybrid", unique_ptr<Distance>(new HybridDistance(instance, 0.3)));

#ifdef DEBUG
        log_out.open(date_folder + '/' + file_name + ".log", ios::out);
#endif
        auto search_start_time = system_clock::now();

        instance._rng.change(seed[random_seed % seed_size]);
        instance.best_total_route_length = DBL_MAX;
        instance.best_sol_time = DBL_MAX;

        LocalSearch local_search(instance, tabu_step);
        local_search.initialize_score_matrix(instance);

        Individual initial_solution = CWScanner(instance)();
//        Individual initial_solution = NearestScanner(instance, *instance.distance_look_tbl["cost"])();
//                        Individual initial_solution = nearest_scanning(instance, vector<int>());

        auto neg = get_negative_coding(initial_solution.sequence);
        local_search.unpack_seq(initial_solution.sequence, instance);
        local_search.trace(instance);
        cout << "Begin Local search..." << endl;

        for (auto iter = 1; iter <= iterations
            && get_time_difference(search_start_time, cur_time) < search_time; iter++) {
            cout << "start "<< iter << "th iterations...\n";
            iteration_start_time = system_clock::now();
            local_search.neighbor_search(instance);
            cur_time = system_clock::now();
            cout << "Finish " << iter << "th iterations, spent: "
                 << get_time_difference(iteration_start_time, cur_time) << 's' << endl;
        }

        cout << "Search finish!\n"
             << "Total spend " << get_time_difference(search_start_time, cur_time) <<"s"<< endl;


        Result global_best = Result(instance.best_total_route_length,
                                           instance.best_sol_time,
                                           instance.best_sol_neg);

        if (global_best.cost == DBL_MAX) {
            cerr << "ERROR! Can't find a solution\n";
            return -1;
        }

#ifdef DEBUG
        log_out.close();
#endif

        /***********************************output part********************************************/
        cout << string(2, '\n') << string(24, '*')
             << "Searching Finish!" << string(24, '*') << endl;
        cout << "Searching result:" << string(2, '\n') << flush;

        result_out.open(result_dir + '/' + file_name + ".res", ios::out);
        result_out << setprecision(2) << fixed;

        print(cout, result_out, "Instance: " + file_name);
        print(cout, result_out, "Best Cost: " + to_string(global_best.cost));
        print(cout, result_out, "Best Solution Search Time: " + to_string(global_best.time) + 's');
        print(cout, result_out, "Best Solution: ");

        string buff = "";
        for (int i = 0; i < global_best.solution.size(); i++) {
            if (i == global_best.solution.size() - 1)
                buff += to_string(global_best.solution[i]);
            else
                buff += (to_string(global_best.solution[i]) + "->");
        }
        print(cout, result_out, buff);

        result_out.close();

    }

    return 0;
}