//
// Created by luke on 2020/9/4.
//

#ifndef CONFIG_H
#define CONFIG_H
#include <string>

extern std::string instance_directory;
extern std::string config_file;

extern int search_time;

extern int iterations;
extern double merge_split_portion;
extern int random_seed;
extern int neighbor_size;
extern double intensive_local_search_portion;
extern double QNDF_weights;

extern std::string neighbor_search_mode;

// Determine whether we arrive at a local minimum
extern int local_minimum_threshold;

extern int tabu_step;

// Infeasible Local Search
extern double infeasible_distance_threshold;


#endif //CONFIG_H
