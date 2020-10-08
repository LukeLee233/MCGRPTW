//
// Created by luke on 2019/11/29.
//

#include <algorithm>
#include "ConstructPolicy.h"
#include "RNG.h"
#include <unordered_map>

using namespace std;

Individual nearest_scanning(const MCGRP &mcgrp, vector<int> unserved_task_set)
{
    int load;
    int trial;
    int min_dist;
    int drive_time; // the earliest time of a vehicle begins to serve the task

    unordered_set<int> unserved_task_id_set;

    std::vector<int> candidate_task_set;

    std::vector<int> nearest_task_set;

    int current_tail_task;
    int chosen_task;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // 统计未服务任务
    int serve_task_num = 0;
    if(!unserved_task_set.empty()){
        serve_task_num = unserved_task_set.size();
        unserved_task_id_set = unordered_set<int>(unserved_task_set.begin(),unserved_task_set.end());
    }
    else{
        serve_task_num = mcgrp.req_arc_num + mcgrp.req_node_num + mcgrp.req_edge_num;
        for (int i = 1; i <= mcgrp.actual_task_num; i++)
            unserved_task_id_set.insert(i);
    }

    load = 0;
    trial = 0;
    drive_time = 0;

    vector<int> solution;
    solution.push_back(DUMMY);

    while (trial < serve_task_num) {
        current_tail_task = solution.back();

        //Find all the tasks that satisfy the capacity constraint
        candidate_task_set.clear();
        for (auto unserved_task : unserved_task_id_set) {
            if (mcgrp.inst_tasks[unserved_task].demand <= mcgrp.capacity - load
                && mcgrp.cal_arrive_time(current_tail_task, unserved_task, drive_time, true)
                    <= mcgrp.inst_tasks[unserved_task].time_window.second) {
                candidate_task_set.push_back(unserved_task);
            }
        }

        if (candidate_task_set.empty()) {
            solution.push_back(DUMMY);
            load = 0;
            drive_time = 0;
            continue;
        }

        min_dist = MAX(min_dist);

        //Find the nearest task from the current candidate task set
        for (auto candidate_task : candidate_task_set) {
            if (mcgrp.min_cost[mcgrp.inst_tasks[current_tail_task].tail_node][mcgrp.inst_tasks[candidate_task]
                .head_node]
                < min_dist) {
                min_dist =
                    mcgrp.min_cost[mcgrp.inst_tasks[current_tail_task].tail_node][mcgrp.inst_tasks[candidate_task]
                        .head_node];
                nearest_task_set.clear();
                nearest_task_set.push_back(candidate_task);
            }
            else if (
                mcgrp.min_cost[mcgrp.inst_tasks[current_tail_task].tail_node][mcgrp.inst_tasks[candidate_task]
                    .head_node]
                    == min_dist) {
                nearest_task_set.push_back(candidate_task);
            }
        }

        //If multiple tasks both satisfy the capacity constraint and are closest, randomly choose one
        int k = (int) mcgrp._rng.Randint(0, nearest_task_set.size() - 1);
        chosen_task = nearest_task_set[k];

        load += mcgrp.inst_tasks[chosen_task].demand;
        drive_time = mcgrp.cal_arrive_time(current_tail_task, chosen_task, drive_time, true);

        trial++;
        solution.push_back(chosen_task);

        unserved_task_id_set.erase(chosen_task);

        if (mcgrp.inst_tasks[chosen_task].inverse != ARC_NO_INVERSE
            && mcgrp.inst_tasks[chosen_task].inverse != NODE_NO_INVERSE) {
            int inverse_task = mcgrp.inst_tasks[chosen_task].inverse;
            if(unserved_task_id_set.find(inverse_task) != unserved_task_id_set.end())
                unserved_task_id_set.erase(inverse_task);
        }
    }

    if (solution.back() != DUMMY) {
        solution.push_back(DUMMY);
    }

    return mcgrp.parse_delimiter_seq(solution);
}


const int max_merge_set_num = 20;

void merge_split(NeighBorSearch &ns, const MCGRP &mcgrp, const int merge_size, const int pseudo_capacity)
{
    My_Assert(!ns.delimiter_coding_sol.empty(), "Neighbor search doesn't hold a sequence!");

    if (ns.routes.size() < merge_size) {
        return;
    }


    COMB combination_manager(ns.routes.size(), merge_size);

    int count = 0;
    vector<int> chosen_routes_id;

    double best_choice = MAX(best_choice);
    vector<int> best_buffer;

    while (count < max_merge_set_num && combination_manager.get_combinations(chosen_routes_id)) {
        count++;

        vector<vector<int>> task_routes(ns.routes.size());

        //generate tasks distribution in routes
        int cursor = -1;
        for (auto task : ns.negative_coding_sol) {
            if (task < 0) {
                cursor++;
                task_routes[cursor].push_back(-task);
            }
            else if (task > 0) {
                task_routes[cursor].push_back(task);
            }
            else {
                My_Assert(false, "dummy task cannot occur in negative coding soluiton!");
            }
        }

        //extract chosen tasks
        vector<int> candidate_tasks;
        for (auto route_id : chosen_routes_id) {
            for (auto task : task_routes[route_id]) {
                candidate_tasks.push_back(task);

                if (mcgrp.is_edge(task)) {
                    candidate_tasks.push_back(mcgrp.inst_tasks[task].inverse);
                }
            }
        }


        vector<int> seq_buffer;
        seq_buffer = split_task(ns.policy, mcgrp, candidate_tasks, pseudo_capacity);
        My_Assert(seq_buffer.back() == DUMMY, "incorrect sequence!");
        seq_buffer.pop_back();

        //concatenate unchosen routes
        for (auto route_id = 0; route_id < ns.routes.size(); route_id++) {
            if (find(chosen_routes_id.begin(), chosen_routes_id.end(), route_id) == chosen_routes_id.end()) {
                seq_buffer.push_back(DUMMY);
                for (auto task : task_routes[route_id]) {
                    seq_buffer.push_back(task);
                }
            }
        }
        My_Assert(seq_buffer.back() != DUMMY, "incorrect sequence!");
        seq_buffer.push_back(DUMMY);

        //Best Accept
        Individual candidate_choice = mcgrp.parse_delimiter_seq(seq_buffer);
        if (candidate_choice.total_cost < best_choice) {
            best_buffer = seq_buffer;
        }
    }

    //Update neighbor search info
    ns.unpack_seq(best_buffer, mcgrp);
    ns.delimiter_coding_sol = get_delimiter_coding(ns.negative_coding_sol);
    //    ns.create_individual(mcgrp, ns.ns_indi);

    //check for missed
    ns.check_missed(mcgrp);
}

void merge_split(class HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int merge_size, const int pseudo_capacity)
{
    if (ns.routes.activated_route_id.size() < merge_size) {
        DEBUG_PRINT("Too few routes to merge!");
        return;
    }
    else if (ns.routes.activated_route_id.size() == merge_size) {
        DEBUG_PRINT("Merge whole routes");
        ns.clear();

        auto candidate_tasks = ns.get_tasks_set();
        auto seq_buffer = split_task(ns.policy, mcgrp, candidate_tasks, pseudo_capacity);

        ns.unpack_seq(seq_buffer, mcgrp);

        return;
    }

    My_Assert(ns.valid_sol(mcgrp), "Wrong state!");


    //generate tasks distribution in routes
    unordered_map<int, vector<int>> task_routes;
    for (auto route_id : ns.routes.activated_route_id) {
        int cur_task = ns.routes[route_id]->start;
        vector<int> buffer;
        while (cur_task != ns.routes[route_id]->end) {
            buffer.push_back(cur_task);
            cur_task = ns.solution[cur_task]->next->ID;
        }
        buffer.push_back(cur_task);
        task_routes.emplace(route_id, buffer);
    }


    vector<int> total_routes_id;
    for (auto i:ns.routes.activated_route_id) {
        total_routes_id.push_back(i);
    }
    COMB combination_manager(ns.routes.activated_route_id.size(), merge_size);
    vector<int> permutation;

    double best_choice = MAX(best_choice);
    vector<int> best_buffer;
    vector<int> best_routes;

    int count = 0;
    while (count < max_merge_set_num && combination_manager.get_combinations(permutation)) {
        count++;
        My_Assert(permutation.size() == merge_size, "Wrong state!");

        vector<int> chosen_routes_id;
        for (auto i: permutation) {
            chosen_routes_id.push_back(total_routes_id[i]);
        }

        //extract chosen tasks
        vector<int> candidate_tasks;
        for (auto route_id : chosen_routes_id) {
            My_Assert(ns.routes.activated_route_id.find(route_id) != ns.routes.activated_route_id.end(),
                      "Invalid route");


            for (auto task : task_routes.at(route_id)) {
                candidate_tasks.push_back(task);

                if (mcgrp.is_edge(task)) {
                    candidate_tasks.push_back(mcgrp.inst_tasks[task].inverse);
                }
            }
        }


        vector<int> seq_buffer = split_task(ns.policy, mcgrp, candidate_tasks, pseudo_capacity);
        My_Assert(seq_buffer.front() == DUMMY && seq_buffer.back() == DUMMY, "Incorrect sequence!");

        //Best Accept
        Individual res = mcgrp.parse_delimiter_seq(seq_buffer);
        double fitness = res.total_cost + ns.policy.beta * res.total_vio_load;
        if (fitness < best_choice) {
            best_choice = fitness;
            best_buffer = seq_buffer;
            best_routes = chosen_routes_id;
        }
    }

    //Update neighbor search info
    My_Assert(!best_buffer.empty() && !best_routes.empty(), "Wrong state");
    ns.update(mcgrp, best_buffer, best_routes);

}

vector<int> split_task(Policy &policy, const MCGRP &mcgrp, const vector<int> &tasks, const int load_constraint)
{
    vector<int> merge_sequence;

    Individual res;
    double fitness;
    double buffer;


    /*-----------------Nearest L2 distance merge policy--------------------------*/
    merge_sequence = nearest_growing(mcgrp, tasks, load_constraint);
    Individual nearest_L2_indi;
    nearest_L2_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    res = nearest_L2_indi;
    fitness = res.total_cost + policy.beta * res.total_vio_load;
    /*-----------------Nearest L2 distance merge policy--------------------------*/

    /*-----------------Furthest L2 distance merge policy--------------------------*/
    merge_sequence = nearest_depot_growing(mcgrp, tasks, load_constraint);
    Individual nearest_depot_indi;
    nearest_depot_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = nearest_depot_indi.total_cost + policy.beta * nearest_depot_indi.total_vio_load;
    if (buffer < fitness) {
        res = nearest_depot_indi;
        fitness = buffer;
    }
    /*-----------------Furthest L2 distance merge policy--------------------------*/

    /*-----------------Max yield merge policy--------------------------*/
    merge_sequence = maximum_yield_growing(mcgrp, tasks, load_constraint);
    Individual maximum_yield_indi;
    maximum_yield_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = maximum_yield_indi.total_cost + policy.beta * maximum_yield_indi.total_vio_load;
    if (buffer < fitness) {
        res = maximum_yield_indi;
        fitness = buffer;
    }
    /*-----------------Max yield merge policy--------------------------*/

    /*-----------------Min yield merge policy--------------------------*/
    merge_sequence = minimum_yield_growing(mcgrp, tasks, load_constraint);
    Individual minimum_yield_indi;
    minimum_yield_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = minimum_yield_indi.total_cost + policy.beta * minimum_yield_indi.total_vio_load;
    if (buffer < fitness) {
        res = minimum_yield_indi;
        fitness = buffer;
    }
    /*-----------------Max yield merge policy--------------------------*/

    /*-----------------mixtured merge policy--------------------------*/
    merge_sequence = mixture_growing(mcgrp, tasks, load_constraint);
    Individual mixtured_indi;
    mixtured_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = mixtured_indi.total_cost + policy.beta * mixtured_indi.total_vio_load;
    if (buffer < fitness) {
        res = mixtured_indi;
        fitness = buffer;
    }
    /*-----------------Max yield merge policy--------------------------*/

    return res.sequence;

}

vector<int> nearest_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_task;
    int current_load = 0;
    int drive_time = 0;
    vector<int> candidate_tasks;
    vector<int> nearest_tasks;

    while (!tasks.empty()) {
        current_task = sequence.back();

        candidate_tasks.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= mcgrp.inst_tasks[task].time_window.second) {
                candidate_tasks.push_back(task);
            }
        }

        //check whether need to create a new route
        if (candidate_tasks.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double min_dist = MAX(min_dist);
        nearest_tasks.clear();
        for (auto task:candidate_tasks) {
            if (mcgrp.min_cost[mcgrp.inst_tasks[current_task].tail_node][mcgrp.inst_tasks[task].head_node] < min_dist) {
                min_dist = mcgrp.min_cost[mcgrp.inst_tasks[current_task].tail_node][mcgrp.inst_tasks[task].head_node];
                nearest_tasks.clear();
                nearest_tasks.push_back(task);
            }
            else if (mcgrp.min_cost[mcgrp.inst_tasks[current_task].tail_node][mcgrp.inst_tasks[task].head_node]
                == min_dist) {
                nearest_tasks.push_back(task);
            }
            else {
                continue;
            }
        }

        My_Assert(!nearest_tasks.empty(), "you cannot have an empty nearest task set!");

        int chosen_task = nearest_tasks[0];
        drive_time = mcgrp.cal_arrive_time(sequence.back(), chosen_task, drive_time, true);
        sequence.push_back(chosen_task);
        current_load += mcgrp.inst_tasks[chosen_task].demand;

        auto ite = find(tasks.begin(), tasks.end(), chosen_task);
        My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
        tasks.erase(ite);

        if (mcgrp.is_edge(chosen_task)) {
            ite = find(tasks.begin(), tasks.end(), mcgrp.inst_tasks[chosen_task].inverse);
            My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
            tasks.erase(ite);
        }
    }

    if (sequence.back() != DUMMY) {
        sequence.push_back(DUMMY);
    }

    return sequence;
}

vector<int> nearest_depot_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> candidate_tasks;
    vector<int> nearest_tasks;

    while (!tasks.empty()) {
        candidate_tasks.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= mcgrp.inst_tasks[task].time_window.second) {
                candidate_tasks.push_back(task);
            }
        }

        //check whether need to create a new route
        if (candidate_tasks.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double min_dist = MAX(min_dist);
        nearest_tasks.clear();

        for (auto task:candidate_tasks) {
            if (mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node] < min_dist) {
                min_dist = mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
                nearest_tasks.clear();
                nearest_tasks.push_back(task);
            }
            else if (mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node] == min_dist) {
                nearest_tasks.push_back(task);
            }
            else {
                continue;
            }
        }

        My_Assert(!nearest_tasks.empty(), "you cannot have an empty nearest task set!");

        int chosen_task = nearest_tasks[0];
        drive_time = mcgrp.cal_arrive_time(sequence.back(), chosen_task, drive_time, true);
        sequence.push_back(chosen_task);
        current_load += mcgrp.inst_tasks[chosen_task].demand;


        auto ite = find(tasks.begin(), tasks.end(), chosen_task);
        My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
        tasks.erase(ite);

        if (mcgrp.is_edge(chosen_task)) {
            ite = find(tasks.begin(), tasks.end(), mcgrp.inst_tasks[chosen_task].inverse);
            My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
            tasks.erase(ite);
        }
    }

    if (sequence.back() != DUMMY) {
        sequence.push_back(DUMMY);
    }

    return sequence;
}

vector<int> maximum_yield_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> candidate_tasks;
    vector<int> maximum_tasks;

    while (!tasks.empty()) {
        candidate_tasks.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= mcgrp.inst_tasks[task].time_window.second) {
                candidate_tasks.push_back(task);
            }
        }

        //check whether need to create a new route
        if (candidate_tasks.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double max_yield = -1;
        maximum_tasks.clear();

        for (auto task:candidate_tasks) {
            if (mcgrp.get_yield(task) > max_yield) {
                max_yield = mcgrp.get_yield(task);
                maximum_tasks.clear();
                maximum_tasks.push_back(task);
            }
            else if (mcgrp.get_yield(task) == max_yield) {
                maximum_tasks.push_back(task);
            }
            else {
                continue;
            }
        }

        My_Assert(!maximum_tasks.empty(), "you cannot have an empty nearest task set!");

        int chosen_task = maximum_tasks[0];
        drive_time = mcgrp.cal_arrive_time(sequence.back(), chosen_task, drive_time, true);
        sequence.push_back(chosen_task);
        current_load += mcgrp.inst_tasks[chosen_task].demand;


        auto ite = find(tasks.begin(), tasks.end(), chosen_task);
        My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
        tasks.erase(ite);

        if (mcgrp.is_edge(chosen_task)) {
            ite = find(tasks.begin(), tasks.end(), mcgrp.inst_tasks[chosen_task].inverse);
            My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
            tasks.erase(ite);
        }
    }

    if (sequence.back() != DUMMY) {
        sequence.push_back(DUMMY);
    }

    return sequence;
}

vector<int> minimum_yield_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> candidate_tasks;
    vector<int> minimum_tasks;

    while (!tasks.empty()) {
        candidate_tasks.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= mcgrp.inst_tasks[task].time_window.second) {
                candidate_tasks.push_back(task);
            }
        }

        //check whether need to create a new route
        if (candidate_tasks.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double min_yield = MAX(min_yield);
        minimum_tasks.clear();

        for (auto task:candidate_tasks) {
            if (mcgrp.get_yield(task) < min_yield) {
                min_yield = mcgrp.get_yield(task);
                minimum_tasks.clear();
                minimum_tasks.push_back(task);
            }
            else if (mcgrp.get_yield(task) == min_yield) {
                minimum_tasks.push_back(task);
            }
            else {
                continue;
            }
        }

        My_Assert(!minimum_tasks.empty(), "you cannot have an empty nearest task set!");

        int chosen_task = minimum_tasks[0];
        drive_time = mcgrp.cal_arrive_time(sequence.back(), chosen_task, drive_time, true);
        sequence.push_back(chosen_task);
        current_load += mcgrp.inst_tasks[chosen_task].demand;


        auto ite = find(tasks.begin(), tasks.end(), chosen_task);
        My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
        tasks.erase(ite);

        if (mcgrp.is_edge(chosen_task)) {
            ite = find(tasks.begin(), tasks.end(), mcgrp.inst_tasks[chosen_task].inverse);
            My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
            tasks.erase(ite);
        }
    }

    if (sequence.back() != DUMMY) {
        sequence.push_back(DUMMY);
    }

    return sequence;
}

vector<int> mixture_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> candidate_tasks;
    vector<int> potential_tasks;

    while (!tasks.empty()) {
        candidate_tasks.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= mcgrp.inst_tasks[task].time_window.second) {
                candidate_tasks.push_back(task);
            }
        }

        //check whether need to create a new route
        if (candidate_tasks.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        if (current_load < constraint / 2) {
            //minimize depot distance
            double min_dist = MAX(min_dist);
            potential_tasks.clear();

            for (auto task:candidate_tasks) {
                if (mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node] < min_dist) {
                    min_dist = mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
                    potential_tasks.clear();
                    potential_tasks.push_back(task);
                }
                else if (mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node]
                    == min_dist) {
                    potential_tasks.push_back(task);
                }
                else {
                    continue;
                }
            }

            My_Assert(!potential_tasks.empty(), "you cannot have an empty nearest task set!");

            int chosen_task = potential_tasks[0];
            drive_time = mcgrp.cal_arrive_time(sequence.back(), chosen_task, drive_time, true);
            sequence.push_back(chosen_task);
            current_load += mcgrp.inst_tasks[chosen_task].demand;


            auto ite = find(tasks.begin(), tasks.end(), chosen_task);
            My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
            tasks.erase(ite);

            if (mcgrp.is_edge(chosen_task)) {
                ite = find(tasks.begin(), tasks.end(), mcgrp.inst_tasks[chosen_task].inverse);
                My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
                tasks.erase(ite);
            }
        }
        else {
            //maximum yield merging
            double max_yield = -1;
            potential_tasks.clear();

            for (auto task:candidate_tasks) {
                if (mcgrp.get_yield(task) > max_yield) {
                    max_yield = mcgrp.get_yield(task);
                    potential_tasks.clear();
                    potential_tasks.push_back(task);
                }
                else if (mcgrp.get_yield(task) == max_yield) {
                    potential_tasks.push_back(task);
                }
                else {
                    continue;
                }
            }

            My_Assert(!potential_tasks.empty(), "you cannot have an empty nearest task set!");

            int chosen_task = potential_tasks[0];
            sequence.push_back(chosen_task);
            current_load += mcgrp.inst_tasks[chosen_task].demand;

            auto ite = find(tasks.begin(), tasks.end(), chosen_task);
            My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
            tasks.erase(ite);

            if (mcgrp.is_edge(chosen_task)) {
                ite = find(tasks.begin(), tasks.end(), mcgrp.inst_tasks[chosen_task].inverse);
                My_Assert(ite != tasks.end(), "Cannot find chosen tasks!");
                tasks.erase(ite);
            }
        }
    }

    if (sequence.back() != DUMMY) {
        sequence.push_back(DUMMY);
    }

    return sequence;
}