//
// Created by luke on 2019/11/29.
//

#include "ConstructPolicy.h"
#include "RNG.h"
#include <bits/stdc++.h>
#include "local_search.h"
#include "SearchPolicy.h"



using namespace std;

const int max_merge_set_num = 20;


NearestScanner::NearestScanner(const MCGRPTW &mcgrp, Distance& distance_)
    : PathConstructor(mcgrp,"nearest scanning"), distance(distance_)
{}

Individual NearestScanner::operator()(const vector<int> &taskList, const string& mode)
{
    int load;
    double min_dist;
    int drive_time; // the earliest time of a vehicle begins to serve the Task

    std::vector<int> FCL; // feasible candidate list
    std::vector<int> nearest_task_set;

    int current_tail_task;
    int chosen_task;

    unordered_set<int> unserved_task_id_set;

    if(taskList.empty()){
        for (int i = 1; i <= mcgrp.actual_task_num; i++)
            unserved_task_id_set.insert(i);
    }else{
        unserved_task_id_set = unordered_set<int>(taskList.begin(), taskList.end());
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    load = 0;
    drive_time = 0;

    vector<int> solution;
    solution.push_back(DUMMY);

    while (!unserved_task_id_set.empty()) {
        current_tail_task = solution.back();

        //Find all the tasks that satisfy the capacity constraint
        if(mode == "feasible"){
            FCL.clear();
            for (auto unserved_task : unserved_task_id_set) {
                if (mcgrp.inst_tasks[unserved_task].demand <= mcgrp.capacity - load
                    && mcgrp.cal_arrive_time(current_tail_task, unserved_task, drive_time, true)
                        <= mcgrp.inst_tasks[unserved_task].time_window.second) {
                    FCL.push_back(unserved_task);
                }
            }
        }else if(mode == "allow_infeasible"){
            FCL = vector<int>(unserved_task_id_set.begin(), unserved_task_id_set.end());
        }else{
            My_Assert(false,"error! unknown mode!");
        }


        if (FCL.empty()) {
            solution.push_back(DUMMY);
            load = 0;
            drive_time = 0;
            continue;
        }

        min_dist = MAX(min_dist);

        //Find the nearest Task from the current candidate Task set
        for (auto candidate_task : FCL) {
            if (distance(mcgrp,current_tail_task,candidate_task) < min_dist) {
                min_dist = distance(mcgrp,current_tail_task,candidate_task);
                nearest_task_set.clear();
                nearest_task_set.push_back(candidate_task);
            }
            else if (
                distance(mcgrp,current_tail_task,candidate_task) == min_dist) {
                nearest_task_set.push_back(candidate_task);
            }
        }

        //If multiple tasks both satisfy the capacity constraint and are closest, randomly choose one
        int k = (int) mcgrp._rng.Randint(0, nearest_task_set.size() - 1);
        chosen_task = nearest_task_set[k];

        load += mcgrp.inst_tasks[chosen_task].demand;
        drive_time = mcgrp.cal_arrive_time(current_tail_task, chosen_task, drive_time, true);

        solution.push_back(chosen_task);

        unserved_task_id_set.erase(chosen_task);

        if (mcgrp.is_edge(chosen_task)) {
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

void
merge_split(LocalSearch &ns,
            const MCGRPTW &mcgrp, const int merge_size)
{
    if (ns.routes.activated_route_id.size() < merge_size) {
        DEBUG_PRINT("Too few routes to merge!");
        return;
    }
    else if (ns.routes.activated_route_id.size() == merge_size) {
        DEBUG_PRINT("Merge whole routes");
        ns.clear();

        auto candidate_tasks = ns.get_tasks_set();
        auto seq_buffer = split_task(mcgrp,ns, candidate_tasks, ns.policy);

        ns.unpack_seq(seq_buffer, mcgrp);

        return;
    }

    My_Assert(ns.valid_sol(mcgrp), "Wrong state!");


    //generate tasks distribution in routes

    unordered_map<int, vector<int>> task_routes;
    for (auto route_id : ns.routes.activated_route_id) {
        int cur_task = ns.solution[ns.routes[route_id]->start]->next->ID;
        vector<int> buffer;
        while (cur_task != ns.solution[ns.routes[route_id]->end]->pre->ID) {
            buffer.push_back(cur_task);
            cur_task = ns.solution[cur_task]->next->ID;
        }
        buffer.push_back(cur_task);

        task_routes.emplace(route_id, buffer);
    }

    auto stable_scores = ns.cal_route_scores(mcgrp);
    vector<int> routes_id;
    for (int i = 0; i < merge_size;i++) {
        routes_id.push_back(stable_scores[i].route_id);
    }

    double best_choice = MAX(best_choice);
    vector<int> best_buffer;
    vector<int> best_routes;

    My_Assert(routes_id.size() == merge_size, "Wrong state!");


    //extract chosen tasks
    vector<int> candidate_tasks;
    for (auto route_id : routes_id) {
        My_Assert(ns.routes.activated_route_id.find(route_id) != ns.routes.activated_route_id.end(),
                  "Invalid route");


        for (auto task : task_routes.at(route_id)) {
            candidate_tasks.push_back(task);

            if (mcgrp.is_edge(task)) {
                candidate_tasks.push_back(mcgrp.inst_tasks[task].inverse);
            }
        }
    }

    vector<int> seq_buffer = split_task(mcgrp,ns, candidate_tasks, ns.policy);
    My_Assert(seq_buffer.front() == DUMMY && seq_buffer.back() == DUMMY, "Incorrect sequence!");

    //Best Accept
    Individual res = mcgrp.parse_delimiter_seq(seq_buffer);
    double fitness = ns.getFitness(mcgrp,ns.policy,res);
    if (fitness < best_choice) {
        best_choice = fitness;
        best_buffer = seq_buffer;
        best_routes = routes_id;
    }

    //Update neighbor search info
    My_Assert(!best_buffer.empty() && !best_routes.empty(), "Wrong state");
    ns.update(mcgrp, best_buffer, best_routes);

}

vector<int> split_task(const MCGRPTW &mcgrp, LocalSearch& ns, const vector<int> &tasks, Policy &policy)
{
    vector<int> merge_sequence;

    Individual res;
    double fitness;
    double buffer;


    /*-----------------Nearest L2 distance merge policy--------------------------*/
    merge_sequence = nearest_growing(mcgrp, tasks, policy);
    Individual nearest_L2_indi;
    nearest_L2_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    res = nearest_L2_indi;
    fitness = ns.getFitness(mcgrp,ns.policy,res);
    /*-----------------Nearest L2 distance merge policy--------------------------*/

    /*-----------------Furthest L2 distance merge policy--------------------------*/
    merge_sequence = nearest_depot_growing(mcgrp, tasks, policy);
    Individual nearest_depot_indi;
    nearest_depot_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = ns.getFitness(mcgrp,ns.policy,nearest_depot_indi);
    if (buffer < fitness) {
        res = nearest_depot_indi;
        fitness = buffer;
    }
    /*-----------------Furthest L2 distance merge policy--------------------------*/

    /*-----------------Max yield merge policy--------------------------*/
    merge_sequence = maximum_yield_growing(mcgrp, tasks, policy);
    Individual maximum_yield_indi;
    maximum_yield_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = ns.getFitness(mcgrp,ns.policy,maximum_yield_indi);
    if (buffer < fitness) {
        res = maximum_yield_indi;
        fitness = buffer;
    }
    /*-----------------Max yield merge policy--------------------------*/

    /*-----------------Min yield merge policy--------------------------*/
    merge_sequence = minimum_yield_growing(mcgrp, tasks, policy);
    Individual minimum_yield_indi;
    minimum_yield_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = ns.getFitness(mcgrp,ns.policy,minimum_yield_indi);
    if (buffer < fitness) {
        res = minimum_yield_indi;
        fitness = buffer;
    }
    /*-----------------Max yield merge policy--------------------------*/

    /*-----------------mixture merge policy--------------------------*/
    merge_sequence = mixture_growing(mcgrp, tasks,policy);
    Individual mixture_indi;
    mixture_indi = mcgrp.parse_delimiter_seq(merge_sequence);
    buffer = ns.getFitness(mcgrp,ns.policy,mixture_indi);
    if (buffer < fitness) {
        res = mixture_indi;
        fitness = buffer;
    }
    /*-----------------Max yield merge policy--------------------------*/

    return res.sequence;

}

vector<int> nearest_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_task;
    int current_load = 0;
    int drive_time = 0;
    vector<int> FCL;
    vector<int> nearest_tasks;

    int constraint = policy.get_pseudo_capacity(mcgrp.capacity);
    while (!tasks.empty()) {
        current_task = sequence.back();

        FCL.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= policy.get_pseudo_time_window(mcgrp.inst_tasks[task].time_window.second)) {
                FCL.push_back(task);
            }
        }

        //check whether need to create a new route
        if (FCL.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double min_dist = MAX(min_dist);
        nearest_tasks.clear();
        for (auto task:FCL) {
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

        My_Assert(!nearest_tasks.empty(), "you cannot have an empty nearest Task set!");

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

vector<int> nearest_depot_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> FCL;
    vector<int> nearest_tasks;

    int constraint = policy.get_pseudo_capacity(mcgrp.capacity);
    while (!tasks.empty()) {
        FCL.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= policy.get_pseudo_time_window(mcgrp.inst_tasks[task].time_window.second)) {
                FCL.push_back(task);
            }
        }

        //check whether need to create a new route
        if (FCL.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double min_dist = MAX(min_dist);
        nearest_tasks.clear();

        for (auto task:FCL) {
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

        My_Assert(!nearest_tasks.empty(), "you cannot have an empty nearest Task set!");

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

vector<int> maximum_yield_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> FCL;
    vector<int> maximum_tasks;

    int constraint = policy.get_pseudo_capacity(mcgrp.capacity);
    while (!tasks.empty()) {
        FCL.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= policy.get_pseudo_time_window(mcgrp.inst_tasks[task].time_window.second)) {
                FCL.push_back(task);
            }
        }

        //check whether need to create a new route
        if (FCL.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double max_yield = -1;
        maximum_tasks.clear();

        for (auto task:FCL) {
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

        My_Assert(!maximum_tasks.empty(), "you cannot have an empty nearest Task set!");

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

vector<int> minimum_yield_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> FCL;
    vector<int> minimum_tasks;

    int constraint = policy.get_pseudo_capacity(mcgrp.capacity);
    while (!tasks.empty()) {
        FCL.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= policy.get_pseudo_time_window(mcgrp.inst_tasks[task].time_window.second)) {
                FCL.push_back(task);
            }
        }

        //check whether need to create a new route
        if (FCL.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        double min_yield = MAX(min_yield);
        minimum_tasks.clear();

        for (auto task:FCL) {
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

        My_Assert(!minimum_tasks.empty(), "you cannot have an empty nearest Task set!");

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

vector<int> mixture_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy)
{
    vector<int> sequence;

    sequence.push_back(DUMMY);
    int current_load = 0;
    int drive_time = 0;
    vector<int> FCL;
    vector<int> potential_tasks;

    int constraint = policy.get_pseudo_capacity(mcgrp.capacity);
    while (!tasks.empty()) {
        FCL.clear();
        for (auto task : tasks) {
            if (mcgrp.inst_tasks[task].demand <= constraint - current_load
                && mcgrp.cal_arrive_time(sequence.back(), task, drive_time, true)
                    <= policy.get_pseudo_time_window(mcgrp.inst_tasks[task].time_window.second)) {
                FCL.push_back(task);
            }
        }

        //check whether need to create a new route
        if (FCL.empty()) {
            sequence.push_back(DUMMY);
            current_load = 0;
            drive_time = 0;
            continue;
        }

        if (current_load < constraint / 2) {
            //minimize depot distance
            double min_dist = MAX(min_dist);
            potential_tasks.clear();

            for (auto task:FCL) {
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

            My_Assert(!potential_tasks.empty(), "you cannot have an empty nearest Task set!");

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

            for (auto task:FCL) {
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

            My_Assert(!potential_tasks.empty(), "you cannot have an empty nearest Task set!");

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

vector<vector<int>> tour_splitting(const MCGRPTW &mcgrp, vector<int>& task_list)
{
    task_list.insert(task_list.begin(),DUMMY);
    int nsize = task_list.size();

    vector<int> W(nsize,0);
    vector<int> P(nsize,DUMMY);

    for(int i = 1;i<W.size();i++)
        W[i] = INT32_MAX;

    for(int i = 1;i<nsize;i++){
        int j = i;
        int load = 0;
        int length = 0;
        int DepartureTime = 0;
        int u = 0;
        while (true){
            int v = j;
            load += mcgrp.inst_tasks[task_list[v]].demand;
            length = length - mcgrp.min_cost[mcgrp.inst_tasks[task_list[u]].tail_node][mcgrp.inst_tasks[DUMMY].head_node]
                + mcgrp.min_cost[mcgrp.inst_tasks[task_list[u]].tail_node][mcgrp.inst_tasks[task_list[v]].head_node]
                + mcgrp.inst_tasks[task_list[v]].serv_cost
                + mcgrp.min_cost[mcgrp.inst_tasks[task_list[v]].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

            int ArrivalTime = DepartureTime
                + mcgrp.min_time[mcgrp.inst_tasks[task_list[u]].tail_node][mcgrp.inst_tasks[task_list[v]].head_node];

            DepartureTime = ArrivalTime + max(0, mcgrp.inst_tasks[task_list[v]].time_window.first - ArrivalTime)
                + mcgrp.inst_tasks[task_list[v]].serve_time;

            if (load <= mcgrp.capacity
                && W[i - 1] + length < W[j]
                && ArrivalTime <= mcgrp.inst_tasks[task_list[v]].time_window.second){
                W[j] = W[i - 1] + length;
                P[j] = i - 1;
            }

            j += 1;
            u = v;
            if(j >= nsize || load > mcgrp.capacity || ArrivalTime > mcgrp.inst_tasks[task_list[u]].time_window.second){
                break;
            }
        }
    }

    vector<vector<int>> routes;
    int cur = nsize - 1;
    int pred = P[cur];
    while(cur != 0){
        vector<int> route;
        for(int i = pred + 1;i <= cur;i++)
            route.push_back(task_list[i]);
        routes.emplace_back(route);
        cur = pred;
        pred = P[cur];
    }

    task_list.erase(task_list.begin());
    return routes;
}

BestServeSeq viterbi_decoding(const MCGRPTW &mcgrp, const vector<int>& task_list, bool allow_infeasible)
{
    BestServeSeq ans;

    if(task_list.empty()){
        ans.cost = 0;
        ans.sequence.clear();
        ans.arrive_time_tbl.clear();
        return ans;
    }

    vector<vector<int>> DAG;
    DAG.push_back(vector<int>(1,DUMMY));
    for(int i = 0;i<task_list.size();i++){
        if(mcgrp.is_edge(task_list[i])){
            DAG.push_back(vector<int>{task_list[i], mcgrp.inst_tasks[task_list[i]].inverse});
        }else{
            DAG.push_back(vector<int>{task_list[i]});
        }
    }
    DAG.push_back(vector<int>(1,DUMMY));


    vector<vector<int>> W(DAG.size(),vector<int>());
    vector<vector<int>> P(DAG.size(),vector<int>());
    vector<vector<int>> ArriveTime(DAG.size(), vector<int>());

    W[0].push_back(0);
    P[0].push_back(0);
    ArriveTime[0].push_back(0);

    function<pair<int,int>(const vector<int>&)> argmin_ = [](const vector<int>& seq){
        int index = -1;
        int val = INT32_MAX;
        for(int ii = 0;ii < seq.size();ii++){
            if (seq[ii] < val){
                index  = ii;
                val = seq[ii];
            }
        }

        return make_pair(index, val);
    };


    for(int i = 1; i < DAG.size();i++){
        for(int j = 0;j < DAG[i].size() ;j++){
            vector<int> distance_;
            vector<int> ArriveTime_;
            int LatestDepartureTime = allow_infeasible ? INT32_MAX : mcgrp.inst_tasks[DAG[i][j]].time_window.second;

            for(int k = 0; k < W[i-1].size();k++){
                if(W[i - 1][k] == INT32_MAX ||
                    ArriveTime[i - 1][k] == INT32_MAX ||
                    ArriveTime[i - 1][k] + mcgrp.inst_tasks[DAG[i - 1][k]].serve_time +
                        mcgrp.min_time[mcgrp.inst_tasks[DAG[i - 1][k]].tail_node][mcgrp.inst_tasks[DAG[i][j]].head_node]
                        > LatestDepartureTime){
                    distance_.push_back(INT32_MAX);
                    ArriveTime_.push_back(INT32_MAX);
                }
                else{
                    distance_.push_back(W[i - 1][k] +
                        mcgrp.inst_tasks[DAG[i - 1][k]].serv_cost +
                        mcgrp.min_cost[mcgrp.inst_tasks[DAG[i - 1][k]].tail_node][mcgrp.inst_tasks[DAG[i][j]].head_node]);
                    ArriveTime_.push_back(max(mcgrp.inst_tasks[DAG[i][j]].time_window.first,
                                              ArriveTime[i - 1][k] + mcgrp.inst_tasks[DAG[i - 1][k]].serve_time +
                                                  mcgrp.min_time[mcgrp.inst_tasks[DAG[i-1][k]].tail_node][mcgrp.inst_tasks[DAG[i][j]].head_node]));
                }
            }

            int idx;
            int val;
            auto tmp = argmin_(distance_);
            idx = tmp.first;
            val = tmp.second;

            P[i].push_back(idx);
            W[i].push_back(val);
            if (idx == -1)
                ArriveTime[i].push_back(INT32_MAX);
            else
                ArriveTime[i].push_back(ArriveTime_[idx]);
        }
    }

    if (P.back() == vector<int>{-1}){
        ans.cost = INT_MAX;
    }else{
        My_Assert(W.back().size() == 1, "Wrong status!");
        ans.cost = W.back().back();
        vector<int> indices{P.back().back()};
        for(int ii = (int)P.size() - 2; ii >= 2 ; ii--){
            indices.push_back(P[ii][indices.back()]);
        }
        My_Assert(indices.size() == P.size() - 2, "Wrong decode sequence!");

        reverse(indices.begin(),indices.end());
        for(int ii = 0;ii < (int)indices.size(); ii++){
            ans.sequence.push_back(DAG[ii+1][indices[ii]]);
            ans.arrive_time_tbl.push_back(ArriveTime[ii+1][indices[ii]]);
        }
    }

    return ans;
}

PathConstructor::PathConstructor(const MCGRPTW &mcgrp_, const string& name_) : mcgrp(mcgrp_)
{
    name = name_;
}


Individual RTFScanner::operator()(const vector<int> &taskList, const string &mode)
{
    int load;
    double min_dist;
    int arrival_time;

    unordered_set<int> unserved_task_id;

    std::vector<int> FCL; // feasible candidate list

    std::vector<int> RCL1; // restricted candidate list which drive away from the depot
    std::vector<int> RCL2; // restricted candidate list which drive towards the depot

    if(!taskList.empty())
        unserved_task_id = unordered_set<int>(taskList.begin(), taskList.end());
    else{
        for (int i = 1; i <= mcgrp.actual_task_num; i++)
            unserved_task_id.insert(i);
    }

    load = 0;
    arrival_time = 0;

    vector<int> solution;
    solution.push_back(DUMMY);

    while (!unserved_task_id.empty()){
        int tail_task = solution.back();

        FCL.clear();
        if(mode == "allow_infeasible"){
            FCL = vector<int>(unserved_task_id.begin(), unserved_task_id.end());
        }
        else if(mode == "feasible"){
            for (auto unserved_task : unserved_task_id) {
                if (mcgrp.inst_tasks[unserved_task].demand <= mcgrp.capacity - load
                    && mcgrp.cal_arrive_time(tail_task, unserved_task, arrival_time, true)
                        <= mcgrp.inst_tasks[unserved_task].time_window.second) {
                    FCL.push_back(unserved_task);
                }
            }
        }else{
            My_Assert(false,"error! unknown mode!");
        }


        if (FCL.empty()) {
            solution.push_back(DUMMY);
            load = 0;
            arrival_time = 0;
            continue;
        }

        min_dist = MAX(min_dist);
        for (const auto& task : FCL) {
            double d_uz = distance(mcgrp,tail_task,task);
            double d_udummy = distance(mcgrp,tail_task,DUMMY);
            double d_zdummy = distance(mcgrp,task,DUMMY);

            if (d_uz > min_dist) continue;

            if (d_uz < min_dist) {
                min_dist = d_uz;
                RCL1.clear();
                RCL2.clear();
            }

            if(d_udummy < d_zdummy) RCL1.push_back(task);
            else RCL2.push_back(task);
        }

        int chosen_task = -1;
        if(RCL1.empty()){
            int seed = (int) mcgrp._rng.Randint(0, RCL2.size() - 1);
            chosen_task = RCL2[seed];
        }
        else if(RCL2.empty()){
            int seed = (int) mcgrp._rng.Randint(0, RCL1.size() - 1);
            chosen_task = RCL1.front();
        }else{
            if((load % mcgrp.capacity) < (mcgrp.capacity / 2)){
                int seed = (int) mcgrp._rng.Randint(0, RCL1.size() - 1);
                chosen_task = RCL1.front();
            }else{
                int seed = (int) mcgrp._rng.Randint(0, RCL2.size() - 1);
                chosen_task = RCL2[seed];
            }
        }

        My_Assert(chosen_task != -1,"error! cannot find a Task");

        load += mcgrp.inst_tasks[chosen_task].demand;
        arrival_time = mcgrp.cal_arrive_time(tail_task, chosen_task, arrival_time, true);

        solution.push_back(chosen_task);

        unserved_task_id.erase(chosen_task);

        if (mcgrp.is_edge(chosen_task)) {
            int inverse_task = mcgrp.inst_tasks[chosen_task].inverse;
            if(unserved_task_id.find(inverse_task) != unserved_task_id.end())
                unserved_task_id.erase(inverse_task);
        }

    }


    if (solution.back() != DUMMY) solution.push_back(DUMMY);

    return mcgrp.parse_delimiter_seq(solution);
}

RTFScanner::RTFScanner(const MCGRPTW &mcgrp, Distance &distance_)
    : PathConstructor(mcgrp, "RTF scanner"), distance(distance_)
{}

SampleScanner::SampleScanner(const MCGRPTW &mcgrp, LearningDistance &distance_, int sample_times_)
    : PathConstructor(mcgrp, "sample scanner"), distance(distance_), sample_times(sample_times_)
{}

Individual SampleScanner::operator()(const vector<int> &taskList, const string &mode)
{
    int load;
    double min_dist;
    int drive_time; // the earliest time of a vehicle begins to serve the Task

    std::unordered_set<int> FCL; // feasible candidate list
    std::vector<int> nearest_task_set;

    int current_tail_task;
    int chosen_task;

    unordered_set<int> unserved_task_id_set;

    if(taskList.empty()){
        for (int i = 1; i <= mcgrp.actual_task_num; i++)
            unserved_task_id_set.insert(i);
    }else{
        unserved_task_id_set = unordered_set<int>(taskList.begin(), taskList.end());
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    load = 0;
    drive_time = 0;

    vector<int> solution;
    solution.push_back(DUMMY);

    while (!unserved_task_id_set.empty()) {
        current_tail_task = solution.back();

        const vector<double>& prob_vector = distance.get_prob_vector(current_tail_task);

        FCL.clear();
        for(int sample_cnt = 0; sample_cnt < sample_times; sample_cnt++){
            //Find all the tasks that satisfy the capacity constraint
            int sample_task = sample::uniform_sample(prob_vector);

            if(unserved_task_id_set.find(sample_task) == unserved_task_id_set.end())
                continue;

            if(mode == "feasible"){
                if (mcgrp.inst_tasks[sample_task].demand <= mcgrp.capacity - load
                    && mcgrp.cal_arrive_time(current_tail_task, sample_task, drive_time, true)
                        <= mcgrp.inst_tasks[sample_task].time_window.second) {
                    FCL.insert(sample_task);
                }
            }else if(mode == "allow_infeasible"){
                FCL.insert(sample_task);
            }else{
                My_Assert(false,"error! unknown mode!");
            }
        }


        if (FCL.empty()) {
            if(solution.back() == DUMMY){
                FCL = unserved_task_id_set;
            }else{
                solution.push_back(DUMMY);
                load = 0;
                drive_time = 0;
                continue;
            }
        }

        min_dist = MAX(min_dist);

        //Find the nearest Task from the current candidate Task set
        for (auto candidate_task : FCL) {
            if (distance(mcgrp,current_tail_task,candidate_task) < min_dist) {
                min_dist = distance(mcgrp,current_tail_task,candidate_task);
                nearest_task_set.clear();
                nearest_task_set.push_back(candidate_task);
            }
            else if (
                distance(mcgrp,current_tail_task,candidate_task) == min_dist) {
                nearest_task_set.push_back(candidate_task);
            }
        }

        //If multiple tasks both satisfy the capacity constraint and are closest, randomly choose one
        int k = (int) mcgrp._rng.Randint(0, nearest_task_set.size() - 1);
        chosen_task = nearest_task_set[k];

        load += mcgrp.inst_tasks[chosen_task].demand;
        drive_time = mcgrp.cal_arrive_time(current_tail_task, chosen_task, drive_time, true);

        solution.push_back(chosen_task);

        unserved_task_id_set.erase(chosen_task);

        if (mcgrp.is_edge(chosen_task)) {
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
