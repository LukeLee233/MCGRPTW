//
// Created by luke on 2020/11/5.
//

#include "insert.h"

vector<vector<RouteInfo::TimeTable>> XPostInsert::expected_time_table(LocalSearch &ns,
                                                                      const MCGRPTW &mcgrp,
                                                                      vector<int> &disturbance_seq,
                                                                      const int i,
                                                                      bool allow_infeasible)
{
    vector<vector<RouteInfo::TimeTable>>
        res(2,vector<RouteInfo::TimeTable>({{-1, -1}}));

    const int i_route = ns.solution[i]->route_id;
    const int u_route = ns.solution[disturbance_seq.front()]->route_id;

    auto intermediate = mcgrp.forecast_time_table(ns.routes[u_route]->time_table,
                                                  disturbance_seq,"remove",-1 ,allow_infeasible);

    if(!mcgrp.isTimeTableFeasible(intermediate)){
        return res;
    }

    vector<RouteInfo::TimeTable> final;
    if(u_route == i_route){
        final = mcgrp.forecast_time_table(intermediate,
                                          disturbance_seq,"insert_after",i,allow_infeasible);
        intermediate = final;
    }else{
        final = mcgrp.forecast_time_table(ns.routes[i_route]->time_table,
                                          disturbance_seq,"insert_after",i,allow_infeasible);
    }

    if(!mcgrp.isTimeTableFeasible(final)){
        return res;
    }

    // 0 for u_route, 1 for i_route
    res[0] = intermediate;
    res[1] = final;

    return res;
}

bool XPostInsert::considerable_move(LocalSearch &ns,
                                    const MCGRPTW &mcgrp,
                                    vector<int> move_seq,
                                    const int u)
{
    // This move used to move the `move_seq` to the back of the task `u`
    // Two policy are contained, 1: no mode choice 2: best mode choice

    move_result.reset();
    call_times++;

    // Task u cannot be dummy Task
    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong Task");
    My_Assert(all_of(move_seq.begin(), move_seq.end(),
                     [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong Task");

    if(u == ns.solution[move_seq.front()]->pre->ID){
        // Nothing to do
        return false;
    }

    const int i_route = ns.solution[move_seq.front()]->route_id;
    const int u_route = ns.solution[u]->route_id;

    // load legality is independent with mode choice
    int load_delta = 0;
    double vio_load_delta = 0;
    if (i_route != u_route) {
        for (auto task : move_seq) {
            load_delta += mcgrp.inst_tasks[task].demand;
        }

        if (ns.policy.has_rule(FEASIBLE)) {
            if (ns.routes[u_route]->load + load_delta > mcgrp.capacity) {
                return false;
            }
        }
        else if (ns.policy.has_rule(INFEASIBLE)) {
            int pseudo_capacity = ns.policy.get_pseudo_capacity(mcgrp.capacity);
            if (ns.routes[u_route]->load + load_delta > pseudo_capacity) {
                return false;
            }

            //u_route vio-load calculate
            if (ns.routes[u_route]->load + load_delta > mcgrp.capacity) {
                //if insert Task to route u and over load
                if (ns.routes[u_route]->load >= mcgrp.capacity) {
                    //if the route u already over loaded
                    vio_load_delta += load_delta;
                }
                else {
                    vio_load_delta += ns.routes[u_route]->load + load_delta - mcgrp.capacity;
                }
            }

            //i_route vio-load calculate
            if (ns.routes[i_route]->load > mcgrp.capacity) {
                //if remove Task from route i and over load
                if (ns.routes[i_route]->load - load_delta >= mcgrp.capacity) {
                    //if still over loaded
                    vio_load_delta -= load_delta;
                }
                else {
                    vio_load_delta -= (ns.routes[i_route]->load - mcgrp.capacity);
                }
            }
        }
    }


    if(ns.policy.has_rule(BEST_ACCEPT) && mcgrp.count_edges(move_seq) > 0){
//    if(false){
        // with mode refine

        vector<int> first_seq;
        vector<int> last_seq;

        if(i_route == u_route){
            // special case
            vector<int> u_route_seq = ns.get_sub_seq(ns.routes[ns.solution[u]->route_id]->start, ns.routes[ns.solution[u]->route_id]->end);

            // try to use offset to generate sequence
            int move_seq_loc = -1;
            int u_loc = -1;
            for(int ii = 0;ii < u_route_seq.size(); ii++){
                if(u_route_seq[ii] == u) u_loc = ii;
                if(u_route_seq[ii] == move_seq.front()) move_seq_loc = ii;
            }

            if(move_seq_loc < u_loc){
                first_seq.insert (first_seq.end(),u_route_seq.begin(),u_route_seq.begin() + move_seq_loc);
                first_seq.insert (first_seq.end(), u_route_seq.begin() + move_seq_loc + move_seq.size(), u_route_seq.begin() + u_loc + 1);

                last_seq.insert (last_seq.end(), u_route_seq.begin() + u_loc + 1, u_route_seq.end());
            }
            else if(move_seq_loc > u_loc){
                first_seq.insert (first_seq.end(),u_route_seq.begin(),u_route_seq.begin() + u_loc + 1);

                last_seq.insert (last_seq.end(), u_route_seq.begin() + u_loc + 1,u_route_seq.begin() + move_seq_loc);
                last_seq.insert( last_seq.end(), u_route_seq.begin() + move_seq_loc + move_seq.size(), u_route_seq.end());
            }else{
                My_Assert(false,"error, wrong location");
            }
        }
        else{
            // trivial case
            first_seq = ns.get_sub_seq(ns.routes[ns.solution[u]->route_id]->start, u);
            last_seq = ns.get_sub_seq(ns.solution[u]->next->ID, ns.routes[ns.solution[u]->route_id]->end);
        }

        viterbi::PseudoTask first_task = Seq2Pseudo(mcgrp,first_seq, 1);
        viterbi::PseudoTask last_task = Seq2Pseudo(mcgrp,last_seq, 2);

        // insert sequence will invoke viterbi result
        auto best_sequence = viterbi::viterbi_decode(&first_task, &last_task,move_seq,mcgrp,ns.policy);


        // use best sequence to update the move result
        if(best_sequence.cost == INT32_MAX){
            return false;
        }

        My_Assert(best_sequence.seq.front() == INT32_MAX && best_sequence.seq.back() == INT32_MAX, "error, wrong sequence");

        move_result.choose_tasks(move_seq.front(), u);

        move_result.move_arguments["input_seq"] = move_seq;
        move_result.move_arguments["output_seq"] = vector<int>(best_sequence.seq.begin() + 1, best_sequence.seq.end() - 1);
        move_result.move_arguments["insert_pos"] = {u};

        double u_route_length_delta = best_sequence.cost - ns.routes[u_route]->length;
        double i_route_length_delta = 0;
        if(i_route != u_route){
            auto routes_cost_delta = cost_delta(ns,mcgrp,move_seq,u);
            i_route_length_delta = routes_cost_delta.second;
        }

        if(i_route == u_route){
            move_result.num_affected_routes = 1;

            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;
            move_result.route_lens.push_back(ns.routes[i_route]->length + move_result.delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

            move_result.considerable = true;

            vector<int> u_route_seq;
            u_route_seq.insert(u_route_seq.end(), first_seq.begin() + 1, first_seq.end());
            u_route_seq.insert(u_route_seq.end(), move_result.move_arguments["output_seq"].begin(), move_result.move_arguments["output_seq"].end());
            u_route_seq.insert(u_route_seq.end(), last_seq.begin(), last_seq.end() - 1);

            move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(u_route_seq, mcgrp.cal_arrive_time(u_route_seq)));


            if(ns.policy.has_rule(INFEASIBLE)){
                auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                move_result.vio_load_delta = vio_load_delta;
                move_result.vio_time_delta = new_time_info_i.second - old_time_info_i.second;
                move_result.vio_time_custom_num_delta = new_time_info_i.first - old_time_info_i.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }
        else{

            // Different routes!
            move_result.num_affected_routes = 2;
            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;

            move_result.route_lens.push_back(ns.routes[i_route]->length + i_route_length_delta);
            move_result.route_lens.push_back(ns.routes[u_route]->length + u_route_length_delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load - load_delta);
            move_result.route_loads.push_back(ns.routes[u_route]->load + load_delta);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers - move_seq.size());
            move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers + move_seq.size());

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

            move_result.considerable = true;

            vector<int> i_route_seq;
            for(int ii = 0; ii < ns.routes[i_route]->time_table.size();){
                if(move_seq.front() == ns.routes[i_route]->time_table[ii].task){
                    ii += move_seq.size();
                    continue;
                }

                i_route_seq.push_back(ns.routes[i_route]->time_table[ii].task);
                ii ++;
            }
            move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(i_route_seq, mcgrp.cal_arrive_time(i_route_seq)));

            vector<int> u_route_seq;
            u_route_seq.insert(u_route_seq.end(), first_seq.begin() + 1, first_seq.end());
            u_route_seq.insert(u_route_seq.end(), move_result.move_arguments["output_seq"].begin(), move_result.move_arguments["output_seq"].end());
            u_route_seq.insert(u_route_seq.end(), last_seq.begin(), last_seq.end() - 1);

            move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(u_route_seq, mcgrp.cal_arrive_time(u_route_seq)));

            if(ns.policy.has_rule(INFEASIBLE)){
                move_result.vio_load_delta = vio_load_delta;

                auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
                auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

                move_result.vio_time_delta =
                    new_time_info_i.second + new_time_info_u.second - old_time_info_i.second - old_time_info_u.second;
                move_result.vio_time_custom_num_delta =
                    new_time_info_i.first + new_time_info_u.first - old_time_info_i.first - old_time_info_u.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }

    }
    else {
        // without mode refine

        bool allow_infeasible = ns.policy.has_rule(INFEASIBLE) ? true : false;
        vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};

        new_time_tbl = expected_time_table(ns, mcgrp, move_seq, u, allow_infeasible);

        if(!mcgrp.isTimeTableFeasible(new_time_tbl[0])){
            return false;
        }

        if(allow_infeasible){
            if(ns.policy.check_time_window(mcgrp,new_time_tbl)){
                return false;
            }
        }

        auto routes_cost_delta = cost_delta(ns,mcgrp,move_seq,u);
        double u_route_length_delta = routes_cost_delta.first;
        double i_route_length_delta = routes_cost_delta.second;

        move_result.choose_tasks(move_seq.front(), u);

        move_result.move_arguments["input_seq"] = move_seq;
        move_result.move_arguments["output_seq"] = move_seq;
        move_result.move_arguments["insert_pos"] = {u};

        if (i_route == u_route) {

            move_result.num_affected_routes = 1;

            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;
            move_result.route_lens.push_back(ns.routes[i_route]->length + move_result.delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

            move_result.considerable = true;

            move_result.route_time_tbl.emplace_back(new_time_tbl[0]);

            if(ns.policy.has_rule(INFEASIBLE)){
                auto old_time_info_i = mcgrp.get_vio_time( ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                move_result.vio_load_delta = vio_load_delta;
                move_result.vio_time_delta = new_time_info_i.second - old_time_info_i.second;
                move_result.vio_time_custom_num_delta = new_time_info_i.first - old_time_info_i.first;
            }
            else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }
        else {
            // Different routes!
            move_result.num_affected_routes = 2;
            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;

            move_result.route_lens.push_back(ns.routes[i_route]->length + i_route_length_delta);
            move_result.route_lens.push_back(ns.routes[u_route]->length + u_route_length_delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load - load_delta);
            move_result.route_loads.push_back(ns.routes[u_route]->load + load_delta);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers - move_seq.size());
            move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers + move_seq.size());

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;


            move_result.considerable = true;

            move_result.route_time_tbl = new_time_tbl;


            if(ns.policy.has_rule(INFEASIBLE)){
                auto old_time_info_i = mcgrp.get_vio_time( ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
                auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

                move_result.vio_load_delta = vio_load_delta;

                move_result.vio_time_delta =
                    new_time_info_i.second + new_time_info_u.second - old_time_info_i.second - old_time_info_u.second;
                move_result.vio_time_custom_num_delta =
                    new_time_info_i.first + new_time_info_u.first - old_time_info_i.first - old_time_info_u.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }

    }

}

void XPostInsert::move(LocalSearch &ns, const MCGRPTW &mcgrp)
{
    DEBUG_PRINT("execute a " + to_string(length) + "-insert:postsert move");
    My_Assert(move_result.considerable,"Invalid predictions");

    //phase 0: extract move arguments
    const int u = move_result.move_arguments.at("insert_pos").front();
    const vector<int>& input_move_seq = move_result.move_arguments.at("input_seq");
    const vector<int>& output_move_seq = move_result.move_arguments.at("output_seq");

    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(all_of(input_move_seq.begin(),input_move_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");
    My_Assert(all_of(output_move_seq.begin(),output_move_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");


    // phase 1: update the route level info
    const int i_route = move_result.route_id[0];
    const int u_route = move_result.route_id[1];

    if(move_result.num_affected_routes == 2){
        int edge_num = mcgrp.count_edges(input_move_seq);

        ns.routes[i_route]->num_edges -= edge_num;
        ns.routes[u_route]->num_edges += edge_num;
    }

    if (move_result.num_affected_routes == 1) {
        // affect one route
        My_Assert(i_route == u_route,"error, wrong arguments");

        ns.routes[i_route]->length = move_result.route_lens[0];

        ns.routes[i_route]->time_table = move_result.route_time_tbl[0];
    }
    else if(move_result.num_affected_routes == 2){
        // affect two routes

        ns.routes[i_route]->length = move_result.route_lens[0];
        ns.routes[u_route]->length = move_result.route_lens[1];

        ns.routes[i_route]->load = move_result.route_loads[0];
        ns.routes[u_route]->load = move_result.route_loads[1];

        ns.routes[i_route]->num_customers = move_result.route_custs_num[0];
        ns.routes[u_route]->num_customers = move_result.route_custs_num[1];

        ns.routes[i_route]->time_table = move_result.route_time_tbl[0];
        ns.routes[u_route]->time_table = move_result.route_time_tbl[1];
    }
    else{
        My_Assert(false, "error, wrong arguments!");
    }


    // phase 3: update the sequence level info

    move_result.move_arguments["pre_input"] = {ns.solution[input_move_seq.front()]->pre->ID};

    // remove from the i route
    ns.solution[input_move_seq.front()]->pre->next = ns.solution[input_move_seq.back()]->next;
    ns.solution[input_move_seq.back()]->next->pre = ns.solution[input_move_seq.front()]->pre;

    // postsert
    if(input_move_seq == output_move_seq){
        ns.solution[input_move_seq.front()]->pre = ns.solution[u];
        ns.solution[input_move_seq.back()]->next = ns.solution[u]->next;
        ns.solution[u]->next->pre = ns.solution[input_move_seq.back()];
        ns.solution[u]->next = ns.solution[input_move_seq.front()];
    }
    else{

        for(auto input_task : input_move_seq){
            ns.solution[input_task]->clear();
        }

        ns.solution[output_move_seq.front()]->pre = ns.solution[u];
        ns.solution[output_move_seq.back()]->next = ns.solution[u]->next;

        ns.solution.connect_tasks(output_move_seq);

        ns.solution[u]->next->pre = ns.solution[output_move_seq.back()];
        ns.solution[u]->next = ns.solution[output_move_seq.front()];

    }

    for(auto task : output_move_seq){
        ns.solution[task]->route_id = u_route;
    }

    // phase 4: compress the empty routes
    if(ns.routes[i_route]->num_customers == 0){
        move_result.move_arguments["empty_route"] = {1};
        ns.compress_empty_route(i_route);
    }

    // phase 5: modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");

    update_score(ns);


    ns.trace(mcgrp);
    move_result.reset();
}


bool XPostInsert::search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

    MoveResult BestM;

    vector<int> chosen_seq = ns.get_successor_tasks(length, chosen_task);
    My_Assert(chosen_seq.size() > 0,"You cannot generate an empty sequence!");

    if (chosen_seq.size() != length) {
        DEBUG_PRINT("The length after chosen Task is not " + to_string(length) + " in " + to_string(length) + "-insert");
        return false;
    }

    int offset = 0;
    if(ns.policy.has_rule(TOLERANCE)){
        offset = mcgrp._rng.Randint(0,mcgrp.actual_task_num);
    }
    ns.create_search_neighborhood(mcgrp, chosen_seq,"predecessor", offset);

    int b = chosen_seq.front();

    for(auto neighbor_task : ns.search_space) {
        My_Assert(neighbor_task != b, "neighbor Task can't be itself!");

        if (neighbor_task != DUMMY) {


            //j can't be dummy and b can't be dummy neither here
            int j = neighbor_task;

            if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                // doesn't overlap
                if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                    if (ns.policy.has_rule(FIRST_ACCEPT)) {
                        move(ns, mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)) {
                        if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                            BestM = move_result;
                    }
                    else {
                        My_Assert(false, "Unknown accept rule!");
                    }
                }
            }
        }
        else {
            DEBUG_PRINT("Neighbor Task is dummy Task");
            //j is dummy here and b can't be dummy neither
            //each start and end location of each route will be considered
            //total 2 x route_nums cases

            int current_start = ns.solution.very_start->next->ID;

            while (ns.solution[current_start]->next != ns.solution.very_end) {


                // Consider the start location
                int j = current_start;

                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.solution[ns.routes[current_route]->end]->pre->ID;
                j = current_end;
                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Advance to next route's starting Task
                if(ns.solution[current_end]->next == ns.solution.very_end){
                    current_start = current_end;
                }else{
                    current_start = ns.solution[current_end]->next->next->ID;
                }
            }
        }
    }

    if (ns.policy.has_rule(FIRST_ACCEPT)) {
        DEBUG_PRINT("No actual move: First Accept Rule");
        return false;
    }
    else if (ns.policy.has_rule(BEST_ACCEPT)) {
        if (BestM.considerable == false) {
            DEBUG_PRINT("No actual move: Best Accept Rule");
            return false;
        }
        else {
            move_result = BestM;
            move(ns, mcgrp);
            return true;
        }
    }
    else {
        My_Assert(false, "Unknown accept rule");
    }

}

bool XPostInsert::update_score(LocalSearch &ns)
{
    if(move_result.delta >= 0) return false;

    double reward = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    const int u = move_result.move_arguments.at("insert_pos").front();

    // intra-sequence
    for(int idx = 1; idx < move_result.move_arguments.at("input_seq").size(); ++idx){
        ns.score_matrix[move_result.move_arguments.at("input_seq")[idx-1]][move_result.move_arguments.at("output_seq")[idx]] -= reward;
    }

    for(int idx = 1; idx < move_result.move_arguments.at("output_seq").size(); ++idx){
        ns.score_matrix[move_result.move_arguments.at("output_seq")[idx-1]][move_result.move_arguments.at("output_seq")[idx]] += reward;
    }

    // ...u-seq()...
    const int seq_front = move_result.move_arguments.at("output_seq").front();
    const int seq_back = move_result.move_arguments.at("output_seq").back();
    ns.score_matrix[u][seq_front] += reward;
    ns.score_matrix[seq_back][max(0,ns.solution[seq_back]->next->ID)] += reward;

    // old sequence
    if(!move_result.move_arguments["empty_route"].size()){
        int before_input = move_result.move_arguments.at("pre_input").front();
        int after_input = ns.solution[before_input]->next->ID;
        ns.score_matrix[max(0,before_input)][max(0,after_input)] -= reward;
    }

    return true;
}

pair<double, double>
XPostInsert::cost_delta(LocalSearch &ns, const MCGRPTW &mcgrp, const vector<int> &move_seq, const int u)
{
    const int v = max(ns.solution[u]->next->ID, 0);

    const int h = max(ns.solution[move_seq.front()]->pre->ID, 0);
    const int l = max(ns.solution[move_seq.back()]->next->ID, 0);

    const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[move_seq.front()].head_node];
    const double kl = mcgrp.min_cost[mcgrp.inst_tasks[move_seq.back()].tail_node][mcgrp.inst_tasks[l].head_node];
    const double hl = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[l].head_node];

    double i_route_length_delta = hl - (hi + kl);
    for (int cursor = 0; cursor < move_seq.size() - 1; cursor++) {
        i_route_length_delta -=
            mcgrp.min_cost[mcgrp.inst_tasks[move_seq[cursor]].tail_node][mcgrp.inst_tasks[move_seq[cursor
                + 1]].head_node];
        i_route_length_delta -= mcgrp.inst_tasks[move_seq[cursor]].serv_cost;
    }
    i_route_length_delta -= mcgrp.inst_tasks[move_seq.back()].serv_cost;


    const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
    double u_route_length_delta = -uv;

    for (auto cursor = 0; cursor < move_seq.size() - 1; cursor++) {
        u_route_length_delta +=
            mcgrp.min_cost[mcgrp.inst_tasks[move_seq[cursor]].tail_node][mcgrp.inst_tasks[move_seq[
                cursor + 1]].head_node];
        u_route_length_delta += mcgrp.inst_tasks[move_seq[cursor]].serv_cost;
    }
    u_route_length_delta += mcgrp.inst_tasks[move_seq.back()].serv_cost;

    //handle ends
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[move_seq.front()].head_node];
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[move_seq.back()].tail_node][mcgrp.inst_tasks[v].head_node];

    return {u_route_length_delta,i_route_length_delta};
}

/*------------------------------------------------------------------------*/


vector<vector<RouteInfo::TimeTable>>
XPreInsert::expected_time_table(LocalSearch &ns,
                                const MCGRPTW &mcgrp,
                                const vector<int> &move_sequence,
                                const int i,
                                bool allow_infeasible){
    vector<vector<RouteInfo::TimeTable>>
        res(2,vector<RouteInfo::TimeTable>({{-1, -1}}));

    const int i_route = ns.solution[i]->route_id;
    const int u_route = ns.solution[move_sequence.front()]->route_id;

    auto intermediate = mcgrp.forecast_time_table(ns.routes[u_route]->time_table,
                                                  move_sequence, "remove", -1 , allow_infeasible);

    if(!mcgrp.isTimeTableFeasible(intermediate)){
        return res;
    }

    vector<RouteInfo::TimeTable> final;
    if(u_route == i_route){
        final = mcgrp.forecast_time_table(intermediate,
                                          move_sequence, "insert_before", i, allow_infeasible);
        intermediate = final;
    }else{
        final = mcgrp.forecast_time_table(ns.routes[i_route]->time_table,
                                          move_sequence, "insert_before", i, allow_infeasible);
    }

    if(!mcgrp.isTimeTableFeasible(final)){
        return res;
    }

    // 0 for u_route, 1 for i_route
    res[0] = intermediate;
    res[1] = final;

    return res;
}


bool
XPreInsert::considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, const vector<int>& move_seq, const int u)
{
    // This move used to move the `move_seq` to the front of the task `u`
    // Two policy are contained, 1: no mode choice 2: best mode choice

    // Task u cannot be dummy Task
    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong Task");
    My_Assert(all_of(move_seq.begin(), move_seq.end(),
                     [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}), "Wrong Task");

    call_times++;
    move_result.reset();

    if(ns.solution[move_seq.back()]->next->ID == u){
        // already at the desired position
        return false;
    }

    const int i_route = ns.solution[move_seq.front()]->route_id;
    const int u_route = ns.solution[u]->route_id;

    // load legality is independent with mode choice
    int load_delta = 0;
    double vio_load_delta = 0;
    if (i_route != u_route) {

        for (auto task : move_seq) {
            load_delta += mcgrp.inst_tasks[task].demand;
        }

        if (ns.policy.has_rule(FEASIBLE) && ns.routes[u_route]->load + load_delta > mcgrp.capacity) {
                return false;
        }
        else if (ns.policy.has_rule(INFEASIBLE)) {

            int pseudo_capacity = ns.policy.get_pseudo_capacity(mcgrp.capacity);
            if (ns.routes[u_route]->load + load_delta > pseudo_capacity) {
                return false;
            }

            //u_route vio-load calculate
            if (ns.routes[u_route]->load + load_delta > mcgrp.capacity) {
                //if insert Task to route u and over load
                if (ns.routes[u_route]->load >= mcgrp.capacity) {
                    //if the route u already over loaded
                    vio_load_delta += load_delta;
                }
                else {
                    vio_load_delta += ns.routes[u_route]->load + load_delta - mcgrp.capacity;
                }
            }

            //i_route vio-load calculate
            if (ns.routes[i_route]->load > mcgrp.capacity) {
                //if remove Task from route i and over load
                if (ns.routes[i_route]->load - load_delta >= mcgrp.capacity) {
                    //if still over loaded
                    vio_load_delta -= load_delta;
                }
                else {
                    vio_load_delta -= (ns.routes[i_route]->load - mcgrp.capacity);
                }
            }
        }
    }


    if(ns.policy.has_rule(BEST_ACCEPT) && mcgrp.count_edges(move_seq) > 0) {
        // with mode refine

        vector<int> first_seq;
        vector<int> last_seq;

        if(i_route == u_route){
            // special case
            vector<int> u_route_seq = ns.get_sub_seq(ns.routes[ns.solution[u]->route_id]->start, ns.routes[ns.solution[u]->route_id]->end);

            // try to use offset to generate sequence
            int move_seq_loc = -1;
            int u_loc = -1;
            for(int ii = 0;ii < u_route_seq.size(); ii++){
                if(u_route_seq[ii] == u) u_loc = ii;
                if(u_route_seq[ii] == move_seq.front()) move_seq_loc = ii;
            }

            if(move_seq_loc < u_loc){
                first_seq.insert (first_seq.end(),u_route_seq.begin(),u_route_seq.begin() + move_seq_loc);
                first_seq.insert (first_seq.end(), u_route_seq.begin() + move_seq_loc + move_seq.size(), u_route_seq.begin() + u_loc);

                last_seq.insert (last_seq.end(), u_route_seq.begin() + u_loc, u_route_seq.end());
            }
            else if(move_seq_loc > u_loc){
                first_seq.insert (first_seq.end(),u_route_seq.begin(),u_route_seq.begin() + u_loc);

                last_seq.insert (last_seq.end(), u_route_seq.begin() + u_loc,u_route_seq.begin() + move_seq_loc);
                last_seq.insert( last_seq.end(), u_route_seq.begin() + move_seq_loc + move_seq.size(), u_route_seq.end());
            }else{
                My_Assert(false,"error, wrong location");
            }

        }
        else{
            // trivial case
            first_seq = ns.get_sub_seq(ns.routes[ns.solution[u]->route_id]->start, ns.solution[u]->pre->ID);
            last_seq = ns.get_sub_seq(u, ns.routes[ns.solution[u]->route_id]->end);
        }

        My_Assert(first_seq.front() < 0 && last_seq.back() < 0, "error, wrong sequence!");

        viterbi::PseudoTask first_task = Seq2Pseudo(mcgrp,first_seq, 1);
        viterbi::PseudoTask last_task = Seq2Pseudo(mcgrp,last_seq, 2);

        // insert sequence will invoke viterbi result
        auto best_sequence = viterbi::viterbi_decode(&first_task, &last_task,move_seq,mcgrp,ns.policy);

        // use best sequence to update the move result
        if(best_sequence.cost == INT32_MAX){
            return false;
        }

        My_Assert(best_sequence.seq.front() == INT32_MAX && best_sequence.seq.back() == INT32_MAX, "error, wrong sequence");

        move_result.choose_tasks(move_seq.front(), u);

        move_result.move_arguments["input_seq"] = move_seq;
        move_result.move_arguments["output_seq"] = vector<int>(best_sequence.seq.begin() + 1, best_sequence.seq.end() - 1);
        move_result.move_arguments["insert_pos"] = {u};

        double u_route_length_delta = best_sequence.cost - ns.routes[u_route]->length;
        double i_route_length_delta = 0;
        if(i_route != u_route){
            auto routes_cost_delta = cost_delta(ns,mcgrp,move_seq,u);
            i_route_length_delta = routes_cost_delta.second;
        }

        if (i_route == u_route) {

            move_result.num_affected_routes = 1;

            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;
            move_result.route_lens.push_back(ns.routes[i_route]->length + move_result.delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

            move_result.considerable = true;

            vector<int> u_route_seq;
            u_route_seq.insert(u_route_seq.end(), first_seq.begin() + 1, first_seq.end());
            u_route_seq.insert(u_route_seq.end(), move_result.move_arguments["output_seq"].begin(), move_result.move_arguments["output_seq"].end());
            u_route_seq.insert(u_route_seq.end(), last_seq.begin(), last_seq.end() - 1);

            move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(u_route_seq, mcgrp.cal_arrive_time(u_route_seq)));


            if(ns.policy.has_rule(INFEASIBLE)){
                auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                move_result.vio_load_delta = vio_load_delta;
                move_result.vio_time_delta = new_time_info_i.second - old_time_info_i.second;
                move_result.vio_time_custom_num_delta = new_time_info_i.first - old_time_info_i.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }
        else{

            // Different routes!
            move_result.num_affected_routes = 2;
            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;

            move_result.route_lens.push_back(ns.routes[i_route]->length + i_route_length_delta);
            move_result.route_lens.push_back(ns.routes[u_route]->length + u_route_length_delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load - load_delta);
            move_result.route_loads.push_back(ns.routes[u_route]->load + load_delta);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers - move_seq.size());
            move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers + move_seq.size());

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

            move_result.considerable = true;

            vector<int> i_route_seq;
            for(int ii = 0; ii < ns.routes[i_route]->time_table.size();){
                if(move_seq.front() == ns.routes[i_route]->time_table[ii].task){
                    ii += move_seq.size();
                    continue;
                }

                i_route_seq.push_back(ns.routes[i_route]->time_table[ii].task);
                ii ++;
            }
            move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(i_route_seq, mcgrp.cal_arrive_time(i_route_seq)));

            vector<int> u_route_seq;
            u_route_seq.insert(u_route_seq.end(), first_seq.begin() + 1, first_seq.end());
            u_route_seq.insert(u_route_seq.end(), move_result.move_arguments["output_seq"].begin(), move_result.move_arguments["output_seq"].end());
            u_route_seq.insert(u_route_seq.end(), last_seq.begin(), last_seq.end() - 1);

            move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(u_route_seq, mcgrp.cal_arrive_time(u_route_seq)));

            if(ns.policy.has_rule(INFEASIBLE)){
                move_result.vio_load_delta = vio_load_delta;

                auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
                auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

                move_result.vio_time_delta =
                    new_time_info_i.second + new_time_info_u.second - old_time_info_i.second - old_time_info_u.second;
                move_result.vio_time_custom_num_delta =
                    new_time_info_i.first + new_time_info_u.first - old_time_info_i.first - old_time_info_u.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;

        }


    }
    else {
        // without mode refine

        bool allow_infeasible = ns.policy.has_rule(INFEASIBLE) ? true : false;
        vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};

        new_time_tbl = expected_time_table(ns, mcgrp, move_seq, u, allow_infeasible);

        if(!mcgrp.isTimeTableFeasible(new_time_tbl[0])){
            return false;
        }

        if(allow_infeasible){
            if(ns.policy.check_time_window(mcgrp,new_time_tbl)){
                return false;
            }
        }

        auto routes_cost_delta = cost_delta(ns,mcgrp,move_seq,u);
        double u_route_length_delta = routes_cost_delta.first;
        double i_route_length_delta = routes_cost_delta.second;

        move_result.choose_tasks(move_seq.front(), u);

        move_result.move_arguments["input_seq"] = move_seq;
        move_result.move_arguments["output_seq"] = move_seq;
        move_result.move_arguments["insert_pos"] = {u};


        if (i_route == u_route) {

            move_result.num_affected_routes = 1;

            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;
            move_result.route_lens.push_back(ns.routes[i_route]->length + move_result.delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;


            move_result.considerable = true;

            move_result.route_time_tbl.emplace_back(new_time_tbl[0]);


            if(ns.policy.has_rule(INFEASIBLE)){
                auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                move_result.vio_load_delta = vio_load_delta;
                move_result.vio_time_delta = new_time_info_i.second - old_time_info_i.second;
                move_result.vio_time_custom_num_delta = new_time_info_i.first - old_time_info_i.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }
        else {
            // Different routes!
            move_result.num_affected_routes = 2;
            move_result.route_id.push_back(i_route);
            move_result.route_id.push_back(u_route);

            move_result.delta = u_route_length_delta + i_route_length_delta;

            move_result.route_lens.push_back(ns.routes[i_route]->length + i_route_length_delta);
            move_result.route_lens.push_back(ns.routes[u_route]->length + u_route_length_delta);

            move_result.route_loads.push_back(ns.routes[i_route]->load - load_delta);
            move_result.route_loads.push_back(ns.routes[u_route]->load + load_delta);

            move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers - move_seq.size());
            move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers + move_seq.size());

            move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

            move_result.considerable = true;

            move_result.route_time_tbl = new_time_tbl;

            if(ns.policy.has_rule(INFEASIBLE)){
                move_result.vio_load_delta = vio_load_delta;

                auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
                auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

                auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
                auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

                move_result.vio_time_delta =
                    new_time_info_i.second + new_time_info_u.second - old_time_info_i.second - old_time_info_u.second;
                move_result.vio_time_custom_num_delta =
                    new_time_info_i.first + new_time_info_u.first - old_time_info_i.first - old_time_info_u.first;
            }else{
                move_result.vio_load_delta = 0;
                move_result.vio_time_delta = 0;
                move_result.vio_time_custom_num_delta = 0;
            }

            return true;
        }

    }

}


void XPreInsert::move(LocalSearch &ns, const MCGRPTW &mcgrp)
{

    success_times++;
    DEBUG_PRINT("execute a " + to_string(length) + "-insert:presert move");
    My_Assert(move_result.considerable,"Invalid predictions");


    //phase 0: extract move arguments
    const int u = move_result.move_arguments.at("insert_pos").front();
    const vector<int>& input_move_seq = move_result.move_arguments.at("input_seq");
    const vector<int>& output_move_seq = move_result.move_arguments.at("output_seq");

    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(all_of(input_move_seq.begin(),input_move_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");
    My_Assert(all_of(output_move_seq.begin(),output_move_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");


    // phase 1: update the route level info
    const int i_route = move_result.route_id[0];
    const int u_route = move_result.route_id[1];

    if(move_result.num_affected_routes == 2){
        int edge_num = mcgrp.count_edges(input_move_seq);

        ns.routes[i_route]->num_edges -= edge_num;
        ns.routes[u_route]->num_edges += edge_num;
    }

    if (move_result.num_affected_routes == 1) {
        // affect one route
        My_Assert(i_route == u_route,"error, wrong arguments");

        ns.routes[i_route]->length = move_result.route_lens[0];

        ns.routes[i_route]->time_table = move_result.route_time_tbl[0];

    }else if(move_result.num_affected_routes == 2){
        // affect two routes

        ns.routes[i_route]->length = move_result.route_lens[0];
        ns.routes[u_route]->length = move_result.route_lens[1];

        ns.routes[i_route]->load = move_result.route_loads[0];
        ns.routes[u_route]->load = move_result.route_loads[1];

        ns.routes[i_route]->num_customers = move_result.route_custs_num[0];
        ns.routes[u_route]->num_customers = move_result.route_custs_num[1];

        ns.routes[i_route]->time_table = move_result.route_time_tbl[0];
        ns.routes[u_route]->time_table = move_result.route_time_tbl[1];


    }else{
        My_Assert(false, "error, wrong arguments!");
    }

    // phase 3: update the sequence level info

    // remove from the i route
    ns.solution[input_move_seq.front()]->pre->next = ns.solution[input_move_seq.back()]->next;
    ns.solution[input_move_seq.back()]->next->pre = ns.solution[input_move_seq.front()]->pre;

    // presert
    move_result.move_arguments["pre_input"] = {ns.solution[input_move_seq.front()]->pre->ID};
    if(input_move_seq == output_move_seq){

        ns.solution[input_move_seq.front()]->pre = ns.solution[u]->pre;
        ns.solution[input_move_seq.back()]->next = ns.solution[u];
        ns.solution[u]->pre->next = ns.solution[input_move_seq.front()];
        ns.solution[u]->pre = ns.solution[input_move_seq.back()];
    }else{

        for(auto input_task : input_move_seq){
            ns.solution[input_task]->clear();
        }

        ns.solution[output_move_seq.front()]->pre = ns.solution[u]->pre;
        ns.solution[output_move_seq.back()]->next = ns.solution[u];

        ns.solution.connect_tasks(output_move_seq);

        ns.solution[u]->pre->next = ns.solution[output_move_seq.front()];
        ns.solution[u]->pre = ns.solution[output_move_seq.back()];
    }

    for(auto task : output_move_seq){
        ns.solution[task]->route_id = u_route;
    }

    // phase 4: compress the empty routes
    if(ns.routes[i_route]->num_customers == 0){
        move_result.move_arguments["empty_route"] = {1};
        ns.compress_empty_route(i_route);
    }

    // phase 5: modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");

    update_score(ns);

    ns.trace(mcgrp);
    move_result.reset();
}


bool XPreInsert::search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

    MoveResult BestM;

    vector<int> chosen_seq = ns.get_successor_tasks(length, chosen_task);
    My_Assert(chosen_seq.size() > 0,"You cannot generate an empty sequence!");

    if (chosen_seq.size() != length) {
        DEBUG_PRINT("warning, the length after chosen task is not " + to_string(length) + " in " + to_string(length) + "-insert");
        return false;
    }

    int offset = 0;
    if(ns.policy.has_rule(TOLERANCE)){
        offset = mcgrp._rng.Randint(0, mcgrp.actual_task_num);
    }
    ns.create_search_neighborhood(mcgrp, chosen_seq,"successor", offset);


    int b = chosen_seq.front();

    for(auto neighbor_task : ns.search_space) {
        My_Assert(neighbor_task != b, "neighbor Task can't be itself!");

        if (neighbor_task != DUMMY) {


            //j can't be dummy and b can't be dummy neither here
            int j = neighbor_task;

            if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                // doesn't overlap
                if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                    if (ns.policy.has_rule(FIRST_ACCEPT)) {
                        move(ns, mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)) {
                        if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                            BestM = move_result;
                    }
                    else {
                        My_Assert(false, "Unknown accept rule!");
                    }
                }
            }
        }
        else {
            DEBUG_PRINT("Neighbor Task is dummy Task");
            //j is dummy here and b can't be dummy neither
            //each start and end location of each route will be considered
            //total 2 x route_nums cases

            int current_start = ns.solution.very_start->next->ID;

            while (ns.solution[current_start]->next != ns.solution.very_end) {

                // Consider the start location
                int j = current_start;

                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.solution[ns.routes[current_route]->end]->pre->ID;
                j = current_end;
                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Advance to next route's starting Task
                if(ns.solution[current_end]->next == ns.solution.very_end){
                    current_start = current_end;
                }else{
                    current_start = ns.solution[current_end]->next->next->ID;
                }
            }
        }
    }

    if (ns.policy.has_rule(FIRST_ACCEPT)) {
        DEBUG_PRINT("No actual move: First Accept Rule");
        return false;
    }
    else if (ns.policy.has_rule(BEST_ACCEPT)) {
        if (BestM.considerable == false) {
            DEBUG_PRINT("No actual move: Best Accept Rule");
            return false;
        }
        else {
            move_result = BestM;
            move(ns, mcgrp);
            return true;
        }
    }
    else {
        My_Assert(false, "Unknown accept rule");
    }

}

bool XPreInsert::update_score(LocalSearch &ns)
{
    if(move_result.delta >= 0) return false;

    double reward = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    const int u = move_result.move_arguments.at("insert_pos").front();

    // intra-sequence
    for(int idx = 1; idx < move_result.move_arguments.at("input_seq").size(); ++idx){
        ns.score_matrix[move_result.move_arguments.at("input_seq")[idx-1]][move_result.move_arguments.at("output_seq")[idx]] -= reward;
    }

    for(int idx = 1; idx < move_result.move_arguments.at("output_seq").size(); ++idx){
        ns.score_matrix[move_result.move_arguments.at("output_seq")[idx-1]][move_result.move_arguments.at("output_seq")[idx]] += reward;
    }

    // ...seq()-u...
    int seq_front = move_result.move_arguments.at("output_seq").front();
    int seq_back = move_result.move_arguments.at("output_seq").back();
    ns.score_matrix[seq_back][u] += reward;
    ns.score_matrix[max(0,ns.solution[seq_front]->pre->ID)][seq_front] += reward;

    // old sequence
    if(!move_result.move_arguments["empty_route"].size()){
        int before_input = move_result.move_arguments.at("pre_input").front();
        int after_input = ns.solution[before_input]->next->ID;
        ns.score_matrix[max(0,before_input)][max(0,after_input)] -= reward;
    }

    return true;
}

pair<double,double> XPreInsert::cost_delta(LocalSearch &ns, const MCGRPTW &mcgrp, const vector<int> &move_seq, const int u)
{
    const int t = max(ns.solution[u]->pre->ID, 0);

    const int h = max(ns.solution[move_seq.front()]->pre->ID, 0);
    const int l = max(ns.solution[move_seq.back()]->next->ID, 0);

    const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[move_seq.front()].head_node];
    const double kl = mcgrp.min_cost[mcgrp.inst_tasks[move_seq.back()].tail_node][mcgrp.inst_tasks[l].head_node];
    const double hl = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[l].head_node];

    double i_route_length_delta = hl - (hi + kl);

    for (int cursor = 0; cursor < move_seq.size() - 1; cursor++) {
        i_route_length_delta -=
            mcgrp.min_cost[mcgrp.inst_tasks[move_seq[cursor]].tail_node][mcgrp.inst_tasks[move_seq[cursor
                + 1]].head_node];
        i_route_length_delta -= mcgrp.inst_tasks[move_seq[cursor]].serv_cost;
    }
    i_route_length_delta -= mcgrp.inst_tasks[move_seq.back()].serv_cost;


    const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
    double u_route_length_delta = -tu;

    for (auto cursor = 0; cursor < move_seq.size() - 1; cursor++) {
        u_route_length_delta +=
            mcgrp.min_cost[mcgrp.inst_tasks[move_seq[cursor]].tail_node][mcgrp.inst_tasks[move_seq[
                cursor + 1]].head_node];
        u_route_length_delta += mcgrp.inst_tasks[move_seq[cursor]].serv_cost;
    }

    u_route_length_delta += mcgrp.inst_tasks[move_seq.back()].serv_cost;
    //handle ends
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[move_seq.front()].head_node];
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[move_seq.back()].tail_node][mcgrp.inst_tasks[u].head_node];

    return {u_route_length_delta,i_route_length_delta};
}
