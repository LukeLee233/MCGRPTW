//
// Created by luke on 2020/11/5.
//

#include "insert.h"

vector<vector<RouteInfo::TimeTable>> XPostInsert::expected_time_table(HighSpeedNeighBorSearch &ns,
                                                                      const MCGRP &mcgrp,
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

bool XPostInsert::considerable_move(HighSpeedNeighBorSearch &ns,
                                    const MCGRP &mcgrp,
                                    vector<int> disturbance_seq,
                                    const int u)
{
    // Task u cannot be dummy Task
    My_Assert(u>=1 && u<=mcgrp.actual_task_num,"Wrong Task");
    My_Assert(all_of(disturbance_seq.begin(),disturbance_seq.end(),
                     [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong Task");

    if(u == ns.solution[disturbance_seq.front()]->pre->ID){
        // Nothing to do
        move_result.reset();
        return false;
    }

    const int i_route = ns.solution[disturbance_seq.front()]->route_id;
    const int u_route = ns.solution[u]->route_id;

    // Check load feasibility easily
    int load_delta = 0;
    double vio_load_delta = 0;
    //not the same route
    if (i_route != u_route) {
        for (auto task : disturbance_seq) {
            load_delta += mcgrp.inst_tasks[task].demand;
        }

        if (ns.policy.has_rule(DELTA_ONLY)) {
            if (ns.routes[u_route]->load + load_delta > mcgrp.capacity) {
                move_result.reset();
                return false;
            }
        }
        else if (ns.policy.has_rule(FITNESS_ONLY)) {
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

    bool allow_infeasible = ns.policy.has_rule(FITNESS_ONLY) ? true : false;
    vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,disturbance_seq,u,allow_infeasible);

    if(!mcgrp.isTimeTableFeasible(new_time_tbl[0])){
        move_result.reset();
        return false;
    }

    const int v = max(ns.solution[u]->next->ID, 0);

    const int h = max(ns.solution[disturbance_seq.front()]->pre->ID, 0);
    const int l = max(ns.solution[disturbance_seq.back()]->next->ID, 0);

    const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[disturbance_seq.front()].head_node];
    const double kl = mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq.back()].tail_node][mcgrp.inst_tasks[l].head_node];
    const double hl = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[l].head_node];

    double i_route_length_delta = hl - (hi + kl);
    for (int cursor = 0; cursor < disturbance_seq.size() - 1; cursor++) {
        i_route_length_delta -=
            mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq[cursor]].tail_node][mcgrp.inst_tasks[disturbance_seq[cursor
                + 1]].head_node];
        i_route_length_delta -= mcgrp.inst_tasks[disturbance_seq[cursor]].serv_cost;
    }
    i_route_length_delta -= mcgrp.inst_tasks[disturbance_seq.back()].serv_cost;


    const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
    double u_route_length_delta = -uv;

    for (auto cursor = 0; cursor < disturbance_seq.size() - 1; cursor++) {
        u_route_length_delta +=
            mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq[cursor]].tail_node][mcgrp.inst_tasks[disturbance_seq[
                cursor + 1]].head_node];
        u_route_length_delta += mcgrp.inst_tasks[disturbance_seq[cursor]].serv_cost;
    }
    u_route_length_delta += mcgrp.inst_tasks[disturbance_seq.back()].serv_cost;
    //handle ends
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[disturbance_seq.front()].head_node];
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq.back()].tail_node][mcgrp.inst_tasks[v].head_node];


    move_result.choose_tasks(disturbance_seq.front(), u);

    move_result.move_arguments = disturbance_seq;
    move_result.move_arguments.push_back(u);

    // Check if we need to reduce the # of routes here
    const int start_i = ns.routes[i_route]->start;
    const int end_i = ns.routes[i_route]->end;
    if (start_i == disturbance_seq.front() && end_i == disturbance_seq.back())
        move_result.total_number_of_routes = ns.routes.activated_route_id.size() - 1;
    else
        move_result.total_number_of_routes = ns.routes.activated_route_id.size();


    if (i_route == u_route) {

        move_result.num_affected_routes = 1;

        move_result.route_id.push_back(i_route);

        move_result.delta = u_route_length_delta + i_route_length_delta;
        move_result.route_lens.push_back(ns.routes[i_route]->length + move_result.delta);

        move_result.route_loads.push_back(ns.routes[i_route]->load);

        move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.vio_load_delta = vio_load_delta;

        move_result.considerable = true;

        move_result.route_time_tbl.emplace_back(new_time_tbl[0]);
        move_result.vio_time_delta =
            mcgrp.get_vio_time(move_result.route_time_tbl[0])
                - mcgrp.get_vio_time(ns.routes[i_route]->time_table);
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

        move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers - disturbance_seq.size());
        move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers + disturbance_seq.size());

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.vio_load_delta = vio_load_delta;

        move_result.considerable = true;

        move_result.route_time_tbl = new_time_tbl;
        move_result.vio_time_delta =
            mcgrp.get_vio_time(move_result.route_time_tbl[0])
                + mcgrp.get_vio_time(move_result.route_time_tbl[1])
                - mcgrp.get_vio_time(ns.routes[i_route]->time_table)
                - mcgrp.get_vio_time(ns.routes[u_route]->time_table);
        return true;
    }

    My_Assert(false,"Cannot reach here!");
}

void XPostInsert::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{

#ifdef DEBUG
    hit_count++;
#endif

    DEBUG_PRINT("execute a " + to_string(length) + "-insert:postsert move");
    My_Assert(move_result.considerable,"Invalid predictions");

    //extract move arguments
    const int u = move_result.move_arguments.back();
    vector<int> disturbance_sequence(move_result.move_arguments.begin(),move_result.move_arguments.end() - 1);

    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(all_of(disturbance_sequence.begin(),disturbance_sequence.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");

    int edge_num = mcgrp.count_edges(disturbance_sequence);
    if(edge_num > 0 && move_result.num_affected_routes == 2){
        const int i_route = move_result.route_id[0];
        const int u_route = move_result.route_id[1];

        ns.routes[i_route]->num_edges -= edge_num;
        ns.routes[u_route]->num_edges += edge_num;
    }

    if (move_result.total_number_of_routes != ns.routes.activated_route_id.size()) {
        //need to eliminate an empty route
        My_Assert(move_result.total_number_of_routes == ns.routes.activated_route_id.size() - 1,"Incorrect routes number");
        My_Assert(move_result.num_affected_routes == 2,"Incorrect routes number");
        const int i_route = move_result.route_id[0];
        const int u_route = move_result.route_id[1];
        My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");
        My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");
        My_Assert(ns.routes[i_route]->start == disturbance_sequence.front() && ns.routes[i_route]->end == disturbance_sequence.back(),"Wrong route");

        //Modify routes info
        ns.routes[u_route]->length = move_result.route_lens[1];

        ns.routes[u_route]->load = move_result.route_loads[1];

        ns.routes[u_route]->num_customers = move_result.route_custs_num[1];

        ns.routes[u_route]->time_table = move_result.route_time_tbl[1];

        for(auto task:disturbance_sequence){
            ns.solution[task]->route_id = u_route;
        }

        if(ns.solution[u]->next->ID < 0){
            ns.routes[u_route]->end = disturbance_sequence.back();
        }

        ns.routes.free_route(i_route);

        //handle solution
        //record next dummy Task of actual u
        int dummy_marker = 0;
        My_Assert(ns.solution[disturbance_sequence.front()]->pre->ID < 0 && ns.solution[disturbance_sequence.back()]->next->ID < 0,"Route is not empty!");
        if(ns.solution[disturbance_sequence.back()]->next == ns.solution.very_end){
            dummy_marker = ns.solution[disturbance_sequence.front()]->pre->ID;
        }
        else{
            dummy_marker = ns.solution[disturbance_sequence.back()]->next->ID;
        }

        //handle solution
        //remove
        ns.solution[disturbance_sequence.front()]->pre->next = ns.solution[disturbance_sequence.back()]->next;
        ns.solution[disturbance_sequence.back()]->next->pre = ns.solution[disturbance_sequence.front()]->pre;

        //postsert
        ns.solution[disturbance_sequence.front()]->pre = ns.solution[u];
        ns.solution[disturbance_sequence.back()]->next = ns.solution[u]->next;
        ns.solution[u]->next->pre = ns.solution[disturbance_sequence.back()];
        ns.solution[u]->next = ns.solution[disturbance_sequence.front()];

        //free dummy marker
        ns.solution[dummy_marker]->pre->next = ns.solution[dummy_marker]->next;
        ns.solution[dummy_marker]->next->pre = ns.solution[dummy_marker]->pre;
        ns.solution.dummypool.free_dummy(dummy_marker);
    }
    else {
        //no routes eliminates
        //Modify routes info
        if (move_result.num_affected_routes == 1) {
            //No need to change route id
            const int route_id = move_result.route_id[0];
            My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[route_id]->length = move_result.route_lens[0];

            ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

            if(ns.solution[disturbance_sequence.front()]->pre->ID < 0){
                ns.routes[route_id]->start = ns.solution[disturbance_sequence.back()]->next->ID;
            }

            if(ns.solution[disturbance_sequence.back()]->next->ID < 0){
                ns.routes[route_id]->end = ns.solution[disturbance_sequence.front()]->pre->ID;
            }

            if(ns.solution[u]->next->ID < 0){
                ns.routes[route_id]->end = disturbance_sequence.back();
            }
        }
        else{
            const int i_route = move_result.route_id[0];
            const int u_route = move_result.route_id[1];

            My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[i_route]->length = move_result.route_lens[0];
            ns.routes[u_route]->length = move_result.route_lens[1];

            ns.routes[i_route]->load = move_result.route_loads[0];
            ns.routes[u_route]->load = move_result.route_loads[1];

            ns.routes[i_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[u_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[i_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[u_route]->time_table = move_result.route_time_tbl[1];

            for(auto task:disturbance_sequence){
                ns.solution[task]->route_id = u_route;
            }

            if(ns.solution[disturbance_sequence.front()]->pre->ID < 0){
                ns.routes[i_route]->start = ns.solution[disturbance_sequence.back()]->next->ID;
            }

            if(ns.solution[disturbance_sequence.back()]->next->ID < 0){
                ns.routes[i_route]->end = ns.solution[disturbance_sequence.front()]->pre->ID;
            }

            if(ns.solution[u]->next->ID < 0){
                ns.routes[u_route]->end = disturbance_sequence.back();
            }
        }

        //handle solution
        //remove
        ns.solution[disturbance_sequence.front()]->pre->next = ns.solution[disturbance_sequence.back()]->next;
        ns.solution[disturbance_sequence.back()]->next->pre = ns.solution[disturbance_sequence.front()]->pre;

        //postsert
        ns.solution[disturbance_sequence.front()]->pre = ns.solution[u];
        ns.solution[disturbance_sequence.back()]->next = ns.solution[u]->next;
        ns.solution[u]->next->pre = ns.solution[disturbance_sequence.back()];
        ns.solution[u]->next = ns.solution[disturbance_sequence.front()];

    }


    //modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");


    if(move_result.delta == 0){
        ns.equal_step++;
    }

    update_score(ns);

    ns.trace(mcgrp);
    move_result.reset();
    ns.search_step++;
}


bool XPostInsert::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

#ifdef DEBUG
    attempt_count = 0;
    hit_count = 0;
#endif

    MoveResult BestM;

    vector<int> chosen_seq = ns.get_successor_tasks(length, chosen_task);

    int offset = 0;
    if(ns.policy.has_rule(TOLERANCE)){
        offset = rand() % mcgrp.actual_task_num;
    }
    ns.create_search_neighborhood(mcgrp, chosen_seq,"predecessor", offset);


    My_Assert(chosen_seq.size()>0,"You cannot generate an empty sequence!");

    if (chosen_seq.size() != length) {
        DEBUG_PRINT("The length after chosen Task is not " + to_string(length) + " in " + to_string(length) + "-insert");
        return false;
    }

    int b = chosen_seq.front();

    for(auto neighbor_task : ns.search_space) {
        My_Assert(neighbor_task != b, "neighbor Task can't be itself!");

        if (neighbor_task != DUMMY) {

#ifdef DEBUG
            attempt_count++;
#endif

            //j can't be dummy and b can't be dummy neither here
            int j = neighbor_task;

            if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                // doesn't overlap
                if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(move_result)) {
                    if (ns.policy.has_rule(FIRST_ACCEPT)) {
                        move(ns, mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)) {
                        if (ns.policy.check_result(move_result, BestM))
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

#ifdef DEBUG
                attempt_count += 2;
#endif

                // Consider the start location
                int j = current_start;

                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.routes[current_route]->end;
                j = current_end;
                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Advance to next route's starting Task
                current_start = ns.solution[current_end]->next->next->ID;
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

bool XPostInsert::update_score(HighSpeedNeighBorSearch &ns)
{
    if(move_result.delta > 0) return false;

    const int u = move_result.move_arguments.back();
    const int seq_front = *(move_result.move_arguments.begin());
    const int seq_back = *(move_result.move_arguments.end()-2);
    double penalty = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    // ...u-seq()...
    ns.score_matrix[u][seq_front] += penalty;
    ns.score_matrix[seq_back][max(0,ns.solution[seq_back]->next->ID)] += penalty;

    return true;
}

/*------------------------------------------------------------------------*/


vector<vector<RouteInfo::TimeTable>>
XPreInsert::expected_time_table(HighSpeedNeighBorSearch &ns,
                                const MCGRP &mcgrp,
                                vector<int> &disturbance_seq,
                                const int i,
                                bool allow_infeasible){
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
                                          disturbance_seq,"insert_before",i,allow_infeasible);
        intermediate = final;
    }else{
        final = mcgrp.forecast_time_table(ns.routes[i_route]->time_table,
                                          disturbance_seq,"insert_before",i,allow_infeasible);
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
XPreInsert::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, vector<int> disturbance_seq, const int u)
{
    // Task u cannot be dummy Task
    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong Task");
    My_Assert(all_of(disturbance_seq.begin(), disturbance_seq.end(),
                     [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}), "Wrong Task");

    if(u == ns.solution[disturbance_seq.back()]->next->ID){
        // Nothing to do
        move_result.reset();
        return false;
    }

    const int i_route = ns.solution[disturbance_seq.front()]->route_id;
    const int u_route = ns.solution[u]->route_id;

    // Check load feasibility easily
    int load_delta = 0;
    double vio_load_delta = 0;
    //not the same route
    if (i_route != u_route) {
        for (auto task : disturbance_seq) {
            load_delta += mcgrp.inst_tasks[task].demand;
        }

        if (ns.policy.has_rule(DELTA_ONLY)) {
            if (ns.routes[u_route]->load + load_delta > mcgrp.capacity) {
                move_result.reset();
                return false;
            }
        }
        else if (ns.policy.has_rule(FITNESS_ONLY)) {
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

    bool allow_infeasible = ns.policy.has_rule(FITNESS_ONLY) ? true : false;
    vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,disturbance_seq,u,allow_infeasible);

    if(!mcgrp.isTimeTableFeasible(new_time_tbl[0])){
        move_result.reset();
        return false;
    }

    const int t = max(ns.solution[u]->pre->ID, 0);

    const int h = max(ns.solution[disturbance_seq.front()]->pre->ID, 0);
    const int l = max(ns.solution[disturbance_seq.back()]->next->ID, 0);

    const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[disturbance_seq.front()].head_node];
    const double kl = mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq.back()].tail_node][mcgrp.inst_tasks[l].head_node];
    const double hl = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[l].head_node];

    double i_route_length_delta = hl - (hi + kl);

    for (int cursor = 0; cursor < disturbance_seq.size() - 1; cursor++) {
        i_route_length_delta -=
            mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq[cursor]].tail_node][mcgrp.inst_tasks[disturbance_seq[cursor
                + 1]].head_node];
        i_route_length_delta -= mcgrp.inst_tasks[disturbance_seq[cursor]].serv_cost;
    }
    i_route_length_delta -= mcgrp.inst_tasks[disturbance_seq.back()].serv_cost;


    const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
    double u_route_length_delta = -tu;

    for (auto cursor = 0; cursor < disturbance_seq.size() - 1; cursor++) {
        u_route_length_delta +=
            mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq[cursor]].tail_node][mcgrp.inst_tasks[disturbance_seq[
                cursor + 1]].head_node];
        u_route_length_delta += mcgrp.inst_tasks[disturbance_seq[cursor]].serv_cost;
    }

    u_route_length_delta += mcgrp.inst_tasks[disturbance_seq.back()].serv_cost;
    //handle ends
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[disturbance_seq.front()].head_node];
    u_route_length_delta +=
        mcgrp.min_cost[mcgrp.inst_tasks[disturbance_seq.back()].tail_node][mcgrp.inst_tasks[u].head_node];


    move_result.choose_tasks(disturbance_seq.front(), u);

    move_result.move_arguments = disturbance_seq;
    move_result.move_arguments.push_back(u);


    // Check if we need to reduce the # of routes here
    const int start_i = ns.routes[i_route]->start;
    const int end_i = ns.routes[i_route]->end;
    if (start_i == disturbance_seq.front() && end_i == disturbance_seq.back())
        move_result.total_number_of_routes = ns.routes.activated_route_id.size() - 1;
    else
        move_result.total_number_of_routes = ns.routes.activated_route_id.size();


    if (i_route == u_route) {

        move_result.num_affected_routes = 1;

        move_result.route_id.push_back(i_route);

        move_result.delta = u_route_length_delta + i_route_length_delta;
        move_result.route_lens.push_back(ns.routes[i_route]->length + move_result.delta);

        move_result.route_loads.push_back(ns.routes[i_route]->load);

        move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.vio_load_delta = vio_load_delta;

        move_result.considerable = true;

        move_result.route_time_tbl.emplace_back(new_time_tbl[0]);
        move_result.vio_time_delta =
            mcgrp.get_vio_time(move_result.route_time_tbl[0])
                - mcgrp.get_vio_time(ns.routes[i_route]->time_table);
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

        move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers - disturbance_seq.size());
        move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers + disturbance_seq.size());

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.vio_load_delta = vio_load_delta;

        move_result.considerable = true;

        move_result.route_time_tbl = new_time_tbl;
        move_result.vio_time_delta =
            mcgrp.get_vio_time(move_result.route_time_tbl[0])
                + mcgrp.get_vio_time(move_result.route_time_tbl[1])
                - mcgrp.get_vio_time(ns.routes[i_route]->time_table)
                - mcgrp.get_vio_time(ns.routes[u_route]->time_table);
        return true;
    }

    My_Assert(false,"Cannot reach here!");
}


void XPreInsert::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{

#ifdef DEBUG
    hit_count++;
#endif

    DEBUG_PRINT("execute a " + to_string(length) + "-insert:presert move");
    My_Assert(move_result.considerable,"Invalid predictions");

    //extract move arguments
    const int u = move_result.move_arguments.back();
    vector<int> disturbance_sequence(move_result.move_arguments.begin(),move_result.move_arguments.end() - 1);

    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(all_of(disturbance_sequence.begin(),disturbance_sequence.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");

    int edge_num = mcgrp.count_edges(disturbance_sequence);
    if(edge_num > 0 && move_result.num_affected_routes == 2){
        const int i_route = move_result.route_id[0];
        const int u_route = move_result.route_id[1];

        ns.routes[i_route]->num_edges -= edge_num;
        ns.routes[u_route]->num_edges += edge_num;
    }

    if (move_result.total_number_of_routes != ns.routes.activated_route_id.size()) {
        //need to eliminate an empty route
        My_Assert(move_result.total_number_of_routes == ns.routes.activated_route_id.size() - 1,"Incorrect routes number");
        My_Assert(move_result.num_affected_routes == 2,"Incorrect routes number");
        const int i_route = move_result.route_id[0];
        const int u_route = move_result.route_id[1];
        My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");
        My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");
        My_Assert(ns.routes[i_route]->start == disturbance_sequence.front() && ns.routes[i_route]->end == disturbance_sequence.back(),"Wrong route");

        //Modify routes info
        ns.routes[u_route]->length = move_result.route_lens[1];

        ns.routes[u_route]->load = move_result.route_loads[1];

        ns.routes[u_route]->num_customers = move_result.route_custs_num[1];

        ns.routes[u_route]->time_table = move_result.route_time_tbl[1];

        for(auto task:disturbance_sequence){
            ns.solution[task]->route_id = u_route;
        }

        if(ns.solution[u]->pre->ID < 0){
            ns.routes[u_route]->start = disturbance_sequence.front();
        }

        ns.routes.free_route(i_route);

        //handle solution
        //record next dummy Task of actual u
        int dummy_marker = 0;
        My_Assert(ns.solution[disturbance_sequence.front()]->pre->ID < 0 && ns.solution[disturbance_sequence.back()]->next->ID < 0,"Route is not empty!");
        if(ns.solution[disturbance_sequence.back()]->next == ns.solution.very_end){
            dummy_marker = ns.solution[disturbance_sequence.front()]->pre->ID;
        }
        else{
            dummy_marker = ns.solution[disturbance_sequence.back()]->next->ID;
        }

        //handle solution
        //remove
        ns.solution[disturbance_sequence.front()]->pre->next = ns.solution[disturbance_sequence.back()]->next;
        ns.solution[disturbance_sequence.back()]->next->pre = ns.solution[disturbance_sequence.front()]->pre;

        //presert
        ns.solution[disturbance_sequence.front()]->pre = ns.solution[u]->pre;
        ns.solution[disturbance_sequence.back()]->next = ns.solution[u];
        ns.solution[u]->pre->next = ns.solution[disturbance_sequence.front()];
        ns.solution[u]->pre = ns.solution[disturbance_sequence.back()];

        //free dummy marker
        ns.solution[dummy_marker]->pre->next = ns.solution[dummy_marker]->next;
        ns.solution[dummy_marker]->next->pre = ns.solution[dummy_marker]->pre;
        ns.solution.dummypool.free_dummy(dummy_marker);
    }
    else{
        //no routes eliminates
        //Modify routes info
        if (move_result.num_affected_routes == 1) {
            //No need to change route id
            const int route_id = move_result.route_id[0];
            My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[route_id]->length = move_result.route_lens[0];

            ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

            if(ns.solution[disturbance_sequence.front()]->pre->ID < 0){
                ns.routes[route_id]->start = ns.solution[disturbance_sequence.back()]->next->ID;
            }

            if(ns.solution[disturbance_sequence.back()]->next->ID < 0){
                ns.routes[route_id]->end = ns.solution[disturbance_sequence.front()]->pre->ID;
            }

            if(ns.solution[u]->pre->ID < 0){
                ns.routes[route_id]->start = disturbance_sequence.front();
            }
        }
        else{
            const int i_route = move_result.route_id[0];
            const int u_route = move_result.route_id[1];

            My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[i_route]->length = move_result.route_lens[0];
            ns.routes[u_route]->length = move_result.route_lens[1];

            ns.routes[i_route]->load = move_result.route_loads[0];
            ns.routes[u_route]->load = move_result.route_loads[1];

            ns.routes[i_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[u_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[i_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[u_route]->time_table = move_result.route_time_tbl[1];

            for(auto task:disturbance_sequence){
                ns.solution[task]->route_id = u_route;
            }

            if(ns.solution[disturbance_sequence.front()]->pre->ID < 0){
                ns.routes[i_route]->start = ns.solution[disturbance_sequence.back()]->next->ID;
            }

            if(ns.solution[disturbance_sequence.back()]->next->ID < 0){
                ns.routes[i_route]->end = ns.solution[disturbance_sequence.front()]->pre->ID;
            }

            if(ns.solution[u]->pre->ID < 0){
                ns.routes[u_route]->start = disturbance_sequence.front();
            }
        }

        //handle solution
        //remove
        ns.solution[disturbance_sequence.front()]->pre->next = ns.solution[disturbance_sequence.back()]->next;
        ns.solution[disturbance_sequence.back()]->next->pre = ns.solution[disturbance_sequence.front()]->pre;

        //presert
        ns.solution[disturbance_sequence.front()]->pre = ns.solution[u]->pre;
        ns.solution[disturbance_sequence.back()]->next = ns.solution[u];
        ns.solution[u]->pre->next = ns.solution[disturbance_sequence.front()];
        ns.solution[u]->pre = ns.solution[disturbance_sequence.back()];

    }


    //modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");


    if(move_result.delta == 0){
        ns.equal_step++;
    }


    update_score(ns);

    ns.trace(mcgrp);
    move_result.reset();
    ns.search_step++;
}


bool XPreInsert::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

#ifdef DEBUG
    attempt_count = 0;
    hit_count = 0;
#endif

    MoveResult BestM;

    vector<int> chosen_seq = ns.get_successor_tasks(length, chosen_task);

    int offset = 0;
    if(ns.policy.has_rule(TOLERANCE)){
        offset = rand() % mcgrp.actual_task_num;
    }
    ns.create_search_neighborhood(mcgrp, chosen_seq,"successor", offset);


    My_Assert(chosen_seq.size()>0,"You cannot generate an empty sequence!");

    if (chosen_seq.size() != length) {
        DEBUG_PRINT("The length after chosen Task is not " + to_string(length) + " in " + to_string(length) + "-insert");
        return false;
    }

    int b = chosen_seq.front();

    for(auto neighbor_task : ns.search_space) {
        My_Assert(neighbor_task != b, "neighbor Task can't be itself!");

        if (neighbor_task != DUMMY) {

#ifdef DEBUG
            attempt_count++;
#endif

            //j can't be dummy and b can't be dummy neither here
            int j = neighbor_task;

            if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                // doesn't overlap
                if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(move_result)) {
                    if (ns.policy.has_rule(FIRST_ACCEPT)) {
                        move(ns, mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)) {
                        if (ns.policy.check_result(move_result, BestM))
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

#ifdef DEBUG
                attempt_count += 2;
#endif

                // Consider the start location
                int j = current_start;

                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.routes[current_route]->end;
                j = current_end;
                if (std::find(chosen_seq.begin(), chosen_seq.end(), j) == chosen_seq.end()) {
                    // doesn't overlap
                    if (considerable_move(ns, mcgrp, chosen_seq, j) && ns.policy.check_move(move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Advance to next route's starting Task
                current_start = ns.solution[current_end]->next->next->ID;
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

bool XPreInsert::update_score(HighSpeedNeighBorSearch &ns)
{
    if(move_result.delta > 0) return false;

    double penalty = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    const int u = move_result.move_arguments.back();
    const int seq_front = move_result.move_arguments.front();
    const int seq_back = *(move_result.move_arguments.end() - 2);

    // ...seq()-u...
    ns.score_matrix[seq_back][u] += penalty;
    ns.score_matrix[max(0,ns.solution[seq_front]->pre->ID)][seq_front] += penalty;

    return true;
}