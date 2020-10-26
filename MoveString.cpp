#include "MoveString.h"
#include "NeighborSearch.h"
#include "PostSert.h"
#include "Presert.h"
#include <algorithm>




bool
PostMoveString::considerable_move(HighSpeedNeighBorSearch &ns,
                                  const MCGRP &mcgrp,
                                  vector<int> disturbance_seq,
                                  const int u)
{
    // task u cannot be dummy task
    My_Assert(u>=1 && u<=mcgrp.actual_task_num,"Wrong task");
    My_Assert(all_of(disturbance_seq.begin(),disturbance_seq.end(),
                     [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong task");

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
                //if insert task to route u and over load
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
                //if remove task from route i and over load
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
    vector<vector<MCGRPRoute::Timetable>> new_time_tbl{{{-1,-1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,disturbance_seq,u,allow_infeasible);

    if(!mcgrp.isTimetableFeasible(new_time_tbl[0])){
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

void PostMoveString::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    DEBUG_PRINT("execute a double insert:postsert move");
    My_Assert(move_result.considerable,"Invalid predictions");

    //extract move arguments
    const int u = move_result.move_arguments.back();
    vector<int> disturbance_sequence(move_result.move_arguments.begin(),move_result.move_arguments.end() - 1);

    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(all_of(disturbance_sequence.begin(),disturbance_sequence.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");

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
        //record next dummy task of actual u
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

    move_result.reset();
    ns.search_step++;
}

vector<vector<MCGRPRoute::Timetable>>
PostMoveString::expected_time_table(HighSpeedNeighBorSearch &ns,
                                    const MCGRP &mcgrp,
                                    vector<int> &disturbance_seq,
                                    const int i,
                                    bool allow_infeasible)
{
    vector<vector<MCGRPRoute::Timetable>>
        res(2,vector<MCGRPRoute::Timetable>({{-1,-1}}));

    const int i_route = ns.solution[i]->route_id;
    const int u_route = ns.solution[disturbance_seq.front()]->route_id;

    auto intermediate = mcgrp.forecast_time_table(ns.routes[u_route]->time_table,
                                                  disturbance_seq,"remove",-1 ,allow_infeasible);

    if(!mcgrp.isTimetableFeasible(intermediate)){
        return res;
    }

    vector<MCGRPRoute::Timetable> final;
    if(u_route == i_route){
        final = mcgrp.forecast_time_table(intermediate,
                                          disturbance_seq,"insert_after",i,allow_infeasible);
        intermediate = final;
    }else{
        final = mcgrp.forecast_time_table(ns.routes[i_route]->time_table,
                                          disturbance_seq,"insert_after",i,allow_infeasible);
    }

    if(!mcgrp.isTimetableFeasible(final)){
        return res;
    }

    // 0 for u_route, 1 for i_route
    res[0] = intermediate;
    res[1] = final;

    return res;
}

//---------------------------------------------------------------------------------------------



bool
PreMoveString::considerable_move(HighSpeedNeighBorSearch &ns,
                                 const class MCGRP &mcgrp,
                                 vector<int> disturbance_seq,
                                 const int u)
{

    // task u cannot be dummy task
    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong task");
    My_Assert(all_of(disturbance_seq.begin(), disturbance_seq.end(),
                     [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}), "Wrong task");

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
                //if insert task to route u and over load
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
                //if remove task from route i and over load
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
    vector<vector<MCGRPRoute::Timetable>> new_time_tbl{{{-1,-1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,disturbance_seq,u,allow_infeasible);

    if(!mcgrp.isTimetableFeasible(new_time_tbl[0])){
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

void PreMoveString::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    DEBUG_PRINT("execute a double insert:postsert move");
    My_Assert(move_result.considerable,"Invalid predictions");

    //extract move arguments
    const int u = move_result.move_arguments.back();
    vector<int> disturbance_sequence(move_result.move_arguments.begin(),move_result.move_arguments.end() - 1);

    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(all_of(disturbance_sequence.begin(),disturbance_sequence.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong arguments");

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
        //record next dummy task of actual u
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

    move_result.reset();
    ns.search_step++;
}

vector<vector<MCGRPRoute::Timetable>> PreMoveString::expected_time_table(HighSpeedNeighBorSearch &ns,
                                                                         const MCGRP &mcgrp,
                                                                         vector<int> &disturbance_seq,
                                                                         const int i,
                                                                         bool allow_infeasible)
{
    vector<vector<MCGRPRoute::Timetable>>
        res(2,vector<MCGRPRoute::Timetable>({{-1,-1}}));

    const int i_route = ns.solution[i]->route_id;
    const int u_route = ns.solution[disturbance_seq.front()]->route_id;

    auto intermediate = mcgrp.forecast_time_table(ns.routes[u_route]->time_table,
                                                  disturbance_seq,"remove",-1 ,allow_infeasible);

    if(!mcgrp.isTimetableFeasible(intermediate)){
        return res;
    }

    vector<MCGRPRoute::Timetable> final;
    if(u_route == i_route){
        final = mcgrp.forecast_time_table(intermediate,
                                          disturbance_seq,"insert_before",i,allow_infeasible);
        intermediate = final;
    }else{
        final = mcgrp.forecast_time_table(ns.routes[i_route]->time_table,
                                          disturbance_seq,"insert_before",i,allow_infeasible);
    }

    if(!mcgrp.isTimetableFeasible(final)){
        return res;
    }

    // 0 for u_route, 1 for i_route
    res[0] = intermediate;
    res[1] = final;

    return res;
}