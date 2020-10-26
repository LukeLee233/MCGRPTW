#include "Invert.h"
#include <algorithm>

using namespace std;



bool Invert::search(HighSpeedNeighBorSearch &ns, const class MCGRP &mcgrp, int chosen_task)
{
    //No search space in Invert operator, No accept rule for invert operator
    My_Assert(chosen_task != DUMMY, "Chosen task can't be dummy");

    if (!mcgrp.is_edge(chosen_task)) {
        return false;
    }
    else if (considerable_move(ns, mcgrp, chosen_task) && ns.policy.check_move(move_result)) {
        move(ns, mcgrp);
        return true;
    }
    else
    {
        move_result.reset();
        return false;
    }
}

bool Invert::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int u)
{
    My_Assert(u != DUMMY, "task u cannot be dummy task!");
    My_Assert(mcgrp.is_edge(u), "task u must be edge task!");

    bool allow_infeasible = ns.policy.has_rule(FITNESS_ONLY) ? true : false;

    const int u_tilde = mcgrp.inst_tasks[u].inverse;

    int t, v;
    t = max(ns.solution[u]->pre->ID, DUMMY);
    v = max(ns.solution[u]->next->ID, DUMMY);

    const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
    const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];

    const double tu_tilde = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u_tilde].head_node];
    const double u_tildev = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[v].head_node];

    const double
        delta = -tu - uv - mcgrp.inst_tasks[u].serv_cost + tu_tilde + u_tildev + mcgrp.inst_tasks[u_tilde].serv_cost;

    int u_route = ns.solution[u]->route_id;
    vector<MCGRPRoute::Timetable> new_time_tbl{{{-1,-1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,u,u_tilde,allow_infeasible);
    if(!mcgrp.isTimetableFeasible(new_time_tbl)){
        move_result.reset();
        return false;
    }


    move_result.choose_tasks(u, u_tilde);
//    move_result.num_affected_routes = 1;
    move_result.delta = delta;
    move_result.route_id.push_back(u_route);
    move_result.route_lens.push_back(ns.routes[u_route]->length + move_result.delta);

//    move_result.route_loads.push_back(ns.routes[u_route]->load);
//    move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers);
    move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;
//    move_result.total_number_of_routes = ns.routes.activated_route_id.size();
    move_result.move_arguments.push_back(u);
    move_result.move_arguments.push_back(u_tilde);
    move_result.considerable = true;

    move_result.route_time_tbl.emplace_back(new_time_tbl);
    move_result.vio_time_delta =
        mcgrp.get_vio_time(move_result.route_time_tbl[0])
        - mcgrp.get_vio_time(ns.routes[u_route]->time_table);
    return true;
}

void Invert::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    DEBUG_PRINT("execute an invert move");

    My_Assert(move_result.considerable,"Invalid predictions");
    My_Assert(move_result.move_arguments.size() == 2, "Incorrect move arguments!");

    const int u = move_result.move_arguments[0];
    const int u_tilde = move_result.move_arguments[1];

    My_Assert(u != DUMMY, "Invert can't handle dummy task!");
    My_Assert(u_tilde == mcgrp.inst_tasks[u].inverse, "Invert can't handle dummy task!");
    My_Assert(ns.solution[u]->next != nullptr && ns.solution[u]->pre != nullptr,"Wrong arguments");
    My_Assert(ns.solution[u_tilde]->next == nullptr && ns.solution[u_tilde]->pre == nullptr,"Wrong arguments");



    const int route_id = move_result.route_id[0];
    My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");

    ns.routes[route_id]->length = move_result.route_lens[0];

    ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

    if(ns.solution[u]->pre->ID < 0){
        ns.routes[route_id]->start = u_tilde;
    }

    if(ns.solution[u]->next->ID < 0){
        ns.routes[route_id]->end = u_tilde;
    }

    ns.solution[u_tilde]->route_id = route_id;

    //handle solution
    //connect
    ns.solution[u]->pre->next = ns.solution[u_tilde];
    ns.solution[u]->next->pre = ns.solution[u_tilde];
    ns.solution[u_tilde]->next = ns.solution[u]->next;
    ns.solution[u_tilde]->pre = ns.solution[u]->pre;


    //clear original u task info
    ns.solution[u]->clear();


    //modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");

    if(move_result.delta == 0){
        ns.equal_step++;
    }

    ns.trace(mcgrp);

    move_result.reset();
    ns.search_step++;
}

void Invert::unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    vector<int> task_set(mcgrp.actual_task_num);
    std::generate(task_set.begin(), task_set.end(), Generator());
    mcgrp._rng.RandPerm(task_set);    //shuffle tasks


    auto original_policy = ns.policy.get();
    ns.policy.set(FIRST_ACCEPT | DOWNHILL | DELTA_ONLY);
    ns.policy.beta = 0.5;
    ns.policy.tolerance = 0.003;
    ns.neigh_size = mcgrp.neigh_size;

    int chosen_task = -1;
    for (int i = 0; i < mcgrp.actual_task_num; i++) {
        chosen_task = task_set[i];

        if (ns.solution[chosen_task]->next == nullptr) {
            if (mcgrp.is_edge(chosen_task)) {
                chosen_task = mcgrp.inst_tasks[chosen_task].inverse;
                My_Assert(ns.solution[chosen_task]->next != nullptr, "An edge task has been missed");
            }
            else {
                My_Assert(false,"A non edge task has been missed!");
            }
        }

        search(ns, mcgrp, chosen_task);
    }

    ns.neigh_size = 0;
    ns.policy.set(original_policy);
    ns.policy.beta = 0;
    ns.policy.tolerance = 0;
}

vector<MCGRPRoute::Timetable>
Invert::expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int u, int u_tilde, bool allow_infeasible)
{
    vector<MCGRPRoute::Timetable> res({{-1,-1}});

    const int u_route = ns.solution[u]->route_id;

    vector<int> route_task;
    for(const auto& node:ns.routes[u_route]->time_table){
        if(node.task == u)
            route_task.push_back(u_tilde);
        else
            route_task.push_back(node.task);
    }

    vector<int> time_tbl = mcgrp.cal_arrive_time(route_task);
    vector<MCGRPRoute::Timetable> intermediate;

    My_Assert(time_tbl.size() == route_task.size(), "wrong time table");
    for(int k = 0; k<time_tbl.size(); k++){
        if(!allow_infeasible && time_tbl[k] > mcgrp.inst_tasks[route_task[k]].time_window.second)
            return res;
        intermediate.push_back({route_task[k],time_tbl[k]});
    }

    res = intermediate;
    return res;
}
