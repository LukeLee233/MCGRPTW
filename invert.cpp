#include "invert.h"
#include <algorithm>

using namespace std;


bool Invert::search(HighSpeedNeighBorSearch &ns, const class MCGRP &mcgrp, int chosen_task)
{
    //No search space in Invert operator, No accept rule for invert operator
    My_Assert(chosen_task != DUMMY, "Chosen Task can't be dummy");

#ifdef DEBUG
    attempt_count++;
#endif

    if (!mcgrp.is_edge(chosen_task)) {

#ifdef DEBUG
        attempt_count--;
#endif

        return false;
    }
    else if (considerable_move(ns, mcgrp, chosen_task) && ns.policy.check_move(mcgrp,ns, move_result)) {
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
    My_Assert(u != DUMMY, "Task u cannot be dummy Task!");
    My_Assert(mcgrp.is_edge(u), "Task u must be edge Task!");

    bool allow_infeasible = ns.policy.has_rule(INFEASIBLE) ? true : false;

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
    vector<RouteInfo::TimeTable> new_time_tbl{{{-1, -1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,u,u_tilde,allow_infeasible);
    if(!mcgrp.isTimeTableFeasible(new_time_tbl)){
        move_result.reset();
        return false;
    }

    if(allow_infeasible){
        vector<vector<RouteInfo::TimeTable>> pack_new_time_tbl;
        pack_new_time_tbl.push_back(new_time_tbl);
        if(ns.policy.check_time_window(mcgrp,pack_new_time_tbl)){
            move_result.reset();
            return false;
        }
    }

    move_result.choose_tasks(u, u_tilde);
    move_result.delta = delta;
    move_result.route_id.push_back(u_route);
    move_result.route_lens.push_back(ns.routes[u_route]->length + move_result.delta);

    move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;
    move_result.move_arguments_bak["input"] = {u};
    move_result.move_arguments_bak["output"] = {u_tilde};
    move_result.considerable = true;

    move_result.route_time_tbl.emplace_back(new_time_tbl);


    if(ns.policy.has_rule(INFEASIBLE)){
        auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
        auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

        move_result.vio_load_delta = 0;
        move_result.vio_time_delta = new_time_info_u.second - old_time_info_u.second;
        move_result.vio_time_custom_num_delta = new_time_info_u.first - old_time_info_u.first;
    }else{
        move_result.vio_load_delta = 0;
        move_result.vio_time_delta = 0;
        move_result.vio_time_custom_num_delta = 0;
    }
    return true;
}

void Invert::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{

#ifdef DEBUG
    hit_count++;
#endif

    DEBUG_PRINT("execute an invert move");

    My_Assert(move_result.considerable,"Invalid predictions");

    const int u = move_result.move_arguments_bak.at("input")[0];
    const int u_tilde = move_result.move_arguments_bak.at("output")[0];

    My_Assert(u != DUMMY, "Invert can't handle dummy Task!");
    My_Assert(u_tilde == mcgrp.inst_tasks[u].inverse, "Invert can't handle dummy Task!");
    My_Assert(ns.solution[u]->next != nullptr && ns.solution[u]->pre != nullptr,"Wrong arguments");
    My_Assert(ns.solution[u_tilde]->next == nullptr && ns.solution[u_tilde]->pre == nullptr,"Wrong arguments");



    const int route_id = move_result.route_id[0];
    My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");

    ns.routes[route_id]->length = move_result.route_lens[0];

    ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

    ns.solution[u_tilde]->route_id = route_id;

    //handle solution
    //connect
    ns.solution[u]->pre->next = ns.solution[u_tilde];
    ns.solution[u]->next->pre = ns.solution[u_tilde];
    ns.solution[u_tilde]->next = ns.solution[u]->next;
    ns.solution[u_tilde]->pre = ns.solution[u]->pre;


    //clear original u Task info
    ns.solution[u]->clear();


    //modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");


    update_stable_likelihood(mcgrp,ns,{u_tilde},move_result);

    update_score(ns);

    ns.trace(mcgrp);

    move_result.reset();
}

vector<RouteInfo::TimeTable>
Invert::expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int u, int u_tilde, bool allow_infeasible)
{
    vector<RouteInfo::TimeTable> res({{-1, -1}});

    const int u_route = ns.solution[u]->route_id;

    vector<int> route_task;
    for(const auto& node:ns.routes[u_route]->time_table){
        if(node.task == u)
            route_task.push_back(u_tilde);
        else
            route_task.push_back(node.task);
    }

    vector<int> time_tbl = mcgrp.cal_arrive_time(route_task);
    vector<RouteInfo::TimeTable> intermediate;

    My_Assert(time_tbl.size() == route_task.size(), "wrong time table");
    for(int k = 0; k<time_tbl.size(); k++){
        if(!allow_infeasible && time_tbl[k] > mcgrp.inst_tasks[route_task[k]].time_window.second)
            return res;
        intermediate.push_back({route_task[k],time_tbl[k]});
    }

    res = intermediate;
    return res;
}

bool Invert::update_score(HighSpeedNeighBorSearch &ns)
{
    if(move_result.delta > 0) return false;

    const int actual_u = move_result.move_arguments_bak.at("output")[0];
    double penalty = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    ns.score_matrix[actual_u][max(0,ns.solution[actual_u]->next->ID)] += penalty;
    ns.score_matrix[max(0,ns.solution[actual_u]->pre->ID)][actual_u] += penalty;

    return true;
}
