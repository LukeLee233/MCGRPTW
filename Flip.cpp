#include "Flip.h"
#include <algorithm>

using namespace std;


/*
 *
 */

bool NewFlip::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int start_task, int end_task){
    vector<int> candidate_seq = get_sequence(ns,start_task, end_task);

    //No need flipping
    if(candidate_seq.size() <= 1){
        move_result.reset();
        return false;
    }

    My_Assert(all_of(candidate_seq.begin(),candidate_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong Task");
    My_Assert(ns.solution[candidate_seq.front()]->route_id == ns.solution[candidate_seq.back()]->route_id,"Flip attempted using different routes!");

    vector<RouteInfo::TimeTable> new_time_tbl{{{-1, -1}}};
    bool allow_infeasible = ns.policy.has_rule(FITNESS_ONLY) ? true : false;

    new_time_tbl = expected_time_table(ns,mcgrp,candidate_seq,allow_infeasible);

    if(!mcgrp.isTimeTableFeasible(new_time_tbl)){
        move_result.reset();
        return false;
    }

    move_result.choose_tasks(start_task, end_task);
    move_result.move_arguments = candidate_seq;

    double delta = 0;
    start_task = max(start_task,0);
    end_task = max(end_task,0);

    delta -= mcgrp.min_cost[mcgrp.inst_tasks[start_task].tail_node][mcgrp.inst_tasks[candidate_seq.front()].head_node];
    for(int cur = 0; cur < candidate_seq.size() - 1; cur++){
        delta -= mcgrp.inst_tasks[candidate_seq[cur]].serv_cost;
        delta -= mcgrp.min_cost[mcgrp.inst_tasks[candidate_seq[cur]].tail_node][mcgrp.inst_tasks[candidate_seq[cur+1]].head_node];
    }
    delta -= mcgrp.inst_tasks[candidate_seq.back()].serv_cost;
    delta -= mcgrp.min_cost[mcgrp.inst_tasks[candidate_seq.back()].tail_node][mcgrp.inst_tasks[end_task].head_node];

    reverse(candidate_seq.begin(),candidate_seq.end());

    delta += mcgrp.min_cost[mcgrp.inst_tasks[start_task].tail_node][mcgrp.inst_tasks[candidate_seq.front()].head_node];

    for(int cur = 0;cur<candidate_seq.size()-1;cur++){
        delta += mcgrp.inst_tasks[candidate_seq[cur]].serv_cost;
        delta += mcgrp.min_cost[mcgrp.inst_tasks[candidate_seq[cur]].tail_node][mcgrp.inst_tasks[candidate_seq[cur + 1]].head_node];
    }

    delta += mcgrp.inst_tasks[candidate_seq.back()].serv_cost;
    delta += mcgrp.min_cost[mcgrp.inst_tasks[candidate_seq.back()].tail_node][mcgrp.inst_tasks[end_task].head_node];


    move_result.num_affected_routes = 1;

    const int route_id = ns.solution[candidate_seq.front()]->route_id;
    move_result.route_id.push_back(route_id);
    move_result.delta = delta;

    move_result.route_lens.push_back(ns.routes[route_id]->length + move_result.delta);
    move_result.route_loads.push_back(ns.routes[route_id]->load);

    move_result.route_custs_num.push_back(ns.routes[route_id]->num_customers);
    move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;
    move_result.considerable = true;

    move_result.route_time_tbl.emplace_back(new_time_tbl);
    move_result.vio_time_delta = mcgrp.get_vio_time(move_result.route_time_tbl[0])
                                - mcgrp.get_vio_time(ns.routes[route_id]->time_table);
    return true;
}

void NewFlip::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp){
    DEBUG_PRINT("execute a 2-opt: flip move");

    My_Assert(move_result.considerable,"Invalid predictions");

    //extract move arguments
    vector<int> seq(move_result.move_arguments.begin(), move_result.move_arguments.end());
    My_Assert(all_of(seq.begin(),seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong Task");
    My_Assert(ns.solution[move_result.task1]->next->ID == seq.front() && ns.solution[move_result.task2]->pre->ID == seq.back(),"Wrong Task");

    My_Assert(move_result.num_affected_routes == 1,"Wrong result info");

    //Modify routes info
    const int route_id = move_result.route_id[0];
    ns.routes[route_id]->length = move_result.route_lens[0];

    ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

    const int start_task = move_result.task1;
    const int end_task = move_result.task2;

    if(ns.solution[start_task]->ID < 0){
        ns.routes[route_id]->start = seq.back();
    }

    if(ns.solution[end_task]->ID < 0){
        ns.routes[route_id]->end = seq.front();
    }

    //handle solution
    for(auto task : seq){
        auto tmp = ns.solution[task]->pre;
        ns.solution[task]->pre = ns.solution[task]->next;
        ns.solution[task]->next = tmp;
    }

    ns.solution[start_task]->next = ns.solution[seq.back()];
    ns.solution[end_task]->pre = ns.solution[seq.front()];
    ns.solution[seq.front()]->next = ns.solution[end_task];
    ns.solution[seq.back()]->pre = ns.solution[start_task];


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

vector<int> NewFlip::get_sequence(HighSpeedNeighBorSearch &ns, const int start, const int end)
{
    vector<int> buffer;
    int cur = ns.solution[start]->next->ID;

    while (cur != end){
        buffer.push_back(cur);
        cur = ns.solution[cur]->next->ID;
    }

    return buffer;
}

vector<RouteInfo::TimeTable>
NewFlip::expected_time_table(HighSpeedNeighBorSearch &ns,
                             const MCGRP &mcgrp,
                             vector<int> &invert_seq,
                             bool allow_infeasible)
{
    vector<RouteInfo::TimeTable>
    res({{-1,-1}});

    const int invert_route = ns.solution[invert_seq.front()]->route_id;

    vector<int> route_task;
    for(int k = 0; k<ns.routes[invert_route]->time_table.size();){
        if(ns.routes[invert_route]->time_table[k].task == invert_seq.front()){
            route_task.insert(route_task.end(),invert_seq.crbegin(),invert_seq.crend());
            k += invert_seq.size();
        }
        else{
            route_task.push_back(ns.routes[invert_route]->time_table[k].task);
            k++;
        }
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
