//
// Created by luke on 2020/1/27.
//


#include "extraction.h"
#include <algorithm>

using namespace std;

bool Extraction::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    // No search space in Extraction operator, No accept rule for invert operator
    return false;
    if(ns.policy.has_rule(INFEASIBLE)) return false;
    if(ns.policy.has_rule(TOLERANCE)) return false;
    if(ns.policy.has_rule(DOWNHILL)) return false;

#ifdef DEBUG
    attempt_count++;
#endif

    My_Assert(chosen_task != DUMMY, "Chosen Task can't be dummy");

    if (considerable_move(ns, mcgrp, chosen_task) && ns.policy.check_move(mcgrp, ns, move_result)) {
        move(ns, mcgrp);
        return true;
    }
    else {
        move_result.reset();
        return false;
    }
}

struct seg_info
{
    int start = 0;
    int end = 0;
    double length = 0;
    double load = 0;
    int num_custs = 0;
};

seg_info get_seg_info(const MCGRP &mcgrp, HighSpeedNeighBorSearch &ns, const int start_task, const int end_task)
{
    My_Assert(start_task >= 1 && start_task <= mcgrp.actual_task_num, "Wrong Task");
    My_Assert(end_task >= 1 && end_task <= mcgrp.actual_task_num, "Wrong Task");

    int current_task = start_task;
    seg_info buffer;

    buffer.start = start_task;
    buffer.end = end_task;

    while (current_task != end_task) {
        int next_task = ns.solution[current_task]->next->ID;
        buffer.length += (mcgrp.inst_tasks[current_task].serv_cost +
            mcgrp.min_cost[mcgrp.inst_tasks[current_task].tail_node][mcgrp.inst_tasks[next_task].head_node]);

        buffer.load += mcgrp.inst_tasks[current_task].demand;

        buffer.num_custs += 1;

        current_task = next_task;
    }

    buffer.length += mcgrp.inst_tasks[current_task].serv_cost;
    buffer.load += mcgrp.inst_tasks[current_task].demand;
    buffer.num_custs += 1;

    return buffer;
}

bool Extraction::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b)
{
    My_Assert(b >= 1 && b <= mcgrp.actual_task_num, "Wrong Task");

    move_result.reset();

    const int a = ns.solution[b]->pre->ID;
    const int c = ns.solution[b]->next->ID;

    if (a < 0 || c < 0) {
        // No need do extraction
        // dummy-b-dummy
        // dummy-b-c...
        // ...a-b-dummy
        move_result.move_type = NeighborOperator::EXTRACTION;
        move_result.reset();
        return false;
    }
    else {
        //...a-b-c...
        My_Assert(a >= 1 && a <= mcgrp.actual_task_num, "Wrong Task");
        My_Assert(c >= 1 && c <= mcgrp.actual_task_num, "Wrong Task");

        int old_route_id = ns.solution[b]->route_id;
        const auto old_route_start = ns.solution[ns.routes[old_route_id]->start]->next->ID;
        const auto old_route_end = ns.solution[ns.routes[old_route_id]->end]->pre->ID;

        auto a_seg = get_seg_info(mcgrp, ns, old_route_start, a);
        auto c_seg = get_seg_info(mcgrp, ns, c, old_route_end);

        My_Assert(a_seg.num_custs > 0 && c_seg.num_custs > 0, "Wrong tasks");


        const auto dummy_old_route_start =
            mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[old_route_start].head_node];
        const auto a_dummy = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
        const auto dummy_b = mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[b].head_node];
        const auto b_dummy = mcgrp.min_cost[mcgrp.inst_tasks[b].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
        const auto dummy_c = mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[c].head_node];
        const auto old_route_end_dummy =
            mcgrp.min_cost[mcgrp.inst_tasks[old_route_end].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

        const auto ab = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[b].head_node];
        const auto bc = mcgrp.min_cost[mcgrp.inst_tasks[b].tail_node][mcgrp.inst_tasks[c].head_node];

        My_Assert(dummy_old_route_start + a_seg.length + ab + mcgrp.inst_tasks[b].serv_cost + bc + c_seg.length
                      + old_route_end_dummy == ns.routes[old_route_id]->length, "Wrong tasks");

        move_result.task1 = b;
        move_result.move_arguments_bak["cut_position"] = vector<int>();
        move_result.move_arguments_bak["cut_position"].push_back(old_route_start);
        move_result.move_arguments_bak["cut_position"].push_back(a);
        move_result.move_arguments_bak["cut_position"].push_back(b);
        move_result.move_arguments_bak["cut_position"].push_back(c);
        move_result.move_arguments_bak["cut_position"].push_back(old_route_end);

        move_result.route_id.push_back(old_route_id);

        move_result.route_loads.push_back(a_seg.load);
        move_result.route_loads.push_back(mcgrp.inst_tasks[b].demand);
        move_result.route_loads.push_back(c_seg.load);


        const double a_route_length = dummy_old_route_start + a_seg.length + a_dummy;
        const double b_route_length = dummy_b + mcgrp.inst_tasks[b].serv_cost + b_dummy;
        const double c_route_length = dummy_c + c_seg.length + old_route_end_dummy;
        move_result.route_lens.push_back(a_route_length);
        move_result.route_lens.push_back(b_route_length);
        move_result.route_lens.push_back(c_route_length);

        const auto delta = -ab - bc + a_dummy + dummy_b + b_dummy + dummy_c;
        move_result.delta = delta;

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.route_custs_num.push_back(a_seg.num_custs);
        move_result.route_custs_num.push_back(1);
        move_result.route_custs_num.push_back(c_seg.num_custs);

        move_result.vio_load_delta = 0;

        move_result.considerable = true;

        move_result.vio_time_delta = 0;
        move_result.route_time_tbl = expected_time_table(ns, mcgrp, b);
        return true;
    }

    My_Assert(false, "Cannot reach here!");
}

void Extraction::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{

#ifdef DEBUG
    hit_count++;
#endif

    DEBUG_PRINT("execute a extraction move");

    My_Assert(move_result.considerable, "Invalid predictions");

    //phase 0: extract move arguments
    const int a_route_start = move_result.move_arguments_bak.at("cut_position")[0];
    const int a_route_end = move_result.move_arguments_bak.at("cut_position")[1];
    const int b = move_result.move_arguments_bak.at("cut_position")[2];
    const int c_route_start = move_result.move_arguments_bak.at("cut_position")[3];
    const int c_route_end = move_result.move_arguments_bak.at("cut_position")[4];

    My_Assert(a_route_start >= 1 && a_route_start <= mcgrp.actual_task_num, "Wrong arguments");
    My_Assert(a_route_end >= 1 && a_route_end <= mcgrp.actual_task_num, "Wrong arguments");
    My_Assert(b >= 1 && b <= mcgrp.actual_task_num, "Wrong arguments");
    My_Assert(c_route_start >= 1 && c_route_start <= mcgrp.actual_task_num, "Wrong arguments");
    My_Assert(c_route_end >= 1 && c_route_end <= mcgrp.actual_task_num, "Wrong arguments");

    // phase 1: update the route level info
    const int a_route = move_result.route_id[0];
    const auto b_route = ns.routes.allocate_route();
    const auto c_route = ns.routes.allocate_route();

    My_Assert(ns.routes.activated_route_id.find(a_route) != ns.routes.activated_route_id.end(), "Invalid route");
    My_Assert(ns.routes.activated_route_id.find(b_route) != ns.routes.activated_route_id.end(), "Invalid route");
    My_Assert(ns.routes.activated_route_id.find(c_route) != ns.routes.activated_route_id.end(), "Invalid route");

    ns.routes[a_route]->length = move_result.route_lens[0];
    ns.routes[b_route]->length = move_result.route_lens[1];
    ns.routes[c_route]->length = move_result.route_lens[2];

    ns.routes[a_route]->num_customers = move_result.route_custs_num[0];
    ns.routes[b_route]->num_customers = move_result.route_custs_num[1];
    ns.routes[c_route]->num_customers = move_result.route_custs_num[2];

    ns.routes[a_route]->load = move_result.route_loads[0];
    ns.routes[b_route]->load = move_result.route_loads[1];
    ns.routes[c_route]->load = move_result.route_loads[2];

    ns.routes[a_route]->time_table = move_result.route_time_tbl[0];
    ns.routes[b_route]->time_table = move_result.route_time_tbl[1];
    ns.routes[c_route]->time_table = move_result.route_time_tbl[2];

    ns.routes[b_route]->num_edges = 0;
    if(mcgrp.is_edge(ns.routes[b_route]->time_table[0].task)){
        ns.routes[b_route]->num_edges = 1;
    }

    ns.routes[c_route]->num_edges = 0;
    for(const auto& item : ns.routes[c_route]->time_table){
        if(mcgrp.is_edge(item.task)) ns.routes[c_route]->num_edges++;
    }

    ns.routes[a_route]->num_edges -= (ns.routes[b_route]->num_edges + ns.routes[c_route]->num_edges);


    ns.solution[b]->route_id = b_route;

    for (auto cur = c_route_start; cur != ns.solution[c_route_end]->next->ID; cur = ns.solution[cur]->next->ID) {
        ns.solution[cur]->route_id = c_route;
    }

    // phase 3: update the sequence level info
    auto new_dummy = ns.solution.dummypool.get_new_dummy();
    new_dummy->pre = ns.solution[a_route_end];
    new_dummy->next = ns.solution[b];
    ns.solution[a_route_end]->next = new_dummy;
    ns.solution[b]->pre = new_dummy;

    new_dummy = ns.solution.dummypool.get_new_dummy();
    new_dummy->pre = ns.solution[b];
    new_dummy->next = ns.solution[c_route_start];
    ns.solution[c_route_start]->pre = new_dummy;
    ns.solution[b]->next = new_dummy;

    // phase 4: set the dummy tasks of routes
    ns.routes[c_route]->end = ns.routes[a_route]->end;
    ns.routes[a_route]->end = ns.solution[a_route_end]->next->ID;

    ns.routes[b_route]->start = ns.solution[b]->pre->ID;
    ns.routes[b_route]->end = ns.solution[b]->next->ID;

    ns.routes[c_route]->start = ns.solution[c_route_start]->pre->ID;


    // phase 5: modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp), "Prediction wrong!");

    update_score(ns);
    ns.trace(mcgrp);

    move_result.reset();
}

vector<vector<RouteInfo::TimeTable>>
Extraction::expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b)
{
    vector<vector<RouteInfo::TimeTable>> res;

    const int b_route = ns.solution[b]->route_id;

    vector<int> original_route;

    for (const auto &task : ns.routes[b_route]->time_table) {
        original_route.push_back(task.task);
    }


    auto b_ite = find(original_route.begin(), original_route.end(), b);
    My_Assert(b_ite != original_route.end(), "Cannot find Task b!");

    vector<int> left_route(original_route.begin(), b_ite);
    vector<int> middle_route{b};
    vector<int> right_route(b_ite + 1, original_route.end());

    vector<int> left_time_tbl = mcgrp.cal_arrive_time(left_route);
    vector<int> middle_time_tbl = mcgrp.cal_arrive_time(middle_route);
    vector<int> right_time_tbl = mcgrp.cal_arrive_time(right_route);

    vector<RouteInfo::TimeTable> intermediate_left;
    vector<RouteInfo::TimeTable> intermediate_middle;
    vector<RouteInfo::TimeTable> intermediate_right;

    for (int i = 0; i < left_route.size(); i++)
        intermediate_left.push_back({left_route[i], left_time_tbl[i]});

    for (int i = 0; i < middle_route.size(); i++)
        intermediate_middle.push_back({middle_route[i], middle_time_tbl[i]});

    for (int i = 0; i < right_route.size(); i++)
        intermediate_right.push_back({right_route[i], right_time_tbl[i]});

    res.emplace_back(intermediate_left);
    res.emplace_back(intermediate_middle);
    res.emplace_back(intermediate_right);

    return res;
}

bool Extraction::update_score(HighSpeedNeighBorSearch &ns)
{
    My_Assert(move_result.delta >= 0, "error, wrong outcome!");
    return false;
}
