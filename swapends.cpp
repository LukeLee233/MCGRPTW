#include "swapends.h"
#include <algorithm>

using namespace std;


/*
 * High Speed
 */

RouteSegment get_segment_info(const MCGRP &mcgrp, HighSpeedNeighBorSearch &ns, const int chosen_task)
{
    // Calculate the length, load, and # of customers of the segment
    // after chosen Task.
    // Example:  chosen_task-a-i-j-b-dummy has
    // length: d(a,i) + d(i,j) + d(j,b)
    // load:   a + i + j + b
    // #:      4


    RouteSegment buffer;

    buffer.segment_start = ns.solution[chosen_task]->next->ID;

    if (buffer.segment_start < 0) {
        buffer.segment_start = 0;
        buffer.segment_end = 0;
        buffer.num_custs = 0;
        buffer.load = 0;
        buffer.len = 0;
        return buffer;
    }

    const int route = ns.solution[buffer.segment_start]->route_id;
    buffer.segment_end = ns.solution[ns.routes[route]->end]->pre->ID;
    buffer.num_custs = 0;
    buffer.load = 0;

    buffer.len = 0;

    // Now calculate the length, load, and # customer on this segment
    int current_task = buffer.segment_start;

    while (current_task != buffer.segment_end) {
        int next_task = ns.solution[current_task]->next->ID;
        buffer.len += (mcgrp.inst_tasks[current_task].serv_cost +
            mcgrp.min_cost[mcgrp.inst_tasks[current_task].tail_node][mcgrp.inst_tasks[next_task].head_node]);

        buffer.load += mcgrp.inst_tasks[current_task].demand;

        buffer.num_custs += 1;

        current_task = next_task;
    }

    //handle last Task in the segment
    buffer.len += mcgrp.inst_tasks[current_task].serv_cost;
    buffer.load += mcgrp.inst_tasks[current_task].demand;
    buffer.num_custs += 1;

    return buffer;
}

bool
NewSwapEnds::considerable_move(HighSpeedNeighBorSearch &ns,
                               const MCGRP &mcgrp,
                               const int chosen_task,
                               const int neighbor_task,
                               const int chosen_route,
                               const int neighbor_route)
{
    // Example: ( a & v input):
    // VRP_DEPOT-i-a-b-j-k-l-VRP_DEPOT and VRP_DEPOT-t-u-v-w-x-y-z-VRP_DEPOT becomes
    // VRP_DEPOT-i-a-w-x-y-z-VRP_DEPOT and VRP_DEPOT-t-u-v-b-j-k-l-VRP_DEPOT

    move_result.reset();

    if (chosen_task < 0 && chosen_task < 0) {
        //means no actual move happens
        //...-(a|0)-[...]-0
        //...-(v|0)-[...]-0
        move_result.reset();
        return false;
    }

    const int a_route = chosen_route;
    const int v_route = neighbor_route;

    My_Assert(a_route != v_route, "swap ends called with a and v in same route!");
    move_result.choose_tasks(chosen_task,neighbor_task);

    //get sequence that need to be swap
    RouteSegment seg_after_a = get_segment_info(mcgrp, ns, chosen_task);
    RouteSegment seg_after_v = get_segment_info(mcgrp, ns, neighbor_task);

    if (seg_after_a.num_custs == 0 && seg_after_v.num_custs == 0) {
        //means no actual move happens
        //...-a-[]-0
        //...-v-[]-0
        move_result.reset();
        return false;
    }

    const auto a_load_delta = -seg_after_a.load + seg_after_v.load;
    const auto v_load_delta = -seg_after_v.load + seg_after_a.load;

    // load legality is independent with mode choice
    double vio_load_delta = 0;
    if (ns.policy.has_rule(FEASIBLE)) {
        if (ns.routes[a_route]->load + a_load_delta > mcgrp.capacity
            || ns.routes[v_route]->load + v_load_delta > mcgrp.capacity) {
            return false;
        }
    }
    else if (ns.policy.has_rule(INFEASIBLE)) {
        My_Assert(v_load_delta == -a_load_delta, "Wrong arguments");

        int pseudo_capacity = ns.policy.get_pseudo_capacity(mcgrp.capacity);
        if (ns.routes[a_route]->load + a_load_delta > pseudo_capacity) {
            move_result.reset();
            return false;
        }

        if (ns.routes[v_route]->load + v_load_delta > pseudo_capacity) {
            move_result.reset();
            return false;
        }

        if (a_load_delta >= 0) {
            //a_route vio-load calculate
            if (ns.routes[a_route]->load + a_load_delta > mcgrp.capacity) {
                //if insert Task to route a and over load
                if (ns.routes[a_route]->load >= mcgrp.capacity) {
                    //if the route a already over loaded
                    vio_load_delta += a_load_delta;
                }
                else {
                    vio_load_delta += ns.routes[a_route]->load + a_load_delta - mcgrp.capacity;
                }
            }

            My_Assert(v_load_delta <= 0, "Wrong arguments!");
            if (ns.routes[v_route]->load > mcgrp.capacity) {
                //if insert Task to route a and over load
                if (ns.routes[v_route]->load + v_load_delta >= mcgrp.capacity) {
                    //if the route v still over loaded
                    vio_load_delta += v_load_delta;
                }
                else {
                    vio_load_delta -= (ns.routes[v_route]->load - mcgrp.capacity);
                }
            }
        }
        else {
            if (ns.routes[a_route]->load > mcgrp.capacity) {
                //if insert Task to route a and over load
                if (ns.routes[a_route]->load + a_load_delta >= mcgrp.capacity) {
                    //if the route v still over loaded
                    vio_load_delta += a_load_delta;
                }
                else {
                    vio_load_delta -= (ns.routes[a_route]->load - mcgrp.capacity);
                }
            }

            My_Assert(v_load_delta >= 0, "Wrong arguments!");
            //v_route vio-load calculate
            if (ns.routes[v_route]->load + v_load_delta > mcgrp.capacity) {
                //if insert Task to route v and over load
                if (ns.routes[v_route]->load >= mcgrp.capacity) {
                    //if the route v already over loaded
                    vio_load_delta += v_load_delta;
                }
                else {
                    vio_load_delta += ns.routes[v_route]->load + v_load_delta - mcgrp.capacity;
                }
            }
        }

    }


    vector<int> input_a_seq;
    int cur_task = ns.solution[chosen_task]->ID;
    for (auto i = 0; i < seg_after_a.num_custs; i++) {
        cur_task = ns.solution[cur_task]->next->ID;
        input_a_seq.push_back(cur_task);
    }

    vector<int> input_v_seq;
    cur_task = ns.solution[neighbor_task]->ID;
    for (auto i = 0; i < seg_after_v.num_custs; i++) {
        cur_task = ns.solution[cur_task]->next->ID;
        input_v_seq.push_back(cur_task);
    }

    if(ns.policy.has_rule(BEST_ACCEPT) && (mcgrp.count_edges(input_a_seq) + mcgrp.count_edges(input_v_seq) > 0)){
//    if(false){
        vector<int> first_seq_a;
        vector<int> first_seq_v;

        vector<int> a_route_seq = ns.get_sub_seq(ns.routes[a_route]->start, ns.routes[a_route]->end);
        vector<int> v_route_seq = ns.get_sub_seq(ns.routes[v_route]->start, ns.routes[v_route]->end);

        if(input_a_seq.empty()){
            first_seq_a.insert(first_seq_a.end(),a_route_seq.begin(),a_route_seq.end() - 1);
        }
        else{
            int a_seq_loc = -1;

            for(int ii = 0;ii < a_route_seq.size(); ii++){
                if(a_route_seq[ii] == input_a_seq.front()) a_seq_loc = ii;
            }

            first_seq_a.insert(first_seq_a.end(),a_route_seq.begin(),a_route_seq.begin() + a_seq_loc);
        }

        if(input_v_seq.empty()){
            first_seq_v.insert(first_seq_v.end(),v_route_seq.begin(),v_route_seq.end() - 1);
        }
        else{
            int v_seq_loc = -1;

            for(int ii = 0;ii < v_route_seq.size(); ii++){
                if(v_route_seq[ii] == input_v_seq.front()) v_seq_loc = ii;
            }

            first_seq_v.insert(first_seq_v.end(), v_route_seq.begin(),v_route_seq.begin() + v_seq_loc);
        }


        My_Assert(first_seq_a.front() < 0 && first_seq_v.front() < 0, "error, wrong sequence!");

        viterbi::PseudoTask first_task_a = Seq2Pseudo(mcgrp,first_seq_a, 1);
        viterbi::PseudoTask first_task_v = Seq2Pseudo(mcgrp,first_seq_v, 1);

        // insert sequence will invoke viterbi result
        auto best_sequence_a = viterbi::viterbi_decode(&first_task_a, nullptr, input_v_seq, mcgrp, ns.policy);
        auto best_sequence_v = viterbi::viterbi_decode(&first_task_v, nullptr, input_a_seq, mcgrp, ns.policy);

        // use best sequence to update the move result
        if(best_sequence_a.cost == INT32_MAX || best_sequence_v.cost == INT32_MAX){
            move_result.reset();
            return false;
        }

        My_Assert(best_sequence_a.seq.front() == INT32_MAX, "error, wrong sequence");
        My_Assert(best_sequence_v.seq.front() == INT32_MAX, "error, wrong sequence");

        double a_route_length_delta = best_sequence_a.cost - ns.routes[a_route]->length;
        double v_route_length_delta = best_sequence_v.cost - ns.routes[v_route]->length;
        move_result.delta = a_route_length_delta + v_route_length_delta;

        move_result.route_lens.push_back(best_sequence_a.cost);
        move_result.route_lens.push_back(best_sequence_v.cost);

        move_result.seq1_cus_num = seg_after_a.num_custs;
        move_result.seq2_cus_num = seg_after_v.num_custs;

        move_result.num_affected_routes = 2;
        My_Assert(move_result.route_id.empty(), "Wrong statement");
        move_result.route_id.push_back(a_route);
        move_result.route_id.push_back(v_route);


        auto new_a_load = ns.routes[a_route]->load + a_load_delta;
        auto new_v_load = ns.routes[v_route]->load + v_load_delta;
        move_result.route_loads.push_back(new_a_load);
        move_result.route_loads.push_back(new_v_load);

        auto new_a_custs_num = ns.routes[a_route]->num_customers - seg_after_a.num_custs + seg_after_v.num_custs;
        auto new_v_custs_num = ns.routes[v_route]->num_customers - seg_after_v.num_custs + seg_after_a.num_custs;
        move_result.route_custs_num.push_back(new_a_custs_num);
        move_result.route_custs_num.push_back(new_v_custs_num);

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.move_arguments_bak["input_a_seq"] = input_a_seq;
        move_result.move_arguments_bak["input_v_seq"] = input_v_seq;
        move_result.move_arguments_bak["output_a_seq"] = vector<int>(best_sequence_v.seq.begin() + 1, best_sequence_v.seq.end());
        move_result.move_arguments_bak["output_v_seq"] = vector<int>(best_sequence_a.seq.begin() + 1, best_sequence_a.seq.end());

        move_result.considerable = true;

        vector<int> new_a_route_seq;
        new_a_route_seq.insert(new_a_route_seq.end(), first_seq_a.begin() + 1, first_seq_a.end());
        new_a_route_seq.insert(new_a_route_seq.end(), move_result.move_arguments_bak["output_v_seq"].begin(), move_result.move_arguments_bak["output_v_seq"].end());
        move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(new_a_route_seq, mcgrp.cal_arrive_time(new_a_route_seq)));

        vector<int> new_v_route_seq;
        new_v_route_seq.insert(new_v_route_seq.end(), first_seq_v.begin() + 1, first_seq_v.end());
        new_v_route_seq.insert(new_v_route_seq.end(), move_result.move_arguments_bak["output_a_seq"].begin(), move_result.move_arguments_bak["output_a_seq"].end());
        move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(new_v_route_seq, mcgrp.cal_arrive_time(new_v_route_seq)));

        if(ns.policy.has_rule(INFEASIBLE)){
            move_result.vio_load_delta = vio_load_delta;

            auto old_time_info_a = mcgrp.get_vio_time(ns.routes[a_route]->time_table);
            auto new_time_info_a = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

            auto old_time_info_v = mcgrp.get_vio_time(ns.routes[v_route]->time_table);
            auto new_time_info_v = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

            move_result.vio_time_delta =
                new_time_info_a.second + new_time_info_v.second - old_time_info_a.second - old_time_info_v.second;
            move_result.vio_time_custom_num_delta =
                new_time_info_a.first + new_time_info_v.first - old_time_info_a.first - old_time_info_v.first;
        }else{
            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = 0;
            move_result.vio_time_custom_num_delta = 0;
        }
    }
    else {

        //a can be dummy and v can be dummy too
        const int a = max(chosen_task, 0);
        const int v = max(neighbor_task, 0);

        bool allow_infeasible = ns.policy.has_rule(INFEASIBLE) ? true : false;
        vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};

        new_time_tbl = expected_time_table(ns,mcgrp,a,v,a_route,v_route,seg_after_a,seg_after_v,allow_infeasible);

        if(!mcgrp.isTimeTableFeasible(new_time_tbl[0])){
            move_result.reset();
            return false;
        }

        if(allow_infeasible){
            if(ns.policy.check_time_window(mcgrp,new_time_tbl)){
                move_result.reset();
                return false;
            }
        }

        double delta;
        if (seg_after_a.num_custs == 0) {
            //a segment is empty
            //0-a-[]-0
            //0-v-[w]-0
            My_Assert(v_route == ns.solution[seg_after_v.segment_start]->route_id, "Wrong routes");
            const int w = max(ns.solution[neighbor_task]->next->ID, 0);

            const auto vw = mcgrp.min_cost[mcgrp.inst_tasks[v].tail_node][mcgrp.inst_tasks[w].head_node];
            const auto aw = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[w].head_node];

            const auto a_dummy = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
            const auto v_dummy = mcgrp.min_cost[mcgrp.inst_tasks[v].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

            delta = -a_dummy + aw - vw + v_dummy;

            const auto seg_after_v_dummy =
                mcgrp.min_cost[mcgrp.inst_tasks[seg_after_v.segment_end].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

            auto new_a_len = ns.routes[a_route]->length - a_dummy + aw + seg_after_v.len + seg_after_v_dummy;
            auto new_v_len = ns.routes[v_route]->length - vw - seg_after_v.len - seg_after_v_dummy + v_dummy;
            move_result.route_lens.push_back(new_a_len);
            move_result.route_lens.push_back(new_v_len);

            move_result.seq1_cus_num = seg_after_a.num_custs;
            move_result.seq2_cus_num = seg_after_v.num_custs;
        }
        else if (seg_after_v.num_custs == 0) {
            //v segment is empty
            //0-a-[b]-0
            //0-v-[]-0
            My_Assert(a_route == ns.solution[seg_after_a.segment_start]->route_id, "Wrong routes");

            const int b = max(ns.solution[chosen_task]->next->ID, 0);

            const auto ab = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[b].head_node];
            const auto vb = mcgrp.min_cost[mcgrp.inst_tasks[v].tail_node][mcgrp.inst_tasks[b].head_node];

            const auto a_dummy = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
            const auto v_dummy = mcgrp.min_cost[mcgrp.inst_tasks[v].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

            delta = -v_dummy + vb - ab + a_dummy;

            const auto seg_after_a_dummy =
                mcgrp.min_cost[mcgrp.inst_tasks[seg_after_a.segment_end].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

            auto new_a_len = ns.routes[a_route]->length - ab - seg_after_a.len - seg_after_a_dummy + a_dummy;
            auto new_v_len = ns.routes[v_route]->length - v_dummy + vb + seg_after_a.len + seg_after_a_dummy;
            move_result.route_lens.push_back(new_a_len);
            move_result.route_lens.push_back(new_v_len);

            move_result.seq1_cus_num = seg_after_a.num_custs;
            move_result.seq2_cus_num = seg_after_v.num_custs;
        }
        else {
            //normal case,two segments are not empty
            //0-a-[b]-0
            //0-v-[w]-0
            My_Assert(a_route == ns.solution[seg_after_a.segment_start]->route_id, "Wrong routes");
            My_Assert(v_route == ns.solution[seg_after_v.segment_start]->route_id, "Wrong routes");

            const int b = max(ns.solution[chosen_task]->next->ID, 0);
            const int w = max(ns.solution[neighbor_task]->next->ID, 0);

            const auto ab = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[b].head_node];
            const auto vw = mcgrp.min_cost[mcgrp.inst_tasks[v].tail_node][mcgrp.inst_tasks[w].head_node];
            const auto aw = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[w].head_node];
            const auto vb = mcgrp.min_cost[mcgrp.inst_tasks[v].tail_node][mcgrp.inst_tasks[b].head_node];

            delta = -ab - vw + aw + vb;

            const auto seg_after_a_dummy =
                mcgrp.min_cost[mcgrp.inst_tasks[seg_after_a.segment_end].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
            const auto seg_after_v_dummy =
                mcgrp.min_cost[mcgrp.inst_tasks[seg_after_v.segment_end].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

            auto new_a_len = ns.routes[a_route]->length - ab - seg_after_a.len - seg_after_a_dummy + aw + seg_after_v.len
                + seg_after_v_dummy;
            auto new_v_len = ns.routes[v_route]->length - vw - seg_after_v.len - seg_after_v_dummy + vb + seg_after_a.len
                + seg_after_a_dummy;
            move_result.route_lens.push_back(new_a_len);
            move_result.route_lens.push_back(new_v_len);

            move_result.seq1_cus_num = seg_after_a.num_custs;
            move_result.seq2_cus_num = seg_after_v.num_custs;
        }

        move_result.delta = delta;

        move_result.num_affected_routes = 2;
        My_Assert(move_result.route_id.empty(), "Wrong statement");
        move_result.route_id.push_back(a_route);
        move_result.route_id.push_back(v_route);


        auto new_a_load = ns.routes[a_route]->load + a_load_delta;
        auto new_v_load = ns.routes[v_route]->load + v_load_delta;
        move_result.route_loads.push_back(new_a_load);
        move_result.route_loads.push_back(new_v_load);

        auto new_a_custs_num = ns.routes[a_route]->num_customers - seg_after_a.num_custs + seg_after_v.num_custs;
        auto new_v_custs_num = ns.routes[v_route]->num_customers - seg_after_v.num_custs + seg_after_a.num_custs;
        move_result.route_custs_num.push_back(new_a_custs_num);
        move_result.route_custs_num.push_back(new_v_custs_num);

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;


        move_result.move_arguments_bak["input_a_seq"] = input_a_seq;
        move_result.move_arguments_bak["input_v_seq"] = input_v_seq;
        move_result.move_arguments_bak["output_a_seq"] = input_a_seq;
        move_result.move_arguments_bak["output_v_seq"] = input_v_seq;

        move_result.considerable = true;

        move_result.route_time_tbl = new_time_tbl;

        if(ns.policy.has_rule(INFEASIBLE)){
            move_result.vio_load_delta = vio_load_delta;

            auto old_time_info_a = mcgrp.get_vio_time(ns.routes[a_route]->time_table);
            auto new_time_info_a = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

            auto old_time_info_v = mcgrp.get_vio_time(ns.routes[v_route]->time_table);
            auto new_time_info_v = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

            move_result.vio_time_delta =
                new_time_info_a.second + new_time_info_v.second - old_time_info_a.second - old_time_info_v.second;
            move_result.vio_time_custom_num_delta =
                new_time_info_a.first + new_time_info_v.first - old_time_info_a.first - old_time_info_v.first;
        }else{
            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = 0;
            move_result.vio_time_custom_num_delta = 0;
        }

        return true;

    }

}

void NewSwapEnds::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    DEBUG_PRINT("execute a 2-opt:swap-ends move");

    My_Assert(move_result.considerable, "Invalid predictions");

    //phase 0: extract move arguments
    const int a = move_result.task1;
    const int v = move_result.task2;

    My_Assert(ns.valid_sol(mcgrp), "Prediction wrong!");


    //Extract the seq which needs swapped
    vector<int> input_a_seq = move_result.move_arguments_bak.at("input_a_seq");
    vector<int> output_a_seq = move_result.move_arguments_bak.at("output_a_seq");
    vector<int> input_v_seq = move_result.move_arguments_bak.at("input_v_seq");
    vector<int> output_v_seq = move_result.move_arguments_bak.at("output_v_seq");


    My_Assert(all_of(input_a_seq.begin(), input_a_seq.end(), [&](int i){ return i >= 1 && i <= mcgrp.actual_task_num; }), "Wrong Task");
    My_Assert(all_of(output_a_seq.begin(), output_a_seq.end(), [&](int i){ return i >= 1 && i <= mcgrp.actual_task_num; }), "Wrong Task");
    My_Assert(all_of(input_v_seq.begin(), input_v_seq.end(), [&](int i){ return i >= 1 && i <= mcgrp.actual_task_num; }), "Wrong Task");
    My_Assert(all_of(output_v_seq.begin(), output_v_seq.end(), [&](int i){ return i >= 1 && i <= mcgrp.actual_task_num; }), "Wrong Task");
    My_Assert(input_a_seq.size() == move_result.seq1_cus_num && input_v_seq.size() == move_result.seq2_cus_num, "Wrong Task");
    My_Assert(move_result.num_affected_routes == 2, "Wrong routes number");

    // phase 1: update the route level info
    const int a_route = move_result.route_id[0];
    const int v_route = move_result.route_id[1];
    My_Assert(a_route != v_route, "Wrong routes!");
    My_Assert(ns.routes.activated_route_id.find(a_route) != ns.routes.activated_route_id.end(), "Invalid route");
    My_Assert(ns.routes.activated_route_id.find(v_route) != ns.routes.activated_route_id.end(), "Invalid route");


    int a_seq_edge_num = mcgrp.count_edges(input_a_seq);
    int v_seq_edge_num = mcgrp.count_edges(input_v_seq);
    ns.routes[a_route]->num_edges += (v_seq_edge_num - a_seq_edge_num);
    ns.routes[v_route]->num_edges += (a_seq_edge_num - v_seq_edge_num);

    ns.routes[a_route]->length = move_result.route_lens[0];
    ns.routes[v_route]->length = move_result.route_lens[1];

    ns.routes[a_route]->load = move_result.route_loads[0];
    ns.routes[v_route]->load = move_result.route_loads[1];

    ns.routes[a_route]->num_customers = move_result.route_custs_num[0];
    ns.routes[v_route]->num_customers = move_result.route_custs_num[1];

    ns.routes[a_route]->time_table = move_result.route_time_tbl[0];
    ns.routes[v_route]->time_table = move_result.route_time_tbl[1];


    // phase 3: update the sequence level info
    const auto a_final_dummy = ns.solution[ns.routes[a_route]->end];
    const auto v_final_dummy = ns.solution[ns.routes[v_route]->end];

    if(input_v_seq == output_v_seq && input_a_seq == output_a_seq){
        if(!input_v_seq.empty()){
            ns.solution[a]->next = ns.solution[input_v_seq.front()];
            ns.solution[input_v_seq.back()]->next = a_final_dummy;
            ns.solution[input_v_seq.front()]->pre = ns.solution[a];
            a_final_dummy->pre = ns.solution[input_v_seq.back()];
        }
        else{
            ns.solution[a]->next = a_final_dummy;
            a_final_dummy->pre = ns.solution[a];
        }

        if(!input_a_seq.empty()){
            ns.solution[v]->next = ns.solution[input_a_seq.front()];
            ns.solution[input_a_seq.back()]->next = v_final_dummy;
            ns.solution[input_a_seq.front()]->pre = ns.solution[v];
            v_final_dummy->pre = ns.solution[input_a_seq.back()];
        }
        else{
            ns.solution[v]->next = v_final_dummy;
            v_final_dummy->pre = ns.solution[v];
        }

    }
    else{
        for(auto input_task : input_a_seq){
            ns.solution[input_task]->clear();
        }

        for(auto input_task : input_v_seq){
            ns.solution[input_task]->clear();
        }

        if(!output_v_seq.empty()){
            ns.solution[a]->next = ns.solution[output_v_seq.front()];
            ns.solution[output_v_seq.back()]->next = a_final_dummy;

            ns.solution.connect_tasks(output_v_seq);

            ns.solution[output_v_seq.front()]->pre = ns.solution[a];
            a_final_dummy->pre = ns.solution[output_v_seq.back()];
        }
        else{
            ns.solution[a]->next = a_final_dummy;
            a_final_dummy->pre = ns.solution[a];
        }

        if(!output_a_seq.empty()){
            ns.solution[v]->next = ns.solution[output_a_seq.front()];
            ns.solution[output_a_seq.back()]->next = v_final_dummy;

            ns.solution.connect_tasks(output_a_seq);

            ns.solution[output_a_seq.front()]->pre = ns.solution[v];
            v_final_dummy->pre = ns.solution[output_a_seq.back()];
        }
        else{
            ns.solution[v]->next = v_final_dummy;
            v_final_dummy->pre = ns.solution[v];
        }
    }

    for(auto task : output_a_seq){
        ns.solution[task]->route_id = v_route;
    }

    for(auto task : output_v_seq){
        ns.solution[task]->route_id = a_route;
    }

    // phase 4: compress the empty routes
    if(ns.routes[a_route]->num_customers == 0){
        ns.compress_empty_route(a_route);
    }

    if(ns.routes[v_route]->num_customers == 0){
        ns.compress_empty_route(v_route);
    }

    // phase 5: modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp), "Prediction wrong!");


    update_score(ns);
    move_result.reset();
}

vector<vector<RouteInfo::TimeTable>>
NewSwapEnds::expected_time_table(HighSpeedNeighBorSearch &ns,
                    const MCGRP &mcgrp,
                    int a,
                    int v,
                    int a_route,
                    int v_route,
                    const RouteSegment& a_seg,
                    const RouteSegment& v_seg,
                    bool allow_infeasible)
{
    vector<vector<RouteInfo::TimeTable>>
        res(2, vector<RouteInfo::TimeTable>({{-1, -1}}));

    My_Assert(a_route != v_route, "Two routes can't be the same route in swap ends!");

    vector<int> route_task_a;
    vector<int> route_task_v;
    int a_pos = -1;
    int v_pos = -1;

    if(a == DUMMY){
        if(a_seg.num_custs == 0)
            a_pos = ns.routes[a_route]->time_table.size();
        else
            a_pos = 0;
    }

    if(v == DUMMY){
        if(v_seg.num_custs == 0)
            v_pos = ns.routes[v_route]->time_table.size();
        else
            v_pos = 0;
    }

    if(a_pos == v_pos && a_pos != -1){
        res[0] = ns.routes[v_route]->time_table;
        res[1] = ns.routes[a_route]->time_table;
        if(a_pos > 0)
            swap(res[0],res[1]);
        return res;
    }

    for (int cur = 0; cur < ns.routes[a_route]->time_table.size(); cur++) {
        if (ns.routes[a_route]->time_table[cur].task == a)
            a_pos = cur + 1;
        route_task_a.push_back(ns.routes[a_route]->time_table[cur].task);
    }

    for (int cur = 0; cur < ns.routes[v_route]->time_table.size(); cur++) {
        if (ns.routes[v_route]->time_table[cur].task == v)
            v_pos = cur + 1;
        route_task_v.push_back(ns.routes[v_route]->time_table[cur].task);
    }

    My_Assert(a_pos != -1 && v_pos != -1, "can't find chosen Task!");

    vector<int> buffer;
    std::move(route_task_a.begin() + a_pos, route_task_a.end(), back_inserter(buffer));
    route_task_a.erase(route_task_a.begin() + a_pos, route_task_a.end());

    std::move(route_task_v.begin() + v_pos, route_task_v.end(), back_inserter(route_task_a));
    route_task_v.erase(route_task_v.begin() + v_pos, route_task_v.end());

    std::move(buffer.begin(), buffer.end(), back_inserter(route_task_v));

    vector<int> time_tbl_u = mcgrp.cal_arrive_time(route_task_a);
    vector<int> time_tbl_i = mcgrp.cal_arrive_time(route_task_v);

    My_Assert(time_tbl_u.size() == route_task_a.size(), "wrong time table");
    My_Assert(time_tbl_i.size() == route_task_v.size(), "wrong time table");

    vector<RouteInfo::TimeTable> intermediate_a;
    vector<RouteInfo::TimeTable> intermediate_v;
    for (int k = 0; k < time_tbl_u.size(); k++) {
        if (!allow_infeasible && time_tbl_u[k] > mcgrp.inst_tasks[route_task_a[k]].time_window.second)
            return res;
        intermediate_a.push_back({route_task_a[k], time_tbl_u[k]});
    }

    for (int k = 0; k < time_tbl_i.size(); k++) {
        if (!allow_infeasible && time_tbl_i[k] > mcgrp.inst_tasks[route_task_v[k]].time_window.second)
            return res;
        intermediate_v.push_back({route_task_v[k], time_tbl_i[k]});
    }

    res[0] = intermediate_a;
    res[1] = intermediate_v;

    return res;
}

bool NewSwapEnds::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    // stub block
    return false;
}

bool NewSwapEnds::update_score(HighSpeedNeighBorSearch &ns)
{
    if(move_result.delta > 0) return false;

    // ...a-(seq_v)...
    // ...v-(seq_a)...
    const int a = move_result.task1;
    const int v = move_result.task2;

    //Extract the seq which needs swapped
    vector<int> output_a_seq = move_result.move_arguments_bak["output_a_seq"];
    vector<int> output_v_seq = move_result.move_arguments_bak["output_v_seq"];

    double penalty = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    if(!output_a_seq.empty())
        ns.score_matrix[max(0,ns.solution[output_a_seq.front()]->pre->ID)][max(0,output_a_seq.front())] += penalty;

    if(!output_v_seq.empty())
        ns.score_matrix[max(0,ns.solution[output_v_seq.front()]->pre->ID)][max(0,output_v_seq.front())] += penalty;

    return true;
}
