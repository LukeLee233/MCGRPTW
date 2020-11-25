#include "swapends.h"
#include <algorithm>

using namespace std;

struct RouteSegment
{
    // Contains information about a particular segment of a route.
    int segment_start;
    int segment_end;
    int num_custs;
    int load;
    double len;
};


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
    buffer.segment_end = ns.routes[route]->end;
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

    const int a_route = chosen_route;
    const int v_route = neighbor_route;

    My_Assert(a_route != v_route, "swap ends called with a and v in same route!");

    //get sequence that need to be swap
    RouteSegment seg_after_a = get_segment_info(mcgrp, ns, chosen_task);
    RouteSegment seg_after_v = get_segment_info(mcgrp, ns, neighbor_task);

    const auto a_load_delta = -seg_after_a.load + seg_after_v.load;
    const auto v_load_delta = -seg_after_v.load + seg_after_a.load;

    //handle load constraint
    double vio_load_delta = 0;
    if (ns.policy.has_rule(DELTA_ONLY)) {
        if (ns.routes[a_route]->load + a_load_delta > mcgrp.capacity
            || ns.routes[v_route]->load + v_load_delta > mcgrp.capacity) {
            return false;
        }
    }
    else if (ns.policy.has_rule(FITNESS_ONLY)) {
        My_Assert(v_load_delta == -a_load_delta, "Wrong arguments");

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


    double delta;

    //a can be dummy and v can be dummy too
    const int a = max(chosen_task, 0);
    const int v = max(neighbor_task, 0);

    bool allow_infeasible = ns.policy.has_rule(FITNESS_ONLY) ? true : false;
    vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};

    new_time_tbl = expected_time_table(ns,mcgrp,a,v,a_route,v_route,seg_after_a,seg_after_v,allow_infeasible);

    if(!mcgrp.isTimeTableFeasible(new_time_tbl[0])){
        move_result.reset();
        return false;
    }

    if (seg_after_a.num_custs == 0 && seg_after_v.num_custs == 0) {
        //means no actual move happens
        //...-a-[]-0
        //...-v-[]-0
        move_result.reset();
        return false;
    }
    else if (seg_after_a.num_custs == 0) {
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

    move_result.total_number_of_routes = ns.routes.activated_route_id.size();

    vector<int> a_seq;
    a_seq.push_back(chosen_task);
    int cur_task = ns.solution[chosen_task]->ID;
    for (auto i = 0; i < seg_after_a.num_custs; i++) {
        cur_task = ns.solution[cur_task]->next->ID;
        a_seq.push_back(cur_task);
    }

    vector<int> v_seq;
    v_seq.push_back(neighbor_task);
    cur_task = ns.solution[neighbor_task]->ID;
    for (auto i = 0; i < seg_after_v.num_custs; i++) {
        cur_task = ns.solution[cur_task]->next->ID;
        v_seq.push_back(cur_task);
    }


    move_result.move_arguments = a_seq;
    move_result.move_arguments.insert(move_result.move_arguments.end(), v_seq.begin(), v_seq.end());

    move_result.task1 = chosen_task;
    move_result.task2 = neighbor_task;

    move_result.vio_load_delta = vio_load_delta;

    move_result.considerable = true;

    move_result.route_time_tbl = new_time_tbl;
    move_result.vio_time_delta =
        mcgrp.get_vio_time(move_result.route_time_tbl[0])
            + mcgrp.get_vio_time(move_result.route_time_tbl[1])
            - mcgrp.get_vio_time(ns.routes[a_route]->time_table)
            - mcgrp.get_vio_time(ns.routes[v_route]->time_table);
    return true;
}

void NewSwapEnds::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    DEBUG_PRINT("execute a 2-opt:swap-ends move");

    My_Assert(move_result.considerable, "Invalid predictions");

    const int a = move_result.task1;
    const int v = move_result.task2;

    //Extract the seq which needs swapped
    vector<int> a_seq;
    vector<int> v_seq;

    My_Assert(a == move_result.move_arguments.front(), "Incorrect arguments");
    int i;
    for (i = 1; move_result.move_arguments[i] != v; i++) {
        a_seq.push_back(move_result.move_arguments[i]);
    }

    My_Assert(v == move_result.move_arguments[i], "Incorrect arguments");
    i++;
    for (; i < move_result.move_arguments.size(); i++) {
        v_seq.push_back(move_result.move_arguments[i]);
    }

    My_Assert(all_of(a_seq.begin(), a_seq.end(), [&](int i)
    { return i >= 1 && i <= mcgrp.actual_task_num; }), "Wrong Task");
    My_Assert(all_of(v_seq.begin(), v_seq.end(), [&](int i)
    { return i >= 1 && i <= mcgrp.actual_task_num; }), "Wrong Task");

    My_Assert(a_seq.size() == move_result.seq1_cus_num && v_seq.size() == move_result.seq2_cus_num, "Wrong Task");

    My_Assert(move_result.num_affected_routes == 2, "Wrong routes number");

    const int a_route = move_result.route_id[0];
    const int v_route = move_result.route_id[1];
    My_Assert(a_route != v_route, "Wrong routes!");
    My_Assert(ns.routes.activated_route_id.find(a_route) != ns.routes.activated_route_id.end(), "Invalid route");
    My_Assert(ns.routes.activated_route_id.find(v_route) != ns.routes.activated_route_id.end(), "Invalid route");

    int a_seq_edge_num = mcgrp.count_edges(a_seq);
    int v_seq_edge_num = mcgrp.count_edges(v_seq);
    ns.routes[a_route]->num_edges += (v_seq_edge_num - a_seq_edge_num);
    ns.routes[v_route]->num_edges += (a_seq_edge_num - v_seq_edge_num);

    if (move_result.seq1_cus_num != 0 && move_result.seq2_cus_num != 0) {
        //Modify routes info
        My_Assert(!a_seq.empty() && !v_seq.empty(), "Wrong arguments");
        My_Assert(ns.solution[a]->next->ID == a_seq.front(), "Wrong Task");
        My_Assert(ns.solution[v]->next->ID == v_seq.front(), "Wrong Task");

        ns.routes[a_route]->length = move_result.route_lens[0];
        ns.routes[v_route]->length = move_result.route_lens[1];

        ns.routes[a_route]->load = move_result.route_loads[0];
        ns.routes[v_route]->load = move_result.route_loads[1];

        ns.routes[a_route]->num_customers = move_result.route_custs_num[0];
        ns.routes[v_route]->num_customers = move_result.route_custs_num[1];

        ns.routes[a_route]->time_table = move_result.route_time_tbl[0];
        ns.routes[v_route]->time_table = move_result.route_time_tbl[1];

        for (auto task : a_seq) {
            ns.solution[task]->route_id = v_route;
        }

        for (auto task: v_seq) {
            ns.solution[task]->route_id = a_route;
        }

        if (ns.solution[a_seq.front()]->pre->ID < 0) {
            ns.routes[a_route]->start = v_seq.front();
        }

        ns.routes[a_route]->end = v_seq.back();

        if (ns.solution[v_seq.front()]->pre->ID < 0) {
            ns.routes[v_route]->start = a_seq.front();
        }

        ns.routes[v_route]->end = a_seq.back();

        //handle solution
        //swap two segments
        const auto a_final_dummy = ns.solution[ns.solution[a_seq.back()]->next->ID];
        const auto v_final_dummy = ns.solution[ns.solution[v_seq.back()]->next->ID];

        ns.solution[a]->next = ns.solution[v_seq.front()];
        ns.solution[v_seq.back()]->next = a_final_dummy;
        ns.solution[v_seq.front()]->pre = ns.solution[a];
        a_final_dummy->pre = ns.solution[v_seq.back()];

        ns.solution[v]->next = ns.solution[a_seq.front()];
        ns.solution[a_seq.back()]->next = v_final_dummy;
        ns.solution[a_seq.front()]->pre = ns.solution[v];
        v_final_dummy->pre = ns.solution[a_seq.back()];

    }
    else if (move_result.seq1_cus_num == 0) {
        My_Assert(move_result.seq2_cus_num != 0, "Wrong arguments");
        My_Assert(a_seq.empty() && !v_seq.empty(), "Wrong arguments");
        My_Assert(ns.solution[v]->next->ID == v_seq.front(), "Wrong Task");

        if (v_seq.front() == ns.routes[v_route]->start && v_seq.back() == ns.routes[v_route]->end) {
            //v route need to be eliminated
            //Modify routes info

            ns.routes[a_route]->length = move_result.route_lens[0];

            ns.routes[a_route]->load = move_result.route_loads[0];

            ns.routes[a_route]->num_customers = move_result.route_custs_num[0];

            ns.routes[a_route]->time_table = move_result.route_time_tbl[0];

            for (auto task: v_seq) {
                ns.solution[task]->route_id = a_route;
            }

            ns.routes[a_route]->end = v_seq.back();

            ns.routes.free_route(v_route);

            //handle solution
            //record next dummy Task of seg v
            int dummy_marker = 0;
            HighSpeedNeighBorSearch::TASK_NODE *reserve_dummy;
            bool very_end_case = false;
            if (ns.solution[v_seq.back()]->next == ns.solution.very_end) {
                dummy_marker = ns.solution[v_seq.front()]->pre->ID;
                reserve_dummy = ns.solution[v_seq.back()]->next;
                very_end_case = true;
            }
            else {
                dummy_marker = ns.solution[v_seq.back()]->next->ID;
                reserve_dummy = ns.solution[v_seq.front()]->pre;
            }

            //handle solution
            //swap two segments
            const auto a_final_dummy = ns.solution[ns.solution[a]->next->ID];


            ns.solution[a]->next = ns.solution[v_seq.front()];
            ns.solution[v_seq.back()]->next = a_final_dummy;
            ns.solution[v_seq.front()]->pre = ns.solution[a];
            a_final_dummy->pre = ns.solution[v_seq.back()];

            //free dummy marker
            if (very_end_case) {
                ns.solution[dummy_marker]->pre->next = reserve_dummy;
                reserve_dummy->pre = ns.solution[dummy_marker]->pre;
            }
            else {
                ns.solution[dummy_marker]->next->pre = reserve_dummy;
                reserve_dummy->next = ns.solution[dummy_marker]->next;
            }
            ns.solution.dummypool.free_dummy(dummy_marker);
        }
        else {
            //Modify routes info
            ns.routes[a_route]->length = move_result.route_lens[0];
            ns.routes[v_route]->length = move_result.route_lens[1];

            ns.routes[a_route]->load = move_result.route_loads[0];
            ns.routes[v_route]->load = move_result.route_loads[1];

            ns.routes[a_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[v_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[a_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[v_route]->time_table = move_result.route_time_tbl[1];

            for (auto task: v_seq) {
                ns.solution[task]->route_id = a_route;
            }

            ns.routes[a_route]->end = v_seq.back();

            ns.routes[v_route]->end = v;


            //handle solution
            //swap two segments
            const auto a_final_dummy = ns.solution[ns.solution[a]->next->ID];
            const auto v_final_dummy = ns.solution[ns.solution[v_seq.back()]->next->ID];


            ns.solution[a]->next = ns.solution[v_seq.front()];
            ns.solution[v_seq.back()]->next = a_final_dummy;
            ns.solution[v_seq.front()]->pre = ns.solution[a];
            a_final_dummy->pre = ns.solution[v_seq.back()];

            ns.solution[v]->next = v_final_dummy;
            v_final_dummy->pre = ns.solution[v];
        }

    }
    else if (move_result.seq2_cus_num == 0) {
        My_Assert(move_result.seq1_cus_num != 0, "Wrong arguments");
        My_Assert(!a_seq.empty() && v_seq.empty(), "Wrong arguments");
        My_Assert(ns.solution[a]->next->ID == a_seq.front(), "Wrong Task");

        if (a_seq.front() == ns.routes[a_route]->start && a_seq.back() == ns.routes[a_route]->end) {
            //a route need to be eliminated
            //Modify routes info

            ns.routes[v_route]->length = move_result.route_lens[1];

            ns.routes[v_route]->load = move_result.route_loads[1];

            ns.routes[v_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[v_route]->time_table = move_result.route_time_tbl[1];

            for (auto task: a_seq) {
                ns.solution[task]->route_id = v_route;
            }

            ns.routes[v_route]->end = a_seq.back();

            ns.routes.free_route(a_route);

            //handle solution
            //record next dummy Task of seg v
            int dummy_marker = 0;
            HighSpeedNeighBorSearch::TASK_NODE *reserve_dummy;
            bool very_end_case = false;
            if (ns.solution[a_seq.back()]->next == ns.solution.very_end) {
                dummy_marker = ns.solution[a_seq.front()]->pre->ID;
                reserve_dummy = ns.solution[a_seq.back()]->next;
                very_end_case = true;
            }
            else {
                dummy_marker = ns.solution[a_seq.back()]->next->ID;
                reserve_dummy = ns.solution[a_seq.front()]->pre;
            }

            //handle solution
            //swap two segments
            const auto v_final_dummy = ns.solution[ns.solution[v]->next->ID];

            ns.solution[v]->next = ns.solution[a_seq.front()];
            ns.solution[a_seq.back()]->next = v_final_dummy;
            ns.solution[a_seq.front()]->pre = ns.solution[v];
            v_final_dummy->pre = ns.solution[a_seq.back()];

            //free dummy marker
            if (very_end_case) {
                ns.solution[dummy_marker]->pre->next = reserve_dummy;
                reserve_dummy->pre = ns.solution[dummy_marker]->pre;
            }
            else {
                ns.solution[dummy_marker]->next->pre = reserve_dummy;
                reserve_dummy->next = ns.solution[dummy_marker]->next;
            }


            ns.solution.dummypool.free_dummy(dummy_marker);
        }
        else {
            //Modify routes info
            ns.routes[a_route]->length = move_result.route_lens[0];
            ns.routes[v_route]->length = move_result.route_lens[1];

            ns.routes[a_route]->load = move_result.route_loads[0];
            ns.routes[v_route]->load = move_result.route_loads[1];

            ns.routes[a_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[v_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[a_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[v_route]->time_table = move_result.route_time_tbl[1];

            for (auto task: a_seq) {
                ns.solution[task]->route_id = v_route;
            }

            ns.routes[v_route]->end = a_seq.back();

            ns.routes[a_route]->end = a;

            //handle solution
            //swap two segments
            const auto a_final_dummy = ns.solution[ns.solution[a_seq.back()]->next->ID];
            const auto v_final_dummy = ns.solution[ns.solution[v]->next->ID];


            ns.solution[v]->next = ns.solution[a_seq.front()];
            ns.solution[a_seq.back()]->next = v_final_dummy;
            ns.solution[a_seq.front()]->pre = ns.solution[v];
            v_final_dummy->pre = ns.solution[a_seq.back()];

            ns.solution[a]->next = a_final_dummy;
            a_final_dummy->pre = ns.solution[a];
        }

    }
    else {
        My_Assert(false, "Wrong arguments");
    }

    //modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp), "Prediction wrong!");

    if (move_result.delta == 0) {
        ns.equal_step++;
    }

    update_score(ns);
    move_result.reset();
    ns.search_step++;
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
    vector<int> a_seq;
    vector<int> v_seq;

    My_Assert(a == move_result.move_arguments.front(), "Incorrect arguments");
    int i;
    for (i = 1; move_result.move_arguments[i] != v; i++) {
        a_seq.push_back(move_result.move_arguments[i]);
    }

    My_Assert(v == move_result.move_arguments[i], "Incorrect arguments");
    i++;
    for (; i < move_result.move_arguments.size(); i++) {
        v_seq.push_back(move_result.move_arguments[i]);
    }

    double penalty = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    if(!a_seq.empty())
        ns.score_matrix[max(0,ns.solution[a_seq.front()]->pre->ID)][max(0,a_seq.front())] += penalty;

    if(!v_seq.empty())
        ns.score_matrix[max(0,ns.solution[v_seq.front()]->pre->ID)][max(0,v_seq.front())] += penalty;

    return true;
}
