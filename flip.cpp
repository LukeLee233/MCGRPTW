#include "flip.h"
#include <algorithm>

using namespace std;


/*
 *
 */

bool Flip::considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, int start_task, int end_task){

    call_times++;
    move_result.reset();

    vector<int> move_seq = get_sequence(ns, start_task, end_task);

    //No need flipping
    if(move_seq.size() <= 1){
        return false;
    }

    My_Assert(all_of(move_seq.begin(), move_seq.end(), [&](int i){return i>=1 && i<=mcgrp.actual_task_num;}), "Wrong Task");
    My_Assert(ns.solution[move_seq.front()]->route_id == ns.solution[move_seq.back()]->route_id, "Flip attempted using different routes!");

    const int route_id = ns.solution[move_seq.front()]->route_id;

    if(ns.policy.has_rule(BEST_ACCEPT) && mcgrp.count_edges(move_seq) > 0){
        // with mode refine
        vector<int> first_seq;
        vector<int> last_seq;

        vector<int> route_seq = ns.get_sub_seq(ns.routes[ns.solution[move_seq.front()]->route_id]->start, ns.routes[ns.solution[move_seq.front()]->route_id]->end);

        int move_seq_loc = -1;
        for(int ii = 0;ii < route_seq.size(); ii++){
            if(route_seq[ii] == move_seq.front()) move_seq_loc = ii;
        }

        first_seq.insert(first_seq.end(),route_seq.begin(),route_seq.begin() + move_seq_loc);
        last_seq.insert(last_seq.end(), route_seq.begin() + move_seq_loc + move_seq.size(),route_seq.end());

        My_Assert(first_seq.front() < 0 && last_seq.back() < 0, "error, wrong sequence!");

        viterbi::PseudoTask first_task = Seq2Pseudo(mcgrp,first_seq, 1);
        viterbi::PseudoTask last_task = Seq2Pseudo(mcgrp,last_seq, 2);

        // insert sequence will invoke viterbi result
        reverse(move_seq.begin(),move_seq.end());
        auto best_sequence = viterbi::viterbi_decode(&first_task, &last_task,move_seq,mcgrp,ns.policy);

        // use best sequence to update the move result
        if(best_sequence.cost == INT32_MAX){
            return false;
        }

        move_result.choose_tasks(start_task, end_task);
        move_result.move_arguments["input_seq"] = vector<int>(move_seq.rbegin(), move_seq.rend());
        move_result.move_arguments["output_seq"] = vector<int>(best_sequence.seq.rbegin() + 1, best_sequence.seq.rend() - 1);

        double delta = best_sequence.cost - ns.routes[route_id]->length;

        move_result.num_affected_routes = 1;

        move_result.route_id.push_back(route_id);
        move_result.delta = delta;

        move_result.route_lens.push_back(ns.routes[route_id]->length + move_result.delta);
        move_result.route_loads.push_back(ns.routes[route_id]->load);

        move_result.route_custs_num.push_back(ns.routes[route_id]->num_customers);
        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;
        move_result.considerable = true;

        vector<int> new_route_seq;
        new_route_seq.insert(new_route_seq.end(),first_seq.begin() + 1,first_seq.end());
        new_route_seq.insert(new_route_seq.end(),best_sequence.seq.begin() + 1,best_sequence.seq.end() - 1);
        new_route_seq.insert(new_route_seq.end(),last_seq.begin(),last_seq.end() - 1);

        move_result.route_time_tbl.emplace_back(RouteInfo::TimeTable::zip(new_route_seq, mcgrp.cal_arrive_time(new_route_seq)));

        if(ns.policy.has_rule(INFEASIBLE)){
            auto old_time_info = mcgrp.get_vio_time(ns.routes[route_id]->time_table);
            auto new_time_info = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = new_time_info.second - old_time_info.second;
            move_result.vio_time_custom_num_delta = new_time_info.first - old_time_info.first;
        }else{
            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = 0;
            move_result.vio_time_custom_num_delta = 0;
        }

        return true;
    }
    else{
        // without mode refine

        vector<RouteInfo::TimeTable> new_time_tbl{{{-1, -1}}};
        bool allow_infeasible = ns.policy.has_rule(INFEASIBLE) ? true : false;

        new_time_tbl = expected_time_table(ns, mcgrp, move_seq, allow_infeasible);

        if(!mcgrp.isTimeTableFeasible(new_time_tbl)){
            return false;
        }

        if(allow_infeasible){
            vector<vector<RouteInfo::TimeTable>> pack_new_time_tbl;
            pack_new_time_tbl.push_back(new_time_tbl);
            if(ns.policy.check_time_window(mcgrp,pack_new_time_tbl)){
                return false;
            }
        }

        move_result.choose_tasks(start_task, end_task);
        move_result.move_arguments["input_seq"] = move_seq;
        move_result.move_arguments["output_seq"] = move_seq;

        double delta = 0;
        start_task = max(start_task,0);
        end_task = max(end_task,0);

        delta -= mcgrp.min_cost[mcgrp.inst_tasks[start_task].tail_node][mcgrp.inst_tasks[move_seq.front()].head_node];
        for(int cur = 0; cur < move_seq.size() - 1; cur++){
            delta -= mcgrp.inst_tasks[move_seq[cur]].serv_cost;
            delta -= mcgrp.min_cost[mcgrp.inst_tasks[move_seq[cur]].tail_node][mcgrp.inst_tasks[move_seq[cur+1]].head_node];
        }
        delta -= mcgrp.inst_tasks[move_seq.back()].serv_cost;
        delta -= mcgrp.min_cost[mcgrp.inst_tasks[move_seq.back()].tail_node][mcgrp.inst_tasks[end_task].head_node];

        reverse(move_seq.begin(), move_seq.end());

        delta += mcgrp.min_cost[mcgrp.inst_tasks[start_task].tail_node][mcgrp.inst_tasks[move_seq.front()].head_node];

        for(int cur = 0; cur<move_seq.size()-1; cur++){
            delta += mcgrp.inst_tasks[move_seq[cur]].serv_cost;
            delta += mcgrp.min_cost[mcgrp.inst_tasks[move_seq[cur]].tail_node][mcgrp.inst_tasks[move_seq[cur + 1]].head_node];
        }

        delta += mcgrp.inst_tasks[move_seq.back()].serv_cost;
        delta += mcgrp.min_cost[mcgrp.inst_tasks[move_seq.back()].tail_node][mcgrp.inst_tasks[end_task].head_node];


        move_result.num_affected_routes = 1;

        move_result.route_id.push_back(route_id);
        move_result.delta = delta;

        move_result.route_lens.push_back(ns.routes[route_id]->length + move_result.delta);
        move_result.route_loads.push_back(ns.routes[route_id]->load);

        move_result.route_custs_num.push_back(ns.routes[route_id]->num_customers);
        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;
        move_result.considerable = true;

        move_result.route_time_tbl.emplace_back(new_time_tbl);

        if(ns.policy.has_rule(INFEASIBLE)){
            auto old_time_info = mcgrp.get_vio_time(ns.routes[route_id]->time_table);
            auto new_time_info = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = new_time_info.second - old_time_info.second;
            move_result.vio_time_custom_num_delta = new_time_info.first - old_time_info.first;
        }else{
            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = 0;
            move_result.vio_time_custom_num_delta = 0;
        }

        return true;
    }

}

void Flip::move(LocalSearch &ns, const MCGRPTW &mcgrp){

    success_times++;

    DEBUG_PRINT("execute a 2-opt: flip move");

    My_Assert(move_result.considerable,"Invalid predictions");

    //phase 0: extract move arguments
    vector<int> input_seq = move_result.move_arguments["input_seq"];
    vector<int> output_seq = move_result.move_arguments["output_seq"];
    My_Assert(all_of(input_seq.begin(),input_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong Task");
    My_Assert(all_of(output_seq.begin(),output_seq.end(),[&](int i){return i>=1 && i<=mcgrp.actual_task_num;}),"Wrong Task");
    My_Assert(ns.solution[move_result.task1]->next->ID == input_seq.front() && ns.solution[move_result.task2]->pre->ID == input_seq.back(),"Wrong Task");


    // phase 1: update the route level info
    My_Assert(move_result.num_affected_routes == 1,"Wrong result info");
    const int route_id = move_result.route_id[0];

    ns.routes[route_id]->length = move_result.route_lens[0];

    ns.routes[route_id]->time_table = move_result.route_time_tbl[0];



    // phase 3: update the sequence level info
    const int start_task = move_result.task1;
    const int end_task = move_result.task2;

    if(input_seq == output_seq){
        //handle solution
        for(auto task : input_seq){
            auto tmp = ns.solution[task]->pre;
            ns.solution[task]->pre = ns.solution[task]->next;
            ns.solution[task]->next = tmp;
        }

        ns.solution[start_task]->next = ns.solution[input_seq.back()];
        ns.solution[end_task]->pre = ns.solution[input_seq.front()];
        ns.solution[input_seq.front()]->next = ns.solution[end_task];
        ns.solution[input_seq.back()]->pre = ns.solution[start_task];
    }
    else{
        for(auto task : input_seq){
            ns.solution[task]->clear();
        }

        // extra reverse needed here
        reverse(output_seq.begin(),output_seq.end());
        ns.solution.connect_tasks(output_seq);

        ns.solution[start_task]->next = ns.solution[output_seq.front()];
        ns.solution[end_task]->pre = ns.solution[output_seq.back()];
        ns.solution[output_seq.back()]->next = ns.solution[end_task];
        ns.solution[output_seq.front()]->pre = ns.solution[start_task];
    }

    for(auto task : output_seq){
        ns.solution[task]->route_id = route_id;
    }

    // phase 4: modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");


    update_score(ns);

    ns.trace(mcgrp);
    move_result.reset();
}

vector<int> Flip::get_sequence(LocalSearch &ns, const int start, const int end)
{
    // return sequence between start and end task, exclude
    vector<int> buffer;
    int cur = ns.solution[start]->next->ID;

    while (cur != end){
        buffer.push_back(cur);
        cur = ns.solution[cur]->next->ID;
    }

    return buffer;
}

vector<RouteInfo::TimeTable>
Flip::expected_time_table(LocalSearch &ns,
                          const MCGRPTW &mcgrp,
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

bool Flip::search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task)
{
    // stub block
    return false;
}

void Flip::_apply_reward(LocalSearch &ns,double reward)
{
    // ...start-(actual flip seq)-end...
    const int pre_seq = ns.solution[move_result.move_arguments.at("output_seq").front()]->pre->ID;
    const int after_seq = ns.solution[move_result.move_arguments.at("output_seq").back()]->next->ID;

    // update old
    const vector<int>& input_seq = move_result.move_arguments.at("input_seq");
    ns.score_matrix[pre_seq][input_seq.front()] -= reward;
    for(int idx = 1; idx < input_seq.size();idx++){
        ns.score_matrix[input_seq[idx-1]][input_seq[idx]] -= reward;
    }
    ns.score_matrix[input_seq.back()][after_seq] -= reward;


    // update new
    vector<int> output_seq = move_result.move_arguments.at("output_seq");
    reverse(output_seq.begin(),output_seq.end());
    ns.score_matrix[pre_seq][output_seq.front()] += reward;
    for(int idx = 1; idx < output_seq.size();idx++){
        ns.score_matrix[output_seq[idx-1]][output_seq[idx]] += reward;
    }
    ns.score_matrix[output_seq.back()][after_seq] += reward;
}
