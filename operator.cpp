//
// Created by luke on 2020/11/5.
//

#include "operator.h"
#include "NeighborSearch.h"
#include "SearchPolicy.h"

void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, MoveOperator &move_operator)
{
    vector<int> task_set(mcgrp.actual_task_num);
    std::generate(task_set.begin(), task_set.end(), Generator());
    mcgrp._rng.RandPerm(task_set);    //shuffle tasks

    auto original_policy = ns.policy.getCurrent_policy();
    ns.policy.setCurrent_policy(BEST_ACCEPT | DOWNHILL | FEASIBLE);
    ns.policy.setBeta(0.5);
    ns.policy.tolerance = 0.003;
    ns.neigh_size = mcgrp.neigh_size;

    int chosen_task = -1;

#ifdef DEBUG
    move_operator.attempt_count = 0;
    move_operator.hit_count = 0;
#endif

    for (int i = 0; i < mcgrp.actual_task_num; i++) {
        chosen_task = task_set[i];

        if (ns.solution[chosen_task]->next == nullptr){
            if (mcgrp.is_edge(chosen_task)) {
                chosen_task = mcgrp.inst_tasks[chosen_task].inverse;
                My_Assert(ns.solution[chosen_task]->next != nullptr,"An edge Task has been missed");
            }
            else {
                My_Assert(false,"A non edge Task has been missed!");
            }
        }

        move_operator.search(ns, mcgrp, chosen_task);
    }

    ns.neigh_size = 0;
    ns.policy.setCurrent_policy(original_policy);
    ns.policy.setBeta(0);
    ns.policy.tolerance = 0;
}

MoveOperator::MoveOperator()
    : hit_count(0), attempt_count(0){}


vector<int> MoveOperator::getHitInfo()
{
    return {attempt_count, hit_count};
}

vector<vector<int>> MoveOperator::generate_possible_seq(const MCGRP &mcgrp, const vector<int>& seq)
{
    vector<vector<int>> total_sequence{{}};
    vector<int> candidate_serve_mode;
    for(const int current_task : seq){
        vector<vector<int>> buffer;
        buffer.swap(total_sequence);

        candidate_serve_mode.clear();
        if(mcgrp.is_edge(current_task)){
            candidate_serve_mode.push_back(current_task);
            candidate_serve_mode.push_back(mcgrp.inst_tasks[current_task].inverse);
        }else{
            candidate_serve_mode.push_back(current_task);
        }

        for(int tail_task : candidate_serve_mode){
            for(const auto& head_seq : buffer){
                total_sequence.push_back(head_seq);
                total_sequence.back().push_back(tail_task);
            }
        }
    }

    return total_sequence;
}

viterbi::PseudoTask MoveOperator::Seq2Pseudo(const MCGRP &mcgrp, const vector<int> &seq, int phase)
{
    vector<int> regularized_seq = _seq_regularization(seq);
    int head_node = mcgrp.inst_tasks[regularized_seq.front()].head_node;
    int tail_node = mcgrp.inst_tasks[regularized_seq.back()].tail_node;
    int demand = 0;
    int serve_cost = 0;

    for(int ii = 0; ii < regularized_seq.size(); ii++){
        demand += mcgrp.inst_tasks[regularized_seq[ii]].demand;

        serve_cost += mcgrp.inst_tasks[regularized_seq[ii]].serv_cost;

        if(ii != 0){
            serve_cost += mcgrp.min_cost[mcgrp.inst_tasks[regularized_seq[ii-1]].tail_node][mcgrp.inst_tasks[regularized_seq[ii]].head_node];
        }
    }

    Task::Window time_window;
    if(phase == 1){
        int depart_time = mcgrp.cal_arrive_time(regularized_seq).back() + mcgrp.inst_tasks[regularized_seq.back()].serve_time;
        time_window.first = time_window.second = depart_time;
    }
    else if(phase == 2){
        time_window.first = mcgrp.inst_tasks[regularized_seq.front()].time_window.first;
        int x = mcgrp.inst_tasks[regularized_seq.front()].time_window.second;
        int constant = 0;

        for(int ii = 1;ii < regularized_seq.size();ii++){
            constant += mcgrp.inst_tasks[regularized_seq[ii-1]].serve_time;
            constant += mcgrp.min_time[mcgrp.inst_tasks[regularized_seq[ii-1]].tail_node][mcgrp.inst_tasks[regularized_seq[ii]].head_node];
            x = min(x, mcgrp.inst_tasks[regularized_seq[ii]].time_window.second - constant);
        }

        time_window.second = x;
    }else{
        My_Assert(false, "error, wrong phase");
    }

    return viterbi::PseudoTask(head_node,tail_node,demand,serve_cost,time_window);
}

vector<int> MoveOperator::_seq_regularization(const vector<int> seq)
{
    vector<int> res;
    for(int task : seq){
        res.push_back(max(DUMMY, task));
    }

    return res;
}

void MoveOperator::update_stable_likelihood(const MCGRP &mcgrp, HighSpeedNeighBorSearch &ns ,const vector<int>& output_seq,const MoveResult& move_result)
{
    if(ns.total_vio_load == 0 && ns.total_vio_time == 0 && ns.policy.has_rule(FEASIBLE)){
        int sign = 1;
        if(move_result.delta < 0){
            sign = -1;
        }

        for(const auto& task_id : output_seq){
            My_Assert(task_id > 0, "error, wrong state!");
            mcgrp.inst_tasks[task_id].move_time.total_time -= move_result.delta;
//            if(sign == 1) mcgrp.inst_tasks[task_id].move_time.up_time += 1;
//            else mcgrp.inst_tasks[task_id].move_time.down_time += 1;
        }

        if(!output_seq.empty()){
            int before = max(ns.solution[output_seq.front()]->pre->ID,0);
            int after = max(ns.solution[output_seq.back()]->next->ID,0);
//            if(sign == 1){
//                if(before > 0)  mcgrp.inst_tasks[before].move_time.up_time += 1;
//                if(after > 0)  mcgrp.inst_tasks[after].move_time.up_time += 1;
//            }
//            else{
//                if(before > 0)  mcgrp.inst_tasks[before].move_time.down_time += 1;
//                if(after > 0)  mcgrp.inst_tasks[after].move_time.down_time += 1;
//            }

            if(before > 0)  mcgrp.inst_tasks[before].move_time.total_time -= move_result.delta;
            if(after > 0)  mcgrp.inst_tasks[after].move_time.up_time -= move_result.delta;

        }
    }

}

viterbi::BestDecode viterbi::viterbi_decode(const viterbi::PseudoTask* first,
                                            const viterbi::PseudoTask* last,
                                            const vector<int> &move_seq,
                                            const MCGRP &mcgrp,
                                            Policy &policy)
{

    viterbi::BestDecode ans;

    vector<vector<const Task *>> DAG;
    DAG.push_back({&mcgrp.inst_tasks[DUMMY]});
    if(first != nullptr) DAG.push_back( {first} );
    for(int i = 0; i< move_seq.size(); i++){
        if(mcgrp.is_edge(move_seq[i])){
            DAG.push_back( {&mcgrp.inst_tasks[move_seq[i]], &mcgrp.inst_tasks[mcgrp.inst_tasks[move_seq[i]].inverse]} );
        }else{
            DAG.push_back( { &mcgrp.inst_tasks[move_seq[i]] } );
        }
    }
    if(last != nullptr) DAG.push_back( {last} );
    DAG.push_back({&mcgrp.inst_tasks[DUMMY]});

    vector<vector<int>> W(DAG.size(),vector<int>());
    vector<vector<int>> P(DAG.size(),vector<int>());
    vector<vector<int>> ArriveTime(DAG.size(), vector<int>());

    W[0].push_back(0);
    P[0].push_back(0);
    ArriveTime[0].push_back(0);

    function<pair<int,int>(const vector<int>&)> argmin_ = [](const vector<int>& seq){
        int index = -1;
        int val = INT32_MAX;
        for(int ii = 0;ii < seq.size();ii++){
            if (seq[ii] < val){
                index  = ii;
                val = seq[ii];
            }
        }

        return make_pair(index, val);
    };

    vector<int> distance_;
    vector<int> ArriveTime_;
    for(int i = 1; i < DAG.size();i++){

        for(int j = 0;j < DAG[i].size() ;j++){
            distance_.clear();
            ArriveTime_.clear();

            int LatestDepartureTime;

            if(policy.has_rule(FEASIBLE)){
                LatestDepartureTime = DAG[i][j]->time_window.second;
            }else if(policy.has_rule(INFEASIBLE)){
                LatestDepartureTime = policy.get_pseudo_time_window(DAG[i][j]->time_window.second);
            }else{
                My_Assert(false, "error, wrong arguments!");
            }

            for(int k = 0; k < W[i-1].size(); k++){
                if(W[i - 1][k] == INT32_MAX ||
                    ArriveTime[i - 1][k] == INT32_MAX ||
                    ArriveTime[i - 1][k] + DAG[i - 1][k]->serve_time +
                        mcgrp.min_time[DAG[i - 1][k]->tail_node][DAG[i][j]->head_node]
                        > LatestDepartureTime){
                    distance_.push_back(INT32_MAX);
                    ArriveTime_.push_back(INT32_MAX);
                }
                else{
                    distance_.push_back(W[i - 1][k] +
                        DAG[i - 1][k]->serv_cost +
                        mcgrp.min_cost[DAG[i - 1][k]->tail_node][DAG[i][j]->head_node]);
                    ArriveTime_.push_back(max(DAG[i][j]->time_window.first,
                                              ArriveTime[i - 1][k] + DAG[i - 1][k]->serve_time +
                                                  mcgrp.min_time[DAG[i-1][k]->tail_node][DAG[i][j]->head_node]));
                }
            }

            int idx;
            int val;
            auto tmp = argmin_(distance_);
            idx = tmp.first;
            val = tmp.second;

            P[i].push_back(idx);
            W[i].push_back(val);
            if (idx == -1)
                ArriveTime[i].push_back(INT32_MAX);
            else
                ArriveTime[i].push_back(ArriveTime_[idx]);
        }
    }

    if (P.back() == vector<int>{-1}){
        ans.cost = INT_MAX;
    }else{
        My_Assert(W.back().size() == 1, "Wrong status!");
        ans.cost = W.back().back();
        vector<int> indices{P.back().back()};
        for(int ii = (int)P.size() - 2; ii >= 2 ; ii--){
            indices.push_back(P[ii][indices.back()]);
        }
        My_Assert(indices.size() == P.size() - 2, "Wrong decode sequence!");

        reverse(indices.begin(),indices.end());
        for(int ii = 0;ii < (int)indices.size(); ii++){
            ans.seq.push_back(DAG[ii+1][indices[ii]]->task_id);
        }
    }

    return ans;
}
