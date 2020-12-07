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
    ns.policy.setCurrent_policy(FIRST_ACCEPT | DOWNHILL | FEASIBLE);
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
viterbi::PseudoTask MoveOperator::Seq2Pseudo(const vector<int> &seq, int phase)
{
    // Auto-generated stub
    auto time_window = Task::Window(0,0);
    return viterbi::PseudoTask(0, 0, 0, 0, 0, time_window);
}


viterbi::BestDecode viterbi::viterbi_decode(const viterbi::PseudoTask *first,
                                            const viterbi::PseudoTask *last,
                                            const vector<int> &move_seq,
                                            const MCGRP &mcgrp,
                                            Policy &policy)
{
    // Auto-generated stub
    return viterbi::BestDecode();
}
