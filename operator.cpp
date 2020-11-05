//
// Created by luke on 2020/11/5.
//

#include "operator.h"
#include "NeighborSearch.h"

void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, MoveOperator &move_operator)
{
    vector<int> task_set(mcgrp.actual_task_num);
    std::generate(task_set.begin(), task_set.end(), Generator());
    mcgrp._rng.RandPerm(task_set);    //shuffle tasks

    auto original_policy = ns.policy.get();
    ns.policy.set(BEST_ACCEPT | DOWNHILL | DELTA_ONLY);
    ns.policy.beta = 0.5;
    ns.policy.tolerance = 0.003;
//    ns.neigh_size = mcgrp.neigh_size;

    int chosen_task = -1;
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
    ns.policy.set(original_policy);
    ns.policy.beta = 0;
    ns.policy.tolerance = 0;
}
