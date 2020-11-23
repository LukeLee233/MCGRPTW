//
// Created by luke on 2020/11/25.
//

#include "attraction.h"
#include "NeighborSearch.h"
#include "swapends.h"

bool Attraction::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

#ifdef DEBUG
    attempt_count = 0;
    hit_count = 0;
#endif

    if(!prob_valid(ns.prob_matrix[chosen_task])){
        DEBUG_PRINT("warning, probability is not available currently.");
        return false;
    }


    MoveResult BestM;
    unordered_set<int> sample_set;

    for(int sample_count = 0; sample_count < sample_times;sample_count++){
        int task = sample::uniform_sample(ns.prob_matrix[chosen_task]);
        if(task == DUMMY || ns.solution[task]->next != nullptr)
            sample_set.insert(task);
    }

    ns.search_space = vector<int>(sample_set.begin(), sample_set.end());

    for(auto neighbor_task : ns.search_space) {
        My_Assert(chosen_task != neighbor_task, "error, neighbor Task can't be itself!");
        if (mcgrp.is_edge(chosen_task)) {
            My_Assert(mcgrp.inst_tasks[chosen_task].inverse != neighbor_task, "error, neighbor Task can't be itself!");
        }

#ifdef DEBUG
        attempt_count++;
#endif

        if (considerable_move(ns, mcgrp, chosen_task, neighbor_task) && ns.policy.check_move(move_result)) {
            if (ns.policy.has_rule(FIRST_ACCEPT)) {
                move(ns, mcgrp);
                return true;
            }
            else if (ns.policy.has_rule(BEST_ACCEPT)) {
                if (ns.policy.check_result(move_result, BestM))
                    BestM = move_result;
            }
            else {
                My_Assert(false, "Unknown accept rule!");
            }
        }
    }

    if (ns.policy.has_rule(FIRST_ACCEPT)) {
        DEBUG_PRINT("No actual move: First Accept Rule");
        return false;
    }
    else if (ns.policy.has_rule(BEST_ACCEPT)) {
        if (BestM.considerable == false) {
            DEBUG_PRINT("No actual move: Best Accept Rule");
            return false;
        }
        else {
            move_result = BestM;
            move(ns, mcgrp);
            return true;
        }
    }
    else {
        My_Assert(false, "Unknown accept rule");
    }
}

/*!
 * attract task u behind task a
 *  ...a-u...
 * @param ns
 * @param mcgrp
 * @param a
 * @param u
 * @return
 */
bool Attraction::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int a, const int u)
{
    // Task a cannot be dummy Task
    My_Assert(a >= 1 && a <= mcgrp.actual_task_num,"Wrong Task");

    if(max(ns.solution[a]->next->ID,0) == u){
        // Nothing to do
        move_result.reset();
        return false;
    }

    if(u == DUMMY){
        if(post_slice.considerable_move(ns,mcgrp,a)){
            move_result = post_slice.move_result;
            return true;
        }else{
            move_result.reset();
            return false;
        }
    }else{
        if(post_insert.considerable_move(ns,mcgrp,{u}, a)){
            move_result = post_insert.move_result;
            return true;
        }else{
            move_result.reset();
            return false;
        }
    }

    My_Assert(false,"Cannot reach here!");
}


Attraction::Attraction(int sampleTimes)
    : sample_times(sampleTimes), post_slice(Postslice()), post_insert(XPostInsert(1))
{
    move_result.move_type = ATTRACTION;
}

void Attraction::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{

#ifdef DEBUG
    hit_count++;
#endif

    DEBUG_PRINT("execute a attraction move");

    My_Assert(move_result.considerable, "Invalid predictions");

    if (move_result.move_type == NeighborOperator::POST_SLICE) {
        post_slice.move_result = move_result;
        post_slice.move(ns, mcgrp);
        ns.trace(mcgrp);
    }else{
        post_insert.move_result = move_result;
        post_insert.move(ns,mcgrp);
    }

    move_result.reset();
    move_result.move_type = NeighborOperator::ATTRACTION;
}
