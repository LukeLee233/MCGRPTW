#include "SingleInsert.h"
#include "Presert.h"
#include "PostSert.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

bool SingleInsert::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

    presert_times = 0;
    postsert_times = 0;

    MoveResult BestM;

    ns.create_search_neighborhood(mcgrp, chosen_task);

    int b = chosen_task;


    for(auto neighbor_task : ns.search_space) {
        My_Assert(neighbor_task != b, "neighbor Task can't be itself!");

        if (neighbor_task != DUMMY) {
            //j can't be dummy and b can't be dummy neither here
            int j = neighbor_task;

            if (considerable_move(ns, mcgrp, b, j)) {
                if(ns.policy.has_rule(FIRST_ACCEPT)){
                    move(ns,mcgrp);
                    return true;
                }
                else if (ns.policy.has_rule(BEST_ACCEPT)){
                    if (ns.policy.check_result(move_result,BestM))
                        BestM = move_result;
                }
                else{
                    My_Assert(false,"Unknown accept rule!");
                }
            }
        }
        else {
            DEBUG_PRINT("Neighbor Task is dummy Task");
            //j is dummy here but b can't be dummy
            //each start and end location of each route will be considered
            //total 2 x route_nums cases
            int current_start = ns.solution.very_start->next->ID;

            while (ns.solution[current_start]->next != ns.solution.very_end) {
                // Consider the start location
                int j = current_start;

                if (b != j) { //not the same Task
                    if (considerable_move(ns, mcgrp, b, j)) {
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

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.routes[current_route]->end;
                j = current_end;
                if (b != j) {
                    if (considerable_move(ns, mcgrp, b, j)) {
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

                // Advance to next route's starting Task
                current_start = ns.solution[current_end]->next->next->ID;
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

    DEBUG_PRINT("Single insert presert times:" + to_string(presert_times));
    DEBUG_PRINT("Single insert postsert times:" + to_string(postsert_times));
}

bool SingleInsert::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int b,const int j)
{
    My_Assert(j>=1 && j<=mcgrp.actual_task_num,"Wrong Task");
    My_Assert(b>=1 && b <=mcgrp.actual_task_num,"Wrong Task");
    My_Assert(b != j,"You cannot insert the same Task!");


    const int a = max(ns.solution[b]->pre->ID , 0);
    const int c = max(ns.solution[b]->next->ID , 0);

    if (c == j) {
        // ...a-b-c(j)-k...
        //no need to presert
        postsert_times++;
        if (postsert.considerable_move(ns, mcgrp, b, j) && ns.policy.check_move(postsert.move_result)) {
            move_result = postsert.move_result;
            return true;
        }
        else {
            move_result.reset();
            move_result.move_type = NeighborOperator::SINGLE_INSERT;
            return false;
        }
    }
    else if (a == j) {
        // ...i-j(a)-b-c
        //no need to do postsert
        presert_times++;
        if (presert.considerable_move(ns, mcgrp, b, j) && ns.policy.check_move(presert.move_result)) {
            move_result = presert.move_result;
            return true;
        }
        else {
            move_result.reset();
            move_result.move_type = NeighborOperator::SINGLE_INSERT;
            return false;
        }
    }
    else{
        //...a-b-c...
        //...i-j-k...

        //check presert and postsert
        presert.considerable_move(ns, mcgrp, b, j);
        postsert.considerable_move(ns, mcgrp, b, j);

        if(presert.move_result.considerable && postsert.move_result.considerable){
            //Both considerable
            if(presert.move_result.delta <= postsert.move_result.delta){
                presert_times++;
                if(ns.policy.check_move(presert.move_result)){
                    move_result = presert.move_result;
                    return true;
                }
                else{
                    move_result.reset();
                    move_result.move_type = NeighborOperator::SINGLE_INSERT;
                    return false;
                }
            }
            else{
                postsert_times++;
                if(ns.policy.check_move(postsert.move_result)){
                    move_result = postsert.move_result;
                    return true;
                }
                else{
                    move_result.reset();
                    move_result.move_type = NeighborOperator::SINGLE_INSERT;
                    return false;
                }
            }
        }
        else if (presert.move_result.considerable){
            //only presert is considerable
            presert_times++;
            if(ns.policy.check_move(presert.move_result)){
                move_result = presert.move_result;
                return true;
            }
            else{
                move_result.reset();
                move_result.move_type = NeighborOperator::SINGLE_INSERT;
                return false;
            }
        }
        else if (postsert.move_result.considerable){
            postsert_times++;
            if(ns.policy.check_move(postsert.move_result)){
                move_result = postsert.move_result;
                return true;
            }
            else{
                move_result.reset();
                move_result.move_type = NeighborOperator::SINGLE_INSERT;
                return false;
            }
        }
        else{
            //No considerable
            move_result.reset();
            move_result.move_type = NeighborOperator::SINGLE_INSERT;
            return false;
        }
    }

    My_Assert(false,"Cannot reach here!");
}

void SingleInsert::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{
    DEBUG_PRINT("execute a single insert move");

    My_Assert(move_result.considerable,"Invalid predictions");

    if (move_result.move_type == NeighborOperator::PRESERT) {
        presert.move_result = move_result;
        presert.move(ns, mcgrp);
    }
    else if (move_result.move_type == NeighborOperator::POSTSERT) {
        postsert.move_result = move_result;
        postsert.move(ns, mcgrp);
    }
    else{
        My_Assert(false,"Unknown operator");
    }

    ns.trace(mcgrp);

//    mcgrp.check_best_infeasible_solution(ns.cur_solution_cost,ns.policy.beta,ns.total_vio_load,ns.negative_coding_sol);

    move_result.reset();
    move_result.move_type = NeighborOperator::SINGLE_INSERT;
}

