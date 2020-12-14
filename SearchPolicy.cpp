//
// Created by luke on 2019/12/31.
//

#include "SearchPolicy.h"
#include "NeighborSearch.h"
#include "MCGRP.h"

//Search Rule
const bitset<RULES> UNDEFINE("00000000");                                //undefined policy

const bitset<RULES> DOWNHILL("00000001");
const bitset<RULES> TOLERANCE("00000010");

const bitset<RULES> FEASIBLE("10000000");                 //This means feasible search
const bitset<RULES> INFEASIBLE("01000000");              //This means infeasible search

//Accept Rule
const bitset<RULES> FIRST_ACCEPT("00010000");
const bitset<RULES> BEST_ACCEPT("00100000");



bool Policy::check_move(const MCGRP& mcgrp, HighSpeedNeighBorSearch& ns, const MoveResult &move_result)
{
    /*----------------------------Feasible search policy------------------------------------*/
    if (has_rule(FEASIBLE)) {
        if ((has_rule(DOWNHILL)))
            return move_result.delta <= 0;
        else if (has_rule(TOLERANCE)) {
            if (move_result.delta <= tolerance * benchmark) {
                tabu_time--;
                if(tabu_time <= 0){
                    continue_time--;
                    if(continue_time < 0){
                        continue_time = continue_threshold;
                        tabu_time = tabu_step_threshold;
                    }

                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
    }

    /*----------------------------Infeasible search policy------------------------------------*/
    else if (has_rule(INFEASIBLE)) {
        double delta_fitness = ns.getFitnessDelta(mcgrp,ns.policy, move_result);

        if ((has_rule(DOWNHILL)))
            return delta_fitness <= 0;
        else if (has_rule(TOLERANCE)) {
            if (move_result.delta <= tolerance * benchmark) {
                tabu_time--;
                if(tabu_time <= 0){
                    continue_time--;
                    if(continue_time < 0){
                        continue_time = continue_threshold;
                        tabu_time = tabu_step_threshold;
                    }

                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }
    }

    else{
        My_Assert(false, "error, unknown policy");
        return false;
    }
}

bool Policy::check_result(const MCGRP& mcgrp, HighSpeedNeighBorSearch& ns, const MoveResult &M1, const MoveResult &M2)
{
    // We are only concerned about the savings (increase/decrease) in total length
    // This is the default approach
    // Decide in terms of total length only
    if (has_rule(FEASIBLE)) {
        return M1.new_total_route_length < M2.new_total_route_length;
    }
    else if (has_rule(INFEASIBLE)) {
        if(!M2.considerable){
            return true;
        }

        return ns.getFitnessDelta(mcgrp,ns.policy, M1) < ns.getFitnessDelta(mcgrp,ns.policy, M2);
    }
}

void Policy::update_beta(double infeasible_distance)
{
    if (infeasible_distance > infeasible_distance_threshold)
        too_far++;
    else
        too_near++;

    update_tick_time--;
    if (update_tick_time <= 0) {
        update_tick_time = update_frequency;

        if (too_near == update_frequency) {
            DEBUG_PRINT("Infeasible search too near");
            beta /= 2;
        }
        else if (too_far == update_frequency) {
            DEBUG_PRINT("Infeasible search too far");
            beta *= 2;
        }

        too_far = 0;
        too_near = 0;
    }
}


bool Policy::check_time_window(const MCGRP &mcgrp, const vector<vector<RouteInfo::TimeTable>> &time_table)
{
    for(const auto& tbl : time_table){
        for(const auto & node : tbl){
            if(node.arrive_time > get_pseudo_time_window(mcgrp.inst_tasks[node.task].time_window.second))
                return false;
        }
    }
    return true;
}

void Policy::setBeta(double beta)
{
    Policy::beta = beta;
}

void Policy::clearContext(string mode)
{
    if(mode == "infeasible"){
        beta = DBL_MAX;
        too_far = 0;
        too_near = 0;
        update_tick_time = update_frequency;
    }
}

double Policy::getBeta() const
{
    return beta;
}

const bitset<RULES> &Policy::getCurrent_policy() const
{
    return current_policy;
}

void Policy::setCurrent_policy(const bitset<RULES> &currentPolicy)
{
    current_policy = currentPolicy;
}

void Policy::setTime_window_multiplier(double timeWindowMultiplier)
{
    time_window_multiplier = timeWindowMultiplier;
}

void Policy::setCapacity_multiplier(double capacityMultiplier)
{
    capacity_multiplier = capacityMultiplier;
}




