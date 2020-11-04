//
// Created by luke on 2019/12/31.
//

#include "SearchPolicy.h"

double Policy::getFitnessDelta(const MOVE &move_result) {
    My_Assert(beta != MAX(beta), "beta is undefined!");
    return move_result.delta + beta * (move_result.vio_load_delta + move_result.vio_time_delta);
}

bool Policy::check_move(const MOVE &move_result)
{
    //To control the frequency
    static int count = 0;
    static bool flag = true;



    /*----------------------------Feasible search policy------------------------------------*/
    if (has_rule(DELTA_ONLY)) {
        if ((has_rule(DOWNHILL)))
            return move_result.delta <= 0;
        else if (has_rule(TOLERANCE)) {
            if(move_result.delta < 0){
                return true;
            }

            else if (flag && move_result.new_total_route_length <= (1.0 + tolerance) * benchmark) {
                    flag = false;
                    return true;
                }
            else{
                count++;
                if (count == tabu_step_threshold) {
                    flag = true;
                    count = 0;
                }
                return false;
            }
        }
    }

    /*----------------------------Infeasible search policy------------------------------------*/
    else if (has_rule(FITNESS_ONLY)) {
        double delta_fitness = getFitnessDelta(move_result);

        if ((has_rule(DOWNHILL)))
            return delta_fitness <= 0;
        else if (has_rule(TOLERANCE)) {
            if (tolerance != 0 && delta_fitness <= tolerance * benchmark) {
                if (flag) {
                    flag = false;
                    return true;
                }
                else {
                    count++;
                    if (count == tabu_step_threshold) {
                        flag = true;
                        count = 0;
                    }
                    return false;
                }
            }
            else {
                return false;
            }
        }
    }
}

bool Policy::check_result(const MOVE &M1, const MOVE &M2)
{
    // We are only concerned about the savings (increase/decrese) in total length
    // This is the default approach
    // Decide in terms of total length only
    if (has_rule(DELTA_ONLY)) {
        return M1.new_total_route_length < M2.new_total_route_length;
    }
    else if (has_rule(FITNESS_ONLY)) {
        if(!M2.considerable){
            return true;
        }

        return getFitnessDelta(M1) < getFitnessDelta(M2);
    }
}