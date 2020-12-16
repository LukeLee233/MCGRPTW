//
// Created by luke on 2019/12/31.
//

#ifndef SEARCHPOLICY_H
#define SEARCHPOLICY_H

#include <bitset>
#include "MCGRP.h"
#include "config.h"


enum search_sequence
{
    RTRIDP, IDPRTR
};

//the number of rules
constexpr int RULES = 8;

//Search Rule
extern const bitset<RULES> UNDEFINE;

extern const bitset<RULES> DOWNHILL;
extern const bitset<RULES> TOLERANCE;

extern const bitset<RULES> FEASIBLE;
extern const bitset<RULES> INFEASIBLE;

//Accept Rule
extern const bitset<RULES> FIRST_ACCEPT;
extern const bitset<RULES> BEST_ACCEPT;

class HighSpeedNeighBorSearch;

class Policy{
    bitset<RULES> current_policy;    //the policy used by improvement search

    // parameters to control beta
    const int update_frequency = 5;
    int update_tick_time = update_frequency;
    int too_far = 0;
    int too_near = 0;

    //infeasible search parameters
    double beta = DBL_MAX;

private:
    double infeasible_distance_threshold;
    double time_window_multiplier = 1.5;
    double capacity_multiplier = 1.5;
public:
    void setTime_window_multiplier(double timeWindowMultiplier);
    void setCapacity_multiplier(double capacityMultiplier);

public:
    double getBeta() const;
    void setBeta(double beta);
    void clearContext(string mode);

    inline int get_pseudo_capacity(int capacity){
        // policy 1
        //    double scale = mcgrp.capacity * policy.beta;
        //    int pseudo_capacity = int(mcgrp.capacity *(1+ (scale / (policy.beta*sqrt(1+scale*scale)))));

        //policy 2
        //    int pseudo_capacity = mcgrp.capacity;

        //policy 3
        return int(capacity * capacity_multiplier);
    }

    inline int get_pseudo_time_window(int timestamp){
        return int(timestamp * time_window_multiplier);
    }

    bool check_time_window(const MCGRP& mcgrp,const vector<vector<RouteInfo::TimeTable>>& time_table);

    // oscillation parameters
    int tabu_time = 0;
    int continue_time = 0;
    int continue_threshold = 15;
    int tabu_step_threshold;
    double tolerance = 0;
    double benchmark = 0;

    void update_beta(double infeasible_distance);

    Policy(int tabu_step_):tabu_step_threshold(tabu_step_){
        current_policy = UNDEFINE;
    };

    const bitset<RULES> &getCurrent_policy() const;
    void setCurrent_policy(const bitset<RULES> &currentPolicy);

    /*!
 * test policy whether has a given rule
 * @param policy
 * @param rule
 * @return
 */
   inline bool has_rule(const bitset<RULES> &rule){
        return current_policy == (current_policy | rule);}

    /*!
    * @details check whether the given move is permissible based on the given policy
    * @param mcgrp
    * @param move_result
    * @param policy
    * @return
    */
    bool check_move(const MCGRP& mcgrp, HighSpeedNeighBorSearch& ns,const MoveResult &move_result);

    /*!
     * @details check if moveresult M1 is better than moveresult M2
     * @param ns
     * @param M1
     * @param M2
     * @return
     */
    bool check_result(const MCGRP& mcgrp, HighSpeedNeighBorSearch& ns, const MoveResult &M1, const MoveResult &M2);
};

#endif //SEARCHPOLICY_H
