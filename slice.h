//
// Created by luke on 2020/1/26.
// doesn't support infeasible and descent search
//

#ifndef SLICE_H
#define SLICE_H

#include "instance.h"
#include "local_search.h"
#include "utils.h"

/*
 *   low level operator
 *   slice will not violate load and time window constraints
 *   if the graph satisfy triangle inequality, then this operator will give a degenerate solution.
 *   Otherwise, this operator may improve the solution
 *
 *   This operator will only be used in the feasible search phase
 */


class Preslice : public MoveOperator
{
public:

    Preslice(){
        move_result = MoveResult(NeighborOperator::PRE_SLICE);
    };

    /*!
 * @details presert from Task u to Task i
 * @param ns
 * @param mcgrp
 * @param u
 * @param i
 * @param move_result
 * @return
 */
    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, const int b);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(LocalSearch &ns, const MCGRPTW &mcgrp,
                                                             const int b);

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool update_score(LocalSearch &ns) override;
};

//low level operator
class Postslice : public MoveOperator
{
public:

    Postslice(){
        move_result = MoveResult(NeighborOperator::POST_SLICE);
    };

    /*!
 * prosert from Task u to Task i
 * @param ns
 * @param mcgrp
 * @param u
 * @param i
 * @param move_result
 * @return
 */
    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, const int b);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(LocalSearch &ns, const MCGRPTW &mcgrp,
                                                             const int b);

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool update_score(LocalSearch &ns) override;
};


//high level operator
class Slice : public MoveOperator
{
    int pre_slice_times;
    int post_slice_times;
    Preslice preslice;
    Postslice postslice;
public:
    Slice(): preslice(Preslice()), postslice(Postslice())
    ,pre_slice_times(0)
    ,post_slice_times(0){
        move_result = MoveResult(NeighborOperator::SLICE);
    };

    /*----------------High speed neighbor search---------------------*/
    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    /*!
     * @details insert Task b to Task j
     * @param ns
     * @param mcgrp
     * @param j
     * @param b
     * @param M
     * @param policy
     * @return
     */
    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, const int b);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    bool update_score(LocalSearch &ns) override;
};


#endif //SLICE_H
