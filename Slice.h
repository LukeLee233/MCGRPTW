//
// Created by luke on 2020/1/26.
//

#ifndef SLICE_H
#define SLICE_H

#include "MCGRP.h"
#include "NeighborSearch.h"
#include "utils.h"

/*
 *   low level operator
 *   slice will not violate load and time window constraints
 *   if the graph satisfy triangle inequality, then this operator will give a degenerate solution.
 *   Otherwise, this operator may improve the solution
 *
 *   This operator will only be used in the feasible search phase
 */


class Preslice
{
public:
    MOVE move_result;
public:

    Preslice(): move_result(MOVE(NeighborOperator::PRE_SLICE)){};

    /*!
 * @details presert from Task u to Task i
 * @param ns
 * @param mcgrp
 * @param u
 * @param i
 * @param move_result
 * @return
 */
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                             const int b);
};

//low level operator
class Postslice
{
public:
    MOVE move_result;
public:

    Postslice(): move_result(MOVE(NeighborOperator::POST_SLICE)){};

    /*!
 * prosert from Task u to Task i
 * @param ns
 * @param mcgrp
 * @param u
 * @param i
 * @param move_result
 * @return
 */
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                             const int b);

};


//high level operator
class Slice
{
    MOVE move_result;
    int pre_slice_times;
    int post_slice_times;
    Preslice preslice;
    Postslice postslice;
public:
    Slice(): preslice(Preslice()), postslice(Postslice())
    ,move_result(MOVE(NeighborOperator::SLICE))
    ,pre_slice_times(0)
    ,post_slice_times(0){};

    /*----------------High speed neighbor search---------------------*/
    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task);

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
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,const int b);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};


#endif //SLICE_H
