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
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                             const int b);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    bool update_score(HighSpeedNeighBorSearch &ns) override;
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
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                             const int b);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    bool update_score(HighSpeedNeighBorSearch &ns) override;
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
    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

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

    bool update_score(HighSpeedNeighBorSearch &ns) override;
};


#endif //SLICE_H
