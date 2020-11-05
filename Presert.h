#pragma once

#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"

//low level operator
class Presert : public MoveOperator
{
public:
    MoveResult move_result;
public:

    Presert():move_result(MoveResult(NeighborOperator::PRESERT)){};


    /*!
 * @details presert from Task u to Task i
 * @param ns
 * @param mcgrp
 * @param u
 * @param i
 * @param move_result
 * @return
 */
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int u, const int i);

    /*!
     *
     * @param ns
     * @param mcgrp
     * @param u
     * @param u_tilde: if u is not equal u_tilde, means inverse case are considered
     * @param i
     * @param allow_infeasible
     * @return
     */
    vector<vector<RouteInfo::TimeTable>> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                             int u, int u_tilde, const int i, bool allow_infeasible);


    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;
};