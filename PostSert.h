#pragma once

#include "utils.h"
#include "NeighborSearch.h"


//low level operator
class Postsert
{
public:
    MOVE move_result;
public:

    Postsert(): move_result(MOVE(NeighborOperator::POSTSERT)){};

    /*!
 * prosert from Task u to Task i
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
};