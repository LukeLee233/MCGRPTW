#pragma once

#include "utils.h"
#include "NeighborSearch.h"


//low level operator
class Postsert
{
public:
    MCGRPMOVE move_result;
public:

    Postsert(): move_result(MCGRPMOVE(NeighborOperator::POSTSERT)){};

    /*!
 * prosert from task u to task i
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
    vector<vector<MCGRPRoute::Timetable>> expected_time_table(HighSpeedNeighBorSearch &ns,const MCGRP &mcgrp,
                                                      int u,int u_tilde, const int i,bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};