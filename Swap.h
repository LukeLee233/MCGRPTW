#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"

//low level operator
class NewSwap
{
public:
    MCGRPMOVE move_result;
public:

    NewSwap() : move_result(MCGRPMOVE(NeighborOperator::SWAP)){};

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task);

    /*!
     * @details swap task j and task b
     * @param ns
     * @param C
     * @param j
     * @param b
     * @param rules
     * @return
     */
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &C, int i,int u);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<MCGRPRoute::Timetable>> expected_time_table(HighSpeedNeighBorSearch &ns,const MCGRP &mcgrp,
                                                              int u,int u_tilde, int i, int i_tilde, bool allow_infeasible);

    void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};