#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"


//low level operator
class Invert
{
public:
    MOVE move_result;

    Invert() : move_result(MOVE(NeighborOperator::INVERT)){};

    /***********************************************************/

    bool search(HighSpeedNeighBorSearch &ns, const class MCGRP &mcgrp, int chosen_task);

    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int u);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<RouteInfo::TimeTable> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                     int u, int u_tilde, bool allow_infeasible);


    void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};