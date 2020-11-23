#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"


//low level operator
class Invert : public MoveOperator
{
public:
    Invert(){
        move_result = MoveResult(NeighborOperator::INVERT);
    };

    /***********************************************************/

    bool search(HighSpeedNeighBorSearch &ns, const class MCGRP &mcgrp, int chosen_task) override;

    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int u);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<RouteInfo::TimeTable> expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                                                     int u, int u_tilde, bool allow_infeasible);

    bool update_score(HighSpeedNeighBorSearch &ns) override;
};