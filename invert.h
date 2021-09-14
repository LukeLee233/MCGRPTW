#pragma once

#include "utils.h"
#include "instance.h"
#include "local_search.h"


//low level operator
class Invert : public MoveOperator
{
public:
    Invert() : MoveOperator() {
        move_result = MoveResult(NeighborOperator::INVERT);
    };

    /***********************************************************/

    bool search(LocalSearch &ns, const class MCGRPTW &mcgrp, int chosen_task) override;

    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, int u);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    vector<RouteInfo::TimeTable> expected_time_table(LocalSearch &ns, const MCGRPTW &mcgrp,
                                                     int u, int u_tilde, bool allow_infeasible);

    void _apply_reward(LocalSearch &ns, double reward) override;
};