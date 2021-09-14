#pragma once
#include "utils.h"
#include "instance.h"
#include "local_search.h"

//low level operator
class Swap : public MoveOperator
{
public:

    Swap() : MoveOperator(){
        move_result = MoveResult(NeighborOperator::SWAP);
    };

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    /*!
     * @details swap Task j and Task b
     * @param ns
     * @param C
     * @param j
     * @param b
     * @param rules
     * @return
     */
    bool considerable_move(LocalSearch &ns, const MCGRPTW &C, int i, int u);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    vector<vector<RouteInfo::TimeTable>> expected_time_table(LocalSearch &ns, const MCGRPTW &mcgrp,
                                                             int u, int u_tilde, int i, int i_tilde, bool allow_infeasible);

    void _apply_reward(LocalSearch &ns, double reward) override;
};