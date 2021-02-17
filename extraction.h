//
// Created by luke on 2020/1/27.
// doesn't support infeasible and descent search
//

#ifndef EXTRACTION_H
#define EXTRACTION_H

#include "instance.h"
#include "local_search.h"
#include "utils.h"

/*
 *   low level operator
 *   extraction will not violate load and time window constraints
 *   if the graph satisfy triangle inequality, then this operator will give a degenerate solution.
 *   Otherwise, this operator may improve the solution
 *
 *   This operator will only be used in the feasible search phase
 */

class Extraction : public MoveOperator
{
public:
    Extraction(){
        move_result = MoveResult(NeighborOperator::EXTRACTION);
    };

    /*----------------High speed neighbor search---------------------*/
    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

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
    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, const int b);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    vector<vector<RouteInfo::TimeTable>>
    expected_time_table(LocalSearch &ns, const MCGRPTW &mcgrp,
                        const int b);

    bool update_score(LocalSearch &ns) override;

};

#endif //EXTRACTION_H
