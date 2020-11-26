//
// Created by luke on 2020/1/27.
//

#ifndef EXTRACTION_H
#define EXTRACTION_H

#include "MCGRP.h"
#include "NeighborSearch.h"
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
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int b);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<RouteInfo::TimeTable>>
    expected_time_table(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp,
                        const int b);

    bool update_score(HighSpeedNeighBorSearch &ns) override;
};

#endif //EXTRACTION_H
