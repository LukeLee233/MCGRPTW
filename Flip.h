#pragma once
#include "utils.h"
#include "NeighborSearch.h"
#include "MCGRP.h"

//low level operator
class NewFlip : public MoveOperator
{
    //This will flip the whole tasks within two ends
public:
    /*!
     * @details get the entire sequence needed to be fliped
     * @param ns
     * @param start
     * @param end
     * @return
     */
    vector<int> get_sequence(HighSpeedNeighBorSearch &ns, const int start, const int end);
public:
    NewFlip(){
        move_result = MoveResult(NeighborOperator::FLIP);
    };

    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int start_task, int end_task);

    vector<RouteInfo::TimeTable>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        vector<int> &invert_seq,
                        bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;
};
