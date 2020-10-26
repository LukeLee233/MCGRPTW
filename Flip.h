#pragma once
#include "utils.h"
#include "NeighborSearch.h"
#include "MCGRP.h"

//low level operator
class NewFlip
{
    //This will flip the whole tasks within two ends
public:
    MCGRPMOVE move_result;

    /*!
     * @details get the entire sequence needed to be fliped
     * @param ns
     * @param start
     * @param end
     * @return
     */
    vector<int> get_sequence(HighSpeedNeighBorSearch &ns, const int start, const int end);
public:
    NewFlip():move_result(MCGRPMOVE(NeighborOperator::FLIP)){};

    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int start_task, int end_task);

    vector<MCGRPRoute::Timetable>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        vector<int> &invert_seq,
                        bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};
