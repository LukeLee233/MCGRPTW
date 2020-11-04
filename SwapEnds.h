#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"


class NewSwapEnds
{
public:
    MOVE move_result;
public:
    NewSwapEnds():move_result(MOVE(NeighborOperator::SWAP_ENDS)){}
    
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int chosen_task, const int neighbor_task,const int chosen_route,const int neighbor_route);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    vector<vector<RouteInfo::TimeTable>>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        int a,
                        int v,
                        int a_route,
                        int v_route,
                        const struct RouteSegment& a_seg,
                        const struct RouteSegment& v_seg,
                        bool allow_infeasible);

};