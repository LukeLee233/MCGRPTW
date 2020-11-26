#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"

struct RouteSegment
{
    // Contains information about a particular segment of a route.
    int segment_start;
    int segment_end;
    int num_custs;
    int load;
    double len;
};

class NewSwapEnds : public MoveOperator
{
public:
    NewSwapEnds(){
        move_result = MoveResult(NeighborOperator::SWAP_ENDS);
    }
    
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

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    bool update_score(HighSpeedNeighBorSearch &ns) override;
};

RouteSegment get_segment_info(const MCGRP &mcgrp, HighSpeedNeighBorSearch &ns, const int chosen_task);
