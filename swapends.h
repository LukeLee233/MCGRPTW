#ifndef _SWAPENDS_H_
#define _SWAPENDS_H_

#include "utils.h"
#include "instance.h"
#include "local_search.h"

struct RouteSegment
{
    // Contains information about a particular segment of a route.
    int segment_start;
    int segment_end;
    int num_custs;
    int load;
    double len;
};

class SwapEnds : public MoveOperator
{
public:
    SwapEnds():MoveOperator(){
        call_times = 0;
        success_times = 0;
        move_result = MoveResult(NeighborOperator::SWAP_ENDS);
    }
    
    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, const int chosen_task, const int neighbor_task, const int chosen_route, const int neighbor_route);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    vector<vector<RouteInfo::TimeTable>>
    expected_time_table(LocalSearch &ns,
                        const MCGRPTW &mcgrp,
                        int a,
                        int v,
                        int a_route,
                        int v_route,
                        const struct RouteSegment& a_seg,
                        const struct RouteSegment& v_seg,
                        bool allow_infeasible);

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool update_score(LocalSearch &ns) override;
};

RouteSegment get_segment_info(const MCGRPTW &mcgrp, LocalSearch &ns, const int chosen_task);

#endif