#ifndef _FLIP_H_
#define _FLIP_H_

#include "utils.h"
#include "local_search.h"
#include "instance.h"

//low level operator
class Flip : public MoveOperator
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
    vector<int> get_sequence(LocalSearch &ns, const int start, const int end);

public:
    Flip():MoveOperator(){
        move_result = MoveResult(NeighborOperator::FLIP);
    };

    bool considerable_move(LocalSearch &ns, const MCGRPTW &mcgrp, int start_task, int end_task);

    vector<RouteInfo::TimeTable>
    expected_time_table(LocalSearch &ns,
                        const MCGRPTW &mcgrp,
                        vector<int> &invert_seq,
                        bool allow_infeasible);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool update_score(LocalSearch &ns) override;
};

#endif