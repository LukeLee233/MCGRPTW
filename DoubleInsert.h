#pragma once
#include "utils.h"
#include "NeighborSearch.h"
#include "MoveString.h"


//high level operator
class DoubleInsert : public MoveOperator
{
    // Here you can set a given length to achieve m-sert
    const int length = 2;

    int presert_times;
    int postsert_times;

    PreMoveString pre_move_string;
    PostMoveString post_move_string;

public:
    DoubleInsert(): pre_move_string(PreMoveString())
    ,post_move_string(PostMoveString())
    ,presert_times(0),postsert_times(0){
        move_result = MoveResult(NeighborOperator::DOUBLE_INSERT);
    };


    /*!
 * @details get the successor of the chosen Task within the same route,exclude dummy Task
 * @param ns
 * @param chosen_task
 * @return
 */
    vector<int> get_successor_tasks(HighSpeedNeighBorSearch &ns, const int chosen_task);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    /*!
     * @details double sert from disturbance_seq to candidate_task
     * @param ns
     * @param mcgrp
     * @param disturbance_seq
     * @param j
     * @param M
     * @param policy
     * @return
     */
    bool considerable_move(HighSpeedNeighBorSearch &ns,
                                  const MCGRP &mcgrp,
                                  vector<int> disturbance_seq,
                                  const int j);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

};