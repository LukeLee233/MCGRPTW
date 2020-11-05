#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"

//low level operator
class PreMoveString : public MoveOperator
{
public:

    PreMoveString(){
        move_result = MoveResult(NeighborOperator::PRE_MOVE_STRING);
    };

    /*!
     * @details Evaluates the move of inserting the string i-j-k(disturbance seq) between u and v (i.e.h-i-j-k-l & t-u-v-w )
	 * @details yielding h-l & t-u-i-j-k-v-w
     * @param ns
     * @param mcgrp
     * @param disturbance_seq
     * @param u
     * @param move_result
     * @return
     */
    bool
    considerable_move(HighSpeedNeighBorSearch &ns,
                           const MCGRP &mcgrp,
                           vector<int> disturbance_seq,
                           const int u);

    /*!
 *
 * @param ns
 * @param mcgrp
 * @param u
 * @param u_tilde: if u is not equal u_tilde, means inverse case are considered
 * @param i
 * @param allow_infeasible
 * @return
 */
    vector<vector<RouteInfo::TimeTable>>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        vector<int> &disturbance_seq,
                        const int i, bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;
};

//low level operator
class PostMoveString : public MoveOperator
{
public:
    PostMoveString(){
        move_result = MoveResult(NeighborOperator::POST_MOVE_STRING);
    };

    /*!
     * @details Evaluates the move of inserting the string i-j-k(disturbance seq) between u and v (i.e.h-i-j-k-l & t-u-v-w )
	 * @details yielding h-l & t-u-i-j-k-v-w
     * @param ns
     * @param mcgrp
     * @param disturbance_seq
     * @param u
     * @param move_result
     * @return
     */
    bool
    considerable_move(HighSpeedNeighBorSearch &ns,
                      const MCGRP &mcgrp,
                      vector<int> disturbance_seq,
                      const int u);

    /*!
 *
 * @param ns
 * @param mcgrp
 * @param u
 * @param u_tilde: if u is not equal u_tilde, means inverse case are considered
 * @param i
 * @param allow_infeasible
 * @return
 */
    vector<vector<RouteInfo::TimeTable>>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        vector<int> &disturbance_seq,
                        const int i, bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;
};