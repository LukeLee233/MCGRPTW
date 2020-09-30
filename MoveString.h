#pragma once
#include "utils.h"
#include "MCGRP.h"
#include "NeighborSearch.h"

//low level operator
class MoveString
{
public:
    MCGRPMOVE move_result;
public:

    MoveString()
        : move_result(MCGRPMOVE(NeighborOperator::MOVE_STRING))
    {};

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
    bool considerable_move(NeighBorSearch &ns, const MCGRP &mcgrp, vector<int> disturbance_seq, const int u);

    void move(NeighBorSearch &ns, const MCGRP &mcgrp);
};

//low level operator
class PreMoveString
{
public:
    MCGRPMOVE move_result;
public:

    PreMoveString()
        : move_result(MCGRPMOVE(NeighborOperator::PRE_MOVE_STRING))
    {};

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
    vector<vector<MCGRPRoute::Timetable>>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        vector<int> &disturbance_seq,
                        const int i, bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};

//low level operator
class PostMoveString
{
public:
    MCGRPMOVE move_result;
public:
    PostMoveString()
        : move_result(MCGRPMOVE(NeighborOperator::POST_MOVE_STRING))
    {};

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
    vector<vector<MCGRPRoute::Timetable>>
    expected_time_table(HighSpeedNeighBorSearch &ns,
                        const MCGRP &mcgrp,
                        vector<int> &disturbance_seq,
                        const int i, bool allow_infeasible);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};