//
// Created by luke on 2020/11/5.
//

#ifndef _INSERT_H_
#define _INSERT_H_

#include "operator.h"
#include "local_search.h"

class XPreInsert : public MoveOperator
{
public:
    const int length;

    XPreInsert(int length): MoveOperator(), length(length){
        move_result = MoveResult(NeighborOperator::PRE_INSERT);
    };

    /*!
     * @details Evaluates the move of inserting the string i-j-k(disturbance seq) between u and v (i.e.h-i-j-k-l & t-u-v-w )
	 * @details yielding h-l & t-u-i-j-k-v-w
     * @param ns
     * @param mcgrp
     * @param move_seq
     * @param u
     * @param move_result
     * @return
     */
    bool
    considerable_move(LocalSearch &ns,
                      const MCGRPTW &mcgrp,
                      const vector<int>& move_seq,
                      const int u);


    pair<double,double> cost_delta(LocalSearch &ns,
                                   const MCGRPTW &mcgrp,
                                   const vector<int>& move_seq,
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
    expected_time_table(LocalSearch &ns,
                        const MCGRPTW &mcgrp,
                        const vector<int> &move_sequence,
                        const int i, bool allow_infeasible);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool update_score(LocalSearch &ns) override;

};

class XPostInsert : public MoveOperator
{
public:
    const int length;

    XPostInsert(int length): length(length){
        move_result = MoveResult(NeighborOperator::POST_INSERT);
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
    considerable_move(LocalSearch &ns,
                      const MCGRPTW &mcgrp,
                      vector<int> move_seq,
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
    expected_time_table(LocalSearch &ns,
                        const MCGRPTW &mcgrp,
                        vector<int> &disturbance_seq,
                        const int i, bool allow_infeasible);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool update_score(LocalSearch &ns) override;


    pair<double,double> cost_delta(LocalSearch &ns,
                                   const MCGRPTW &mcgrp,
                                   const vector<int>& move_seq,
                                   const int u);
};


#endif //INSERT_H
