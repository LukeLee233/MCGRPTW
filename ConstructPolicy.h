//
// Created by luke on 2019/11/29.
//


#ifndef MCGRP_SEARCHPOLICY_H
#define MCGRP_SEARCHPOLICY_H

#include "MCGRP.h"
#include "utils.h"
#include "NeighborSearch.h"

/*!
 * nearest_scanning constructive algorithm
 * @param mcgrp: formatted instance
 * @param rs_indi: solution
 */
void nearest_scanning(const MCGRP &mcgrp, Individual &rs_indi);

vector<vector<int>>
nearest_scanning(const MCGRP &mcgrp, vector<int> unserved_task_set);

/*!
 * execute merge and split to the current solution
 * @param ns
 * @param mcgrp
 */
void merge_split(class NeighBorSearch &ns, const MCGRP &mcgrp, const int merge_size, const int pseudo_capacity);
void merge_split(class HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int merge_size, const int pseudo_capacity);

/*!
 * @details use different policies to merge task
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> split_task(Policy& policy,const MCGRP &mcgrp, const vector<int> &tasks, const int constraint);

/*!
 * @details use nearest L2 distance policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> nearest_growing(const MCGRP &mcgrp, vector<int> tasks,const int constraint);

/*!
 * @details use nearest L2 distance to depot policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> nearest_depot_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint);

/*!
 * @details use max yield  policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> maximum_yield_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint);

/*!
 * @details use min yield  policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> minimum_yield_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint);

/*!
 * @details use mixture policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> mixture_growing(const MCGRP &mcgrp, vector<int> tasks, const int constraint);

#endif //MCGRP_SEARCHPOLICY_H
