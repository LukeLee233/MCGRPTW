//
// Created by luke on 2019/11/29.
//


#ifndef MCGRP_SEARCHPOLICY_H
#define MCGRP_SEARCHPOLICY_H

#include "MCGRP.h"
#include "utils.h"
#include "NeighborSearch.h"

// different Distance metric between tasks
class Distance{
protected:
    string name;

    vector<vector<double>> distance_matrix;
public:
    Distance(const MCGRP& mcgrp, const string &name);

    virtual double operator()(const MCGRP& mcgrp, const int task_a, const int task_b) = 0;
};

class CostDistance: public Distance{

public:
    CostDistance(const MCGRP& mcgrp);

    double operator()(const MCGRP& mcgrp, const int task_a, const int task_b) override;
};


class PathConstructor{
protected:
    string name;

public:
    PathConstructor(const MCGRP& mcgrp,const string& name_);

    // mode : [allow_infeasible, feasible]
    virtual Individual
    operator()(const MCGRP& ns,const vector<int> &taskList = vector<int>(), const string& mode="feasible") = 0;
};

class NearestScanner: public PathConstructor{
private:
    Distance& distance;
public:
    NearestScanner(const MCGRP &mcgrp, Distance& distance_);

    Individual
    operator()(const MCGRP& mcgrp,const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};


/*!
 *
 * @param mcgrp
 * @param unserved_task_set if not empty, means execute a partial construction
 * @return
 */
Individual nearest_scanning(const MCGRP &mcgrp, vector<int> unserved_task_set);

/*!
 * Random Task Flower construction policy
 * @param mcgrp
 * @param task_list
 * @param giant if true, a giant tour is going to be built
 * @return
 */
Individual RTF(const MCGRP &mcgrp, const vector<int>& task_list, bool giant);

/*!
 * execute merge and split to the current solution
 * @param ns
 * @param mcgrp
 */
void merge_split(class HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, const int merge_size, const int pseudo_capacity);

/*!
 * @details use different policies to merge Task
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> split_task(Policy& policy,const MCGRP &mcgrp, const vector<int> &tasks, const int load_constraint);

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

vector<vector<int>> tour_splitting(const MCGRP &mcgrp, vector<int>& task_list);


struct BestServeSeq{
    vector<int> sequence;
    vector<int> arrive_time_tbl;
    int cost = INT32_MAX;
};

BestServeSeq viterbi_decoding(const MCGRP &mcgrp, const vector<int>& task_list, bool allow_infeasible);

#endif //MCGRP_SEARCHPOLICY_H
