//
// Created by luke on 2019/11/29.
//


#ifndef MCGRP_SEARCHPOLICY_H
#define MCGRP_SEARCHPOLICY_H

#include "instance.h"
#include "utils.h"


class LocalSearch;
class Policy;

class PathConstructor{
protected:
    string name;
    const MCGRPTW& mcgrp;

public:
    PathConstructor(const MCGRPTW& mcgrp_, const string& name_);

    // mode : [allow_infeasible, feasible]
    virtual Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") = 0;
};

class NearestScanner: public PathConstructor{
private:
    Distance& distance;
public:
    NearestScanner(const MCGRPTW &mcgrp, Distance& distance_);

    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};

class RTFScanner: public PathConstructor{
private:
    Distance& distance;
public:
    RTFScanner(const MCGRPTW &mcgrp, Distance &distance_);

    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};

class SampleScanner: public PathConstructor{
private:
    LearningDistance& distance;
    int sample_times;
public:
    SampleScanner(const MCGRPTW &mcgrp, LearningDistance &distance_, int sample_times_ = 50);

    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};


/*!
 * execute merge and split to the current solution
 * @param ns
 * @param mcgrp
 */
void merge_split(LocalSearch &ns, const MCGRPTW &mcgrp,
                 const int merge_size);

/*!
 * @details use different policies to merge Task
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> split_task(const MCGRPTW &mcgrp, LocalSearch& ns, const vector<int> &tasks, Policy& policy);

/*!
 * @details use nearest L2 distance policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> nearest_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use nearest L2 distance to depot policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> nearest_depot_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use max yield  policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> maximum_yield_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use min yield  policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> minimum_yield_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use mixture policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> mixture_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

vector<vector<int>> tour_splitting(const MCGRPTW &mcgrp, vector<int>& task_list);


struct BestServeSeq{
    vector<int> sequence;
    vector<int> arrive_time_tbl;
    int cost = INT32_MAX;
};

BestServeSeq viterbi_decoding(const MCGRPTW &mcgrp, const vector<int>& task_list, bool allow_infeasible);

#endif //MCGRP_SEARCHPOLICY_H
