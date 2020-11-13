//
// Created by luke on 2020/1/30.
//


#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "MCGRP.h"

/*!
     * @details get hamming distance between two solution
     * @param mcgrp
     * @param neg_a
     * @param neg_b
     * @return
     */
double hamming_dist(const MCGRP &mcgrp, const vector<int> &neg_a, const vector<int> &neg_b);

double edit_distance(const vector<int>& neg_sol_a, const vector<int>& neg_sol_b);

#endif //SIMILARITY_H
