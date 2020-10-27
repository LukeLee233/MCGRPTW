//
// Created by luke on 2020/1/30.
//


#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "MCGRP.h"

/*!
     * @details get hamming distance between two solution
     * @param mcgrp
     * @param a
     * @param b
     * @return
     */
double hamming_dist(const MCGRP &mcgrp, const vector<int> &a, const vector<int> &b);

double edit_distance(const vector<int>& neg_sol_a, const vector<int>& neg_sol_b);

#endif //SIMILARITY_H
