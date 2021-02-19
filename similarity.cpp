//
// Created by luke on 2020/1/30.
//

#include "similarity.h"
#include <algorithm>

double hamming_dist(const MCGRPTW &mcgrp, const vector<int> &neg_a, const vector<int> &neg_b)
{
    vector<vector<int>> atom_link;
    //include offset 1
    atom_link.resize(mcgrp.node_num + 1);

    My_Assert(atom_link[0].empty(), "Incorrect initialization!");

    int pre_task, cur_task, next_task;
    int source_node, target_node;

    //parse solution neg_a
    for (int i = 0; i < neg_a.size(); ++i) {
        cur_task = abs(neg_a[i]);

        //Special case: start of the route
        if (neg_a[i] < 0) {
            //Inside & outside will be considered meantime
            pre_task = DUMMY;
            source_node = mcgrp.inst_tasks[pre_task].head_node;
            target_node = mcgrp.inst_tasks[cur_task].tail_node;
            atom_link[source_node].push_back(target_node);

            source_node = mcgrp.inst_tasks[pre_task].tail_node;
            target_node = mcgrp.inst_tasks[cur_task].head_node;
            atom_link[source_node].push_back(target_node);

        }

        if (i == neg_a.size() - 1)
            next_task = DUMMY;
        else
            next_task = max(DUMMY, neg_a[i + 1]);

        //Inside & outside will be considered meantime
        source_node = mcgrp.inst_tasks[cur_task].head_node;
        target_node = mcgrp.inst_tasks[next_task].tail_node;
        atom_link[source_node].push_back(target_node);

        source_node = mcgrp.inst_tasks[cur_task].tail_node;
        target_node = mcgrp.inst_tasks[next_task].head_node;
        atom_link[source_node].push_back(target_node);
    }


    int same_count = 0, num_rts = 0;


    //parse route neg_b
    for (int i = 0; i < neg_b.size(); ++i) {
        cur_task = abs(neg_b[i]);

        //Special case: start of the route
        if (neg_b[i] < 0) {
            num_rts++;
            pre_task = DUMMY;
            source_node = mcgrp.inst_tasks[pre_task].head_node;
            target_node = mcgrp.inst_tasks[cur_task].tail_node;

            auto ite = find(atom_link[source_node].begin(), atom_link[source_node].end(), target_node);
            if (ite != atom_link[source_node].end()) {
                atom_link[source_node].erase(ite);
                same_count++;
            }

            source_node = mcgrp.inst_tasks[pre_task].tail_node;
            target_node = mcgrp.inst_tasks[cur_task].head_node;

            ite = find(atom_link[source_node].begin(), atom_link[source_node].end(), target_node);
            if (ite != atom_link[source_node].end()) {
                atom_link[source_node].erase(ite);
                same_count++;
            }
        }


        if (i == neg_b.size() - 1)
            next_task = DUMMY;
        else
            next_task = max(DUMMY, neg_b[i + 1]);

        source_node = mcgrp.inst_tasks[cur_task].head_node;
        target_node = mcgrp.inst_tasks[next_task].tail_node;

        auto ite = find(atom_link[source_node].begin(), atom_link[source_node].end(), target_node);
        if (ite != atom_link[source_node].end()) {
            atom_link[source_node].erase(ite);
            same_count++;
        }

        source_node = mcgrp.inst_tasks[cur_task].tail_node;
        target_node = mcgrp.inst_tasks[next_task].head_node;

        ite = find(atom_link[source_node].begin(), atom_link[source_node].end(), target_node);
        if (ite != atom_link[source_node].end()) {
            atom_link[source_node].erase(ite);
            same_count++;
        }
    }

    double dist = 1.0 - double(same_count) / double(2 * (mcgrp.req_node_num + mcgrp.req_arc_num + mcgrp.req_edge_num + num_rts));

    return dist;
}

double edit_distance(const vector<int>& neg_sol_a, const vector<int>& neg_sol_b) {
    vector<vector<int>> shortestDist(neg_sol_a.size() + 1, vector<int>(neg_sol_b.size() + 1, 0));

    shortestDist[0][0] = 0;

    for(int i = 1;i<=neg_sol_a.size();i++){
        shortestDist[i][0] = i;
    }

    for(int j = 1; j<= neg_sol_b.size();j++){
        shortestDist[0][j] = j;
    }

    for(int i = 1;i<= neg_sol_a.size();i++){
        for(int j = 1; j <= neg_sol_b.size();j++){
            if(neg_sol_a[i - 1] == neg_sol_b[j - 1]){
                shortestDist[i][j] = shortestDist[i-1][j-1];
            }else{
                shortestDist[i][j] = 1 + min({shortestDist[i][j-1], shortestDist[i-1][j], shortestDist[i-1][j-1]});
            }
        }
    }

    return double(shortestDist.back().back()) / double(max(neg_sol_a.size(),neg_sol_b.size()));
}