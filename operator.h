//
// Created by luke on 2020/11/5.
//

#ifndef OPERATOR_H
#define OPERATOR_H

#include "utils.h"

class MCGRP;
class HighSpeedNeighBorSearch;
class Policy;

namespace viterbi{

struct BestDecode{
    int cost;
    vector<int> seq;
    vector<int> arrive_time;
};

struct PseudoTask: public Task{
    PseudoTask(int head_node,
               int tail_node,
               int demand,
               int serve_cost,
               int trave_time,
               Window& time_window): Task(){
        task_name = "pseudo";
        task_id = INT32_MAX;

        this->head_node = head_node;
        this->tail_node = tail_node;
        this->demand = demand;
        this->serv_cost = serve_cost;
        this->trave_time = trave_time;
    };
};


BestDecode viterbi_decode(const PseudoTask* first,const PseudoTask* last,
                          const vector<int>& move_seq, const MCGRP& mcgrp,Policy& policy);

}

class MoveOperator{
public:
    MoveResult move_result;

    int hit_count;
    int attempt_count;

    MoveOperator();

    virtual bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) = 0;
    virtual bool update_score(HighSpeedNeighBorSearch &ns) = 0;


    viterbi::PseudoTask Seq2Pseudo(const vector<int>& seq, int phase);
    vector<int> getHitInfo();
    vector<vector<int>> generate_possible_seq(const MCGRP &mcgrp,const vector<int>& seq);

};

void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, MoveOperator& move_operator);


#endif //OPERATOR_H
