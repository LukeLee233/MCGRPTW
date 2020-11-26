//
// Created by luke on 2020/11/25.
//

#ifndef ATTRACTION_H
#define ATTRACTION_H

#include "operator.h"
#include "slice.h"
#include "insert.h"

class Attraction : public MoveOperator
{
private:
    int sample_times;
    Postslice post_slice;
    XPostInsert post_insert;

    inline bool prob_valid(const vector<double> probs){
        return !all_of(probs.begin(),probs.end(),[](double i){ return i == 0;});
    }

public:
    Attraction(int sampleTimes = 50);

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
                      const int a,
                      const int u);


    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    bool update_score(HighSpeedNeighBorSearch &ns) override{
        // This action will be executed by post slice and post insert operator
        return true;
    };
};


#endif //ATTRACTION_H
