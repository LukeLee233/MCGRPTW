//
// Created by luke on 2020/11/5.
//

#ifndef OPERATOR_H
#define OPERATOR_H

#include "utils.h"

class MCGRP;
class HighSpeedNeighBorSearch;

class MoveOperator{
public:
    MoveResult move_result;

    virtual bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) = 0;
};

void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, MoveOperator& move_operator);


#endif //OPERATOR_H
