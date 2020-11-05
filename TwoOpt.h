#pragma once

#include "NeighborSearch.h"
#include "MCGRP.h"
#include "Flip.h"
#include "SwapEnds.h"
#include <vector>


class NewTwoOpt : public MoveOperator{
//    static const int length = 2;
    int flip_times;
    int swapends_times;

    NewFlip flip;
    NewSwapEnds swap_ends;

    bool before(const int a,const int b, HighSpeedNeighBorSearch &ns);
public:

    NewTwoOpt():flip(NewFlip())
    ,swap_ends(NewSwapEnds())
    ,flip_times(0)
    ,swapends_times(0){
        move_result = MoveResult(NeighborOperator::TWO_OPT);
    };


    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    bool considerable_move(HighSpeedNeighBorSearch &ns,
                                  const MCGRP &mcgrp,
                                  int a, int b, int c, int d);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

};


