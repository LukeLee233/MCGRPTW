#pragma once

#include "NeighborSearch.h"
#include "MCGRP.h"
#include "Flip.h"
#include "SwapEnds.h"
#include <vector>


class NewTwoOpt{
//    static const int length = 2;
    MCGRPMOVE move_result;
    int flip_times;
    int swapends_times;

    NewFlip flip;
    NewSwapEnds swap_ends;

    bool before(const int a,const int b, HighSpeedNeighBorSearch &ns);
public:

    NewTwoOpt():flip(NewFlip())
    ,swap_ends(NewSwapEnds())
    ,move_result(MCGRPMOVE(NeighborOperator::TWO_OPT))
    ,flip_times(0)
    ,swapends_times(0){};


    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task);

    bool considerable_move(HighSpeedNeighBorSearch &ns,
                                  const MCGRP &mcgrp,
                                  int a, int b, int c, int d);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

    void unit_test(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);
};


