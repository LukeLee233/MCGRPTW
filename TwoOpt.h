#pragma once

#include "NeighborSearch.h"
#include "MCGRP.h"
#include "flip.h"
#include "swapends.h"
#include <vector>


class NewTwoOpt : public MoveOperator{
    int flip_times;
    int swapends_times;

    NewFlip flip;
    NewSwapEnds swap_ends;

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

    bool update_score(HighSpeedNeighBorSearch &ns) override;
};


