#ifndef _TWOOPT_H_
#define _TWOOPT_H_

#include "local_search.h"
#include "instance.h"
#include "flip.h"
#include "swapends.h"
#include <vector>


class TwoOpt : public MoveOperator{
    int flip_times;
    int swapends_times;

    Flip flip;
    SwapEnds swap_ends;

public:

    TwoOpt(): MoveOperator(),
              flip(Flip()),
              swap_ends(SwapEnds()){
        flip_times = 0;
        swapends_times = 0;
        move_result = MoveResult(NeighborOperator::TWO_OPT);
    };

    bool search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task) override;

    bool considerable_move(LocalSearch &ns,
                           const MCGRPTW &mcgrp,
                           int a, int b, int c, int d);

    void move(LocalSearch &ns, const MCGRPTW &mcgrp);


    bool update_score(LocalSearch &ns) override;
};


#endif