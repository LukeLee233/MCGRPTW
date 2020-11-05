#pragma once

#include "MCGRP.h"
#include "NeighborSearch.h"
#include "utils.h"
#include "Presert.h"
#include "PostSert.h"

//high level operator
class SingleInsert : public MoveOperator
{
    int presert_times;
    int postsert_times;
    Presert presert;
    Postsert postsert;
public:
    SingleInsert()
    : presert(Presert()), postsert(Postsert()), presert_times(0), postsert_times(0){
        move_result=MoveResult(NeighborOperator::SINGLE_INSERT);
    };


    /*----------------High speed neighbor search---------------------*/
    bool search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task) override;

    /*!
     * @details insert Task b to Task j
     * @param ns
     * @param mcgrp
     * @param j
     * @param b
     * @param M
     * @param policy
     * @return
     */
    bool considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int b, const int j);

    void move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp);

};
