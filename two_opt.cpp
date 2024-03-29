#include "two_opt.h"

using namespace std;

bool TwoOpt::search(LocalSearch &ns, const MCGRPTW &mcgrp, int chosen_task){

    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"error, Wrong Task");

    MoveResult BestM;

    int offset = 0;
    if(ns.policy.has_rule(TOLERANCE)){
        offset = mcgrp._rng.Randint(0,mcgrp.actual_task_num);
    }
    ns.create_search_neighborhood(mcgrp, {chosen_task},"basic", offset);

    int b = chosen_task;
    int a = ns.solution[b]->pre->ID;
    int c = ns.solution[b]->next->ID;

    for(auto neighbor_task : ns.search_space){
        My_Assert(neighbor_task != b,"Neighbor Task can't be itself at all!");

        if(neighbor_task != DUMMY){

            //j can't be dummy and b can't be dummy neither here
            int j = neighbor_task;
            My_Assert(j >= 1 && j <= mcgrp.actual_task_num,"Wrong Task");

            int i = ns.solution[j]->pre->ID;
            int k = ns.solution[j]->next->ID;

            //total having four cases:
            //(a,b)<->(i,j)
            //(a,b)<->(j,k)
            //(b,c)<->(i,j)
            //(b,c)<->(j,k)
            if (considerable_move(ns,mcgrp, a, b, i, j)){
                if(ns.policy.has_rule(FIRST_ACCEPT)){
                    move(ns,mcgrp);
                    return true;
                }
                else if (ns.policy.has_rule(BEST_ACCEPT)){
                    if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                        BestM = move_result;
                }
                else{
                    My_Assert(false,"Unknown accept rule!");
                }
            }

            if (considerable_move(ns,mcgrp, a, b, j, k)){
                if(ns.policy.has_rule(FIRST_ACCEPT)){
                    move(ns,mcgrp);
                    return true;
                }
                else if (ns.policy.has_rule(BEST_ACCEPT)){
                    if (ns.policy.check_result(mcgrp,ns, move_result,BestM))
                        BestM = move_result;
                }
                else{
                    My_Assert(false,"Unknown accept rule!");
                }
            }

            if (considerable_move(ns, mcgrp, b, c, i, j)){
                if(ns.policy.has_rule(FIRST_ACCEPT)){
                    move(ns,mcgrp);
                    return true;
                }
                else if (ns.policy.has_rule(BEST_ACCEPT)){
                    if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                        BestM = move_result;
                }
                else{
                    My_Assert(false,"Unknown accept rule!");
                }
            }

            if (considerable_move(ns,mcgrp, b, c, j, k)){
                if(ns.policy.has_rule(FIRST_ACCEPT)){
                    move(ns,mcgrp);
                    return true;
                }
                else if (ns.policy.has_rule(BEST_ACCEPT)){
                    if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                        BestM = move_result;
                }
                else{
                    My_Assert(false,"Unknown accept rule!");
                }
            }
        }
        else{
            DEBUG_PRINT("Neighbor Task is dummy Task");
            //j is dummy here and b can't be dummy neither
            //each start and end location of each route will be considered
            //total 4 x route_nums cases
            int current_start = ns.solution.very_start->next->ID;

            while (ns.solution[current_start]->next != ns.solution.very_end){


                // Consider the start location
                int j = current_start;
                My_Assert(ns.solution[j]->pre->ID < 0,"Wrong Task");

                //case: (a,b)<->(dummy,j)
                if (considerable_move(ns, mcgrp, a, b, ns.solution[j]->pre->ID, j)){
                    if(ns.policy.has_rule(FIRST_ACCEPT)){
                        move(ns,mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)){
                        if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                            BestM = move_result;
                    }
                    else{
                        My_Assert(false,"Unknown accept rule!");
                    }
                }

                //case: (b,c)<->(dummy,j)
                if (considerable_move(ns, mcgrp, b, c, ns.solution[j]->pre->ID, j)){
                    if(ns.policy.has_rule(FIRST_ACCEPT)){
                        move(ns,mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)){
                        if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                            BestM = move_result;
                    }
                    else{
                        My_Assert(false,"Unknown accept rule!");
                    }
                }

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.solution[ns.routes[current_route]->end]->pre->ID;
                j = current_end;

                My_Assert(ns.solution[j]->next->ID < 0,"Wrong Task");

                //case: (a,b)<->(j,dummy)
                if (considerable_move(ns, mcgrp, a, b, j, ns.solution[j]->next->ID)){
                    if(ns.policy.has_rule(FIRST_ACCEPT)){
                        move(ns,mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)){
                        if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                            BestM = move_result;
                    }
                    else{
                        My_Assert(false,"Unknown accept rule!");
                    }
                }

                //case: (b,c)<->(j,dummy)
                if (considerable_move(ns, mcgrp, b, c, j, ns.solution[j]->next->ID)){
                    if(ns.policy.has_rule(FIRST_ACCEPT)){
                        move(ns,mcgrp);
                        return true;
                    }
                    else if (ns.policy.has_rule(BEST_ACCEPT)){
                        if (ns.policy.check_result(mcgrp,ns,move_result,BestM))
                            BestM = move_result;
                    }
                    else{
                        My_Assert(false,"Unknown accept rule!");
                    }
                }

                // Advance to next route's starting Task
                if(ns.solution[current_end]->next == ns.solution.very_end){
                    current_start = current_end;
                }else{
                    current_start = ns.solution[current_end]->next->next->ID;
                }
            }
        }
    }

    if (ns.policy.has_rule(FIRST_ACCEPT))
    {
        DEBUG_PRINT("No actual move: First Accept Rule");
        return false;
    }
    else if (ns.policy.has_rule(BEST_ACCEPT)){
        if(BestM.considerable == false){
            DEBUG_PRINT("No actual move: Best Accept Rule");
            return false;
        }
        else{
            move_result = BestM;
            move(ns,mcgrp);
            return true;
        }
    }
    else{
        My_Assert(false,"Unknown accept rule");
    }

    DEBUG_PRINT("2-opt Flip times:" + to_string(flip_times));
    DEBUG_PRINT("2-opt SwapEnds times:" + to_string(swapends_times));
}


bool
TwoOpt::considerable_move(
    LocalSearch &ns,
    const MCGRPTW &mcgrp,
    int a,
    int b,
    int c,
    int d){

    call_times++;
    // A easy prediction for overlapping
    if ((a == c) || (a == d) || (b == c) || (b == d)) {
        return false;
    }


    int a_route, c_route;

    if (a > 0)
        a_route = ns.solution[a]->route_id;
    else{
        My_Assert(b > 0,"Task A is dummy and Task B is dummy too? Come on! ");
        a_route = ns.solution[b]->route_id;
    }


    if (c > 0)
        c_route = ns.solution[c]->route_id;
    else
    {
        My_Assert(d > 0,"Task C is dummy and Task D is dummy too? Come on! ");
        c_route = ns.solution[d]->route_id;
    }


    if (a_route == c_route) {
        // This is flip operation
        // same route: the 2opt move here corresponds to reversing the sequence of
        // the sub route that is between (a|b) and (c|d), open interval
        DEBUG_PRINT("Flip operator in 2-opt");
        flip_times++;

        if(ns.before(a, c)){      //...ab...cd...
            if(flip.considerable_move(ns, mcgrp, a, d) && ns.policy.check_move(mcgrp,ns,flip.move_result)){
                move_result = flip.move_result;
                return true;
            }
            else{
                move_result.reset();
                move_result.move_type = NeighborOperator::TWO_OPT;
                return false;
            }
        }
        else{        //...cd...ab...
            if(flip.considerable_move(ns,mcgrp,c,b) && ns.policy.check_move(mcgrp,ns,flip.move_result)){
                move_result = flip.move_result;
                return true;
            }
            else{
                move_result.reset();
                move_result.move_type = NeighborOperator::TWO_OPT;
                return false;
            }
        }
    }
    else // This is SwapEnds operation
    {
        // Cross routes, feasibility need checked
        // different routes: the 2-opt move here corresponds to swap the sequence of
        // the sub-route that is after a and c, open interval

        // Example: ( a & v input): DEPOT-i-a-b-j-k-l-DEPOT and DEPOT-t-u-v-w-x-y-z-DEPOT becomes
        // DEPOT-i-a-w-x-y-z-DEPOT and DEPOT-t-u-v-b-j-k-l-DEPOT
        DEBUG_PRINT("SwapEnds operator in 2-opt");
        swapends_times++;

        if(swap_ends.considerable_move(ns, mcgrp, a, c, a_route, c_route)
        && ns.policy.check_move(mcgrp,ns,swap_ends.move_result)){
            move_result = swap_ends.move_result;
            return true;
        }
        else{
            move_result.reset();
            move_result.move_type = NeighborOperator::TWO_OPT;
            return false;
        }

    }

    My_Assert(false,"Can't reach here!");
}


void TwoOpt::move(LocalSearch &ns, const MCGRPTW &mcgrp){

    DEBUG_PRINT("execute a 2-opt move");

    My_Assert(move_result.considerable,"Invalid predictions");

    if(move_result.move_type == NeighborOperator::FLIP){
        flip.move_result = move_result;
        flip.move(ns,mcgrp);
    }
    else if (move_result.move_type == NeighborOperator::SWAP_ENDS){
        swap_ends.move_result = move_result;
        swap_ends.move(ns,mcgrp);
    }
    else{
        My_Assert(false,"error, unknown operator");
    }


    move_result.reset();
    move_result.move_type = NeighborOperator::TWO_OPT;
}

