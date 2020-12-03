#include "swap.h"
#include <iostream>
#include <vector>


using namespace std;

/*
 * High speed
 */

bool NewSwap::search(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int chosen_task)
{
    My_Assert(chosen_task >= 1 && chosen_task <= mcgrp.actual_task_num,"Wrong Task");

#ifdef DEBUG
    attempt_count = 0;
    hit_count = 0;
#endif

    MoveResult BestM;

    int offset = 0;
    if(ns.policy.has_rule(TOLERANCE)){
        offset = rand() % mcgrp.actual_task_num;
    }
    ns.create_search_neighborhood(mcgrp, {chosen_task},"basic",offset);

    int b = chosen_task;


    for(auto neighbor_task : ns.search_space) {
        My_Assert(neighbor_task != b, "neighbor Task can't be itself!");

        if (neighbor_task != DUMMY) {
            //j can't be dummy and b can't be dummy neither here

#ifdef DEBUG
            attempt_count++;
#endif

            int j = neighbor_task;
            if (considerable_move(ns, mcgrp, b, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                if (ns.policy.has_rule(FIRST_ACCEPT)) {
                    move(ns, mcgrp);
                    return true;
                }
                else if (ns.policy.has_rule(BEST_ACCEPT)) {
                    if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                        BestM = move_result;
                }
                else {
                    My_Assert(false, "Unknown accept rule!");
                }
            }
        }
        else {
            DEBUG_PRINT("Neighbor Task is dummy Task");
            //j is dummy here but b can't be dummy
            //each start and end location of each route will be considered
            //total 2 x route_nums cases

            int current_start = ns.solution.very_start->next->ID;

            while (current_start != DUMMY) {

#ifdef DEBUG
                attempt_count += 2;
#endif

                // Consider the start location
                int j = current_start;

                if (b != j) { //not the same Task
                    if (considerable_move(ns, mcgrp, b, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Consider the end location
                const int current_route = ns.solution[current_start]->route_id;
                const int current_end = ns.routes[current_route]->end;
                j = current_end;
                if (b != j) {
                    if (considerable_move(ns, mcgrp, b, j) && ns.policy.check_move(mcgrp,ns,move_result)) {
                        if (ns.policy.has_rule(FIRST_ACCEPT)) {
                            move(ns, mcgrp);
                            return true;
                        }
                        else if (ns.policy.has_rule(BEST_ACCEPT)) {
                            if (ns.policy.check_result(mcgrp,ns,move_result, BestM))
                                BestM = move_result;
                        }
                        else {
                            My_Assert(false, "Unknown accept rule!");
                        }
                    }
                }

                // Advance to next route's starting Task
                current_start = ns.solution[current_end]->next->next->ID;
            }
        }
    }

    if (ns.policy.has_rule(FIRST_ACCEPT)) {
        DEBUG_PRINT("No actual move: First Accept Rule");
        return false;
    }
    else if (ns.policy.has_rule(BEST_ACCEPT)) {
        if (BestM.considerable == false) {
            DEBUG_PRINT("No actual move: Best Accept Rule");
            return false;
        }
        else {
            move_result = BestM;
            move(ns, mcgrp);
            return true;
        }
    }
    else {
        My_Assert(false, "Unknown accept rule");
        return false;
    }
}

bool NewSwap::considerable_move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp, int i, int u)
{
    My_Assert(u != i, "two Task need to be different!");
    My_Assert(i >= 1 && i <= mcgrp.actual_task_num,"Wrong Task");
    My_Assert(u >= 1 && u <= mcgrp.actual_task_num,"Wrong Task");


    const int original_u = u;
    const int original_i = i;


    int u_route, i_route;
    u_route = ns.solution[u]->route_id;
    i_route = ns.solution[i]->route_id;

    // check load legality
    int u_load_delta = 0;
    int i_load_delta = 0;
    double vio_load_delta = 0;
    if (u_route != i_route) {

        u_load_delta = mcgrp.inst_tasks[i].demand - mcgrp.inst_tasks[u].demand;
        i_load_delta = mcgrp.inst_tasks[u].demand - mcgrp.inst_tasks[i].demand;

        if (ns.policy.has_rule(FEASIBLE)) {
            if (ns.routes[u_route]->load + u_load_delta > mcgrp.capacity) {
                move_result.reset();
                return false;
            }

            if (ns.routes[i_route]->load + i_load_delta > mcgrp.capacity) {
                move_result.reset();
                return false;
            }
        }
        else if (ns.policy.has_rule(INFEASIBLE)) {
            int pseudo_capacity = ns.policy.get_pseudo_capacity(mcgrp.capacity);
            if (ns.routes[u_route]->load + u_load_delta > pseudo_capacity) {
                move_result.reset();
                return false;
            }

            if (ns.routes[i_route]->load + i_load_delta > pseudo_capacity) {
                move_result.reset();
                return false;
            }

            //u_route vio-load calculate
            if (u_load_delta >= 0) {
                //vio_load_delta may get bigger
                if (ns.routes[u_route]->load + u_load_delta > mcgrp.capacity) {
                    //if result is overload
                    if (ns.routes[u_route]->load >= mcgrp.capacity) {
                        //if the route u already over loaded
                        vio_load_delta += u_load_delta;
                    }
                    else {
                        vio_load_delta += ns.routes[u_route]->load - mcgrp.capacity + u_load_delta;
                    }
                }
            }
            else {
                //vio_load_delta may get less
                if (ns.routes[u_route]->load + u_load_delta > mcgrp.capacity) {
                    //if result is overload
                    vio_load_delta += u_load_delta;
                }
                else {
                    //the result will be not overload
                    if (ns.routes[u_route]->load > mcgrp.capacity) {
                        //the route used to be overload
                        vio_load_delta += mcgrp.capacity - ns.routes[u_route]->load;
                    }
                }
            }

            if (i_load_delta >= 0) {
                //vio_load_delta may get bigger
                if (ns.routes[i_route]->load + i_load_delta > mcgrp.capacity) {
                    //if result is overload
                    if (ns.routes[i_route]->load >= mcgrp.capacity) {
                        //if the route u already over loaded
                        vio_load_delta += i_load_delta;
                    }
                    else {
                        vio_load_delta += ns.routes[i_route]->load - mcgrp.capacity + i_load_delta;
                    }
                }
            }
            else {
                //vio_load_delta may get less
                if (ns.routes[i_route]->load + i_load_delta > mcgrp.capacity) {
                    //if result is overload
                    vio_load_delta += i_load_delta;
                }
                else {
                    //the result will be not overload
                    if (ns.routes[i_route]->load > mcgrp.capacity) {
                        //the route used to be overload
                        vio_load_delta += mcgrp.capacity - ns.routes[i_route]->load;
                    }
                }
            }
        }
    }

    const int t = max(ns.solution[u]->pre->ID, 0);
    const int v = max(ns.solution[u]->next->ID, 0);

    const int h = max(ns.solution[i]->pre->ID, 0);
    const int j = max(ns.solution[i]->next->ID, 0);


    //These two info are used when two tasks are not in the same route
    double u_delta = 0;
    double i_delta = 0;

    double delta = 0;
    vector<vector<RouteInfo::TimeTable>> new_time_tbl{{{-1, -1}}};
    bool allow_infeasible = ns.policy.has_rule(INFEASIBLE) ? true : false;

    if(j == u){
        //...h-i-u-v...
        DEBUG_PRINT("Swap:two tasks are neighbor");
        My_Assert(i_route == u_route,"Close tasks should be in the same route!");
        My_Assert(t == i,"t Task must be equal with i Task!");

        double deltas[4];
        vector<vector<RouteInfo::TimeTable>> time_tbl[4];
        if (mcgrp.is_edge(u) && mcgrp.is_edge(i)) {
            const int u_tilde = mcgrp.inst_tasks[u].inverse;
            const int i_tilde = mcgrp.inst_tasks[i].inverse;

            const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iv = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hu = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u].head_node];

            const double i_tildev = mcgrp.min_cost[mcgrp.inst_tasks[i_tilde].tail_node][mcgrp.inst_tasks[v].head_node];
            const double ui_tilde = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i_tilde].head_node];
            const double hu_tilde = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u_tilde].head_node];
            const double u_tildei = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[i].head_node];


            const double u_tildei_tilde = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[i_tilde].head_node];

            //remove the original connections firstly
            double base = -hi - mcgrp.inst_tasks[i].serv_cost - iu - mcgrp.inst_tasks[u].serv_cost - uv;

            //add new connections secondly
            //...h-u-i-v...
            deltas[0] = hu + mcgrp.inst_tasks[u].serv_cost + ui + mcgrp.inst_tasks[i].serv_cost + iv;
            time_tbl[0] = expected_time_table(ns,mcgrp,u,u,i,i,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0]))
                deltas[0] = 1e20;

            //...h-u-i_tilde-v...
            deltas[1] = hu + mcgrp.inst_tasks[u].serv_cost + ui_tilde + mcgrp.inst_tasks[i_tilde].serv_cost + i_tildev;
            time_tbl[1] = expected_time_table(ns,mcgrp,u,u,i,i_tilde,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[1][0]))
                deltas[1] = 1e20;

            //...h-u_tilde-i-v...
            deltas[2] = hu_tilde + mcgrp.inst_tasks[u_tilde].serv_cost + u_tildei + mcgrp.inst_tasks[i].serv_cost + iv;
            time_tbl[2] = expected_time_table(ns,mcgrp,u,u_tilde,i,i,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[2][0]))
                deltas[2] = 1e20;

            //...h-u_tilde-i_tilde-v...
            deltas[3] = hu_tilde + mcgrp.inst_tasks[u_tilde].serv_cost + u_tildei_tilde + mcgrp.inst_tasks[i_tilde].serv_cost + i_tildev;
            time_tbl[3] = expected_time_table(ns,mcgrp,u,u_tilde,i,i_tilde,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[3][0]))
                deltas[3] = 1e20;

            int cases = -1;
            double best_delta = 1e20;
            for(int cursor = 0;cursor<4 ;cursor++){
                if(deltas[cursor] < best_delta){
                    cases = cursor;
                    best_delta = deltas[cursor];
                    new_time_tbl = time_tbl[cursor];
                }
            }

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }

            //decide which Task should be swapped finally
            switch (cases){
                case -1:
                case 0:
                    u = u;
                    i = i;
                    break;
                case 1:
                    u = u;
                    i = i_tilde;
                    break;
                case 2:
                    u = u_tilde;
                    i = i;
                    break;
                case 3:
                    u = u_tilde;
                    i = i_tilde;
                    break;
                default:
                    My_Assert(false,"unknown cases!");
            }

        }
        else if (mcgrp.is_edge(u)){

            const int u_tilde = mcgrp.inst_tasks[u].inverse;

            const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iv = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hu = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u].head_node];

            const double hu_tilde = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u_tilde].head_node];
            const double u_tildei = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[i].head_node];


            //remove the original connections firstly
            double base = -hi - mcgrp.inst_tasks[i].serv_cost - iu - mcgrp.inst_tasks[u].serv_cost - uv;

            //add new connections secondly
            //...h-u-i-v...
            deltas[0] = hu + mcgrp.inst_tasks[u].serv_cost + ui + mcgrp.inst_tasks[i].serv_cost + iv;
            time_tbl[0] = expected_time_table(ns,mcgrp,u,u,i,i,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0]))
                deltas[0] = 1e20;


            //...h-u_tilde-i-v...
            deltas[1] = hu_tilde + mcgrp.inst_tasks[u_tilde].serv_cost + u_tildei + mcgrp.inst_tasks[i].serv_cost + iv;
            time_tbl[1] = expected_time_table(ns,mcgrp,u,u_tilde,i,i,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[1][0]))
                deltas[1] = 1e20;


            int cases = -1;
            double best_delta = 1e20;
            for(int cursor = 0;cursor<2;cursor++){
                if(deltas[cursor] < best_delta){
                    cases = cursor;
                    best_delta = deltas[cursor];
                    new_time_tbl = time_tbl[cursor];
                }
            }

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }

            //decide which Task should be swapped finally
            switch (cases){
                case -1:
                case 0:
                    u = u;
                    i = i;
                    break;
                case 1:
                    u = u_tilde;
                    i = i;
                    break;
                default:
                    My_Assert(false,"unknown cases!");
            }
        }
        else if (mcgrp.is_edge(i)){
            const int i_tilde = mcgrp.inst_tasks[i].inverse;

            const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iv = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hu = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u].head_node];

            const double i_tildev = mcgrp.min_cost[mcgrp.inst_tasks[i_tilde].tail_node][mcgrp.inst_tasks[v].head_node];
            const double ui_tilde = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i_tilde].head_node];



            //remove the original connections firstly
            double base = -hi - mcgrp.inst_tasks[i].serv_cost - iu - mcgrp.inst_tasks[u].serv_cost - uv;

            //add new connections secondly
            //...h-u-i-v...
            deltas[0] = hu + mcgrp.inst_tasks[u].serv_cost + ui + mcgrp.inst_tasks[i].serv_cost + iv;
            time_tbl[0] = expected_time_table(ns,mcgrp,u,u,i,i,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0]))
                deltas[0] = 1e20;


            //...h-u-i_tilde-v...
            deltas[1] = hu + mcgrp.inst_tasks[u].serv_cost + ui_tilde + mcgrp.inst_tasks[i_tilde].serv_cost + i_tildev;
            time_tbl[1] = expected_time_table(ns,mcgrp,u,u,i,i_tilde,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[1][0]))
                deltas[1] = 1e20;


            int cases = -1;
            double best_delta = 1e20;
            for(int cursor = 0;cursor<2;cursor++){
                if(deltas[cursor] < best_delta){
                    cases = cursor;
                    best_delta = deltas[cursor];
                    new_time_tbl = time_tbl[cursor];
                }
            }

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }

            //decide which Task should be swapped finally
            switch (cases){
                case -1:
                case 0:
                    u = u;
                    i = i;
                    break;
                case 1:
                    u = u;
                    i = i_tilde;
                    break;
                default:
                    My_Assert(false,"unknown cases!");
            }
        }
        else{
            const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double iv = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[v].head_node];
            const double hu = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u].head_node];

            //remove the original connections firstly
            double base = -hi - mcgrp.inst_tasks[i].serv_cost - iu - mcgrp.inst_tasks[u].serv_cost - uv;

            //add new connections secondly
            //...h-u-i-v...
            deltas[0] = hu + mcgrp.inst_tasks[u].serv_cost + ui + mcgrp.inst_tasks[i].serv_cost + iv;
            time_tbl[0] = expected_time_table(ns,mcgrp,u,u,i,i,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0])){
                deltas[0] = 1e20;
            }


            double best_delta = deltas[0];
            if(best_delta != 1e20)
                new_time_tbl = time_tbl[0];

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }
        }
    }
    else if (v == i){
        //...t-u-i-j...
        DEBUG_PRINT("Swap:two tasks are neighbor");
        My_Assert(i_route == u_route,"Close tasks should be in the same route!");
        My_Assert(u == h,"u Task must be equal with h Task!");

        double deltas[4];
        vector<vector<RouteInfo::TimeTable>> time_tbl[4];
        if (mcgrp.is_edge(u) && mcgrp.is_edge(i)) {
            const int u_tilde = mcgrp.inst_tasks[u].inverse;
            const int i_tilde = mcgrp.inst_tasks[i].inverse;

            const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ij = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[j].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double ti = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i].head_node];
            const double uj = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[j].head_node];

            const double ti_tilde = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i_tilde].head_node];
            const double i_tildeu = mcgrp.min_cost[mcgrp.inst_tasks[i_tilde].tail_node][mcgrp.inst_tasks[u].head_node];
            const double iu_tilde = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u_tilde].head_node];
            const double u_tildej = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[j].head_node];

            const double i_tildeu_tilde =
                mcgrp.min_cost[mcgrp.inst_tasks[i_tilde].tail_node][mcgrp.inst_tasks[u_tilde].head_node];

            //remove the original connections firstly
            double base = -tu - mcgrp.inst_tasks[u].serv_cost - ui - mcgrp.inst_tasks[i].serv_cost - ij;

            //add new connections secondly
            //...t-i-u-j...
            deltas[0] = ti + mcgrp.inst_tasks[i].serv_cost + iu + mcgrp.inst_tasks[u].serv_cost + uj;
            time_tbl[0] = expected_time_table(ns,mcgrp,i,i,u,u,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0]))
                deltas[0] = 1e20;

            //...t-i-u_tilde-j...
            deltas[1] = ti + mcgrp.inst_tasks[i].serv_cost + iu_tilde + mcgrp.inst_tasks[u_tilde].serv_cost + u_tildej;
            time_tbl[1] = expected_time_table(ns,mcgrp,i,i,u,u_tilde,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[1][0]))
                deltas[1] = 1e20;

            //...t-i_tilde-u-j...
            deltas[2] = ti_tilde + mcgrp.inst_tasks[i_tilde].serv_cost + i_tildeu + mcgrp.inst_tasks[u].serv_cost + uj;
            time_tbl[2] = expected_time_table(ns,mcgrp,i,i_tilde,u,u,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[2][0]))
                deltas[2] = 1e20;

            //...t-i_tilde-u_tilde-j...
            deltas[3] = ti_tilde + mcgrp.inst_tasks[i_tilde].serv_cost + i_tildeu_tilde + mcgrp.inst_tasks[u_tilde].serv_cost + u_tildej;
            time_tbl[3] = expected_time_table(ns,mcgrp,i,i_tilde,u,u_tilde,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[3][0]))
                deltas[3] = 1e20;

            int cases = -1;
            double best_delta = 1e20;
            for (int cursor = 0; cursor < 4; cursor++) {
                if (deltas[cursor] < best_delta) {
                    cases = cursor;
                    best_delta = deltas[cursor];
                    new_time_tbl = time_tbl[cursor];
                }
            }

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }

            //decide which Task should be swapped finally
            switch (cases) {
                case -1:
                case 0:i = i;
                    u = u;
                    break;
                case 1:i = i;
                    u = u_tilde;
                    break;
                case 2:i = i_tilde;
                    u = u;
                    break;
                case 3:i = i_tilde;
                    u = u_tilde;
                    break;
                default:My_Assert(false, "unknown cases!");
            }
        }
        else if (mcgrp.is_edge(i)){
            const int i_tilde = mcgrp.inst_tasks[i].inverse;

            const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ij = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[j].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double ti = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i].head_node];
            const double uj = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[j].head_node];

            const double ti_tilde = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i_tilde].head_node];
            const double i_tildeu = mcgrp.min_cost[mcgrp.inst_tasks[i_tilde].tail_node][mcgrp.inst_tasks[u].head_node];

            //remove the original connections firstly
            double base = -tu - mcgrp.inst_tasks[u].serv_cost - ui - mcgrp.inst_tasks[i].serv_cost - ij;

            //add new connections secondly
            //...t-i-u-j...
            deltas[0] = ti + mcgrp.inst_tasks[i].serv_cost + iu + mcgrp.inst_tasks[u].serv_cost + uj;
            time_tbl[0] = expected_time_table(ns,mcgrp,i,i,u,u,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0]))
                deltas[0] = 1e20;

            //...t-i_tilde-u-j...
            deltas[1] = ti_tilde + mcgrp.inst_tasks[i_tilde].serv_cost + i_tildeu + mcgrp.inst_tasks[u].serv_cost + uj;
            time_tbl[1] = expected_time_table(ns,mcgrp,i,i_tilde,u,u,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[1][0]))
                deltas[1] = 1e20;

            int cases = -1;
            double best_delta = 1e20;
            for (int cursor = 0; cursor < 2; cursor++) {
                if (deltas[cursor] < best_delta) {
                    cases = cursor;
                    best_delta = deltas[cursor];
                    new_time_tbl = time_tbl[cursor];
                }
            }

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }

            //decide which Task should be swapped finally
            switch (cases) {
                case -1:
                case 0:
                    i = i;
                    u = u;
                    break;
                case 1:
                    i = i_tilde;
                    u = u;
                    break;
                default:My_Assert(false, "unknown cases!");
            }

        }
        else if (mcgrp.is_edge(u)){
            const int u_tilde = mcgrp.inst_tasks[u].inverse;

            const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ij = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[j].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double ti = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i].head_node];
            const double uj = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[j].head_node];

            const double iu_tilde = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u_tilde].head_node];
            const double u_tildej = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[j].head_node];

            //remove the original connections firstly
            double base = -tu - mcgrp.inst_tasks[u].serv_cost - ui - mcgrp.inst_tasks[i].serv_cost - ij;

            //add new connections secondly
            //...t-i-u-j...
            deltas[0] = ti + mcgrp.inst_tasks[i].serv_cost + iu + mcgrp.inst_tasks[u].serv_cost + uj;
            time_tbl[0] = expected_time_table(ns,mcgrp,i,i,u,u,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0]))
                deltas[0] = 1e20;

            //...h-i-u_tilde-j...
            deltas[1] = ti + mcgrp.inst_tasks[i].serv_cost + iu_tilde + mcgrp.inst_tasks[u_tilde].serv_cost + u_tildej;
            time_tbl[1] = expected_time_table(ns,mcgrp,i,i,u,u_tilde,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[1][0]))
                deltas[1] = 1e20;

            int cases = -1;
            double best_delta = 1e20;
            for (int cursor = 0; cursor < 2; cursor++) {
                if (deltas[cursor] < best_delta) {
                    cases = cursor;
                    best_delta = deltas[cursor];
                    new_time_tbl = time_tbl[cursor];
                }
            }

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }

            //decide which Task should be swapped finally
            switch (cases) {
                case -1:
                case 0:i = i;
                    u = u;
                    break;
                case 1:i = i;
                    u = u_tilde;
                    break;
                default:My_Assert(false, "unknown cases!");
            }
        }
        else{
            const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ij = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[j].head_node];
            const double iu = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[u].head_node];
            const double ui = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[i].head_node];
            const double ti = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i].head_node];
            const double uj = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[j].head_node];


            //remove the original connections firstly
            double base = -tu - mcgrp.inst_tasks[u].serv_cost - ui - mcgrp.inst_tasks[i].serv_cost - ij;

            //add new connections secondly
            //...t-i-u-j...
            deltas[0] = ti + mcgrp.inst_tasks[i].serv_cost + iu + mcgrp.inst_tasks[u].serv_cost + uj;
            time_tbl[0] = expected_time_table(ns,mcgrp,i,i,u,u,allow_infeasible);
            if(!mcgrp.isTimeTableFeasible(time_tbl[0][0])){
                deltas[0] = 1e20;
            }

            double best_delta = deltas[0];
            if(best_delta != 1e20)
                new_time_tbl = time_tbl[0];

            if(best_delta == 1e20){
                delta = best_delta;
            }else{
                delta = best_delta + base;
            }
        }
    }
    else{
        //...t-u-v...h-i-j...
        // OR
        //...h-i-j...t-u-v...
        int u_tilde = u;
        int i_tilde = i;
        //Total four situations need to be considered here
        if (mcgrp.is_edge(u)) {
            u_tilde = mcgrp.inst_tasks[u].inverse;
        }
        if (mcgrp.is_edge(i)) {
            i_tilde = mcgrp.inst_tasks[i].inverse;
        }

        vector<vector<RouteInfo::TimeTable>> time_tbl[4];
        const double tu = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[u].head_node];
        const double uv = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[v].head_node];
        const double hi = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[i].head_node];
        const double ij = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[j].head_node];
        const double ti = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i].head_node];
        const double iv = mcgrp.min_cost[mcgrp.inst_tasks[i].tail_node][mcgrp.inst_tasks[v].head_node];
        const double hu = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u].head_node];
        const double uj = mcgrp.min_cost[mcgrp.inst_tasks[u].tail_node][mcgrp.inst_tasks[j].head_node];

        const double ti_tilde = mcgrp.min_cost[mcgrp.inst_tasks[t].tail_node][mcgrp.inst_tasks[i_tilde].head_node];
        const double i_tildev = mcgrp.min_cost[mcgrp.inst_tasks[i_tilde].tail_node][mcgrp.inst_tasks[v].head_node];

        const double hu_tilde = mcgrp.min_cost[mcgrp.inst_tasks[h].tail_node][mcgrp.inst_tasks[u_tilde].head_node];
        const double u_tildej = mcgrp.min_cost[mcgrp.inst_tasks[u_tilde].tail_node][mcgrp.inst_tasks[j].head_node];

        double u_delta_base = -tu - uv - mcgrp.inst_tasks[u].serv_cost;
        double i_delta_base = -hi - ij - mcgrp.inst_tasks[i].serv_cost;
        double delta_tmp[4][2]; // [...][0] for u_route, [...][1] for i_route

        // u,i
        delta_tmp[0][0] = u_delta_base + ti + iv + mcgrp.inst_tasks[i].serv_cost;
        delta_tmp[0][1] = i_delta_base + hu + uj + mcgrp.inst_tasks[u].serv_cost;
        time_tbl[0] = expected_time_table(ns,mcgrp,u,u,i,i,allow_infeasible);
        if(!mcgrp.isTimeTableFeasible(time_tbl[0][0])){
            delta_tmp[0][0] = 1e20;
            delta_tmp[0][1] = 1e20;
        }

        // u, i_tilde
        delta_tmp[1][0] = u_delta_base + ti_tilde + i_tildev + mcgrp.inst_tasks[i_tilde].serv_cost;
        delta_tmp[1][1] = i_delta_base + hu + uj + mcgrp.inst_tasks[u].serv_cost;
        time_tbl[1] = expected_time_table(ns,mcgrp,u,u,i,i_tilde,allow_infeasible);
        if(!mcgrp.isTimeTableFeasible(time_tbl[1][0])){
            delta_tmp[1][0] = 1e20;
            delta_tmp[1][1] = 1e20;
        }

        // u_tilde, i
        delta_tmp[2][0] = u_delta_base + ti + iv + mcgrp.inst_tasks[i].serv_cost;
        delta_tmp[2][1] = i_delta_base + hu_tilde + u_tildej + mcgrp.inst_tasks[u_tilde].serv_cost;
        time_tbl[2] = expected_time_table(ns,mcgrp,u,u_tilde,i,i,allow_infeasible);
        if(!mcgrp.isTimeTableFeasible(time_tbl[2][0])){
            delta_tmp[2][0] = 1e20;
            delta_tmp[2][1] = 1e20;
        }

        // u_tilde,i_tilde
        delta_tmp[3][0] = u_delta_base + ti_tilde + i_tildev + mcgrp.inst_tasks[i_tilde].serv_cost;
        delta_tmp[3][1] = i_delta_base + hu_tilde + u_tildej + mcgrp.inst_tasks[u_tilde].serv_cost;
        time_tbl[3] = expected_time_table(ns,mcgrp,u,u_tilde,i,i_tilde,allow_infeasible);
        if(!mcgrp.isTimeTableFeasible(time_tbl[3][0])){
            delta_tmp[3][0] = 1e20;
            delta_tmp[3][1] = 1e20;
        }

        delta = 1e20;
        int cur = -1;
        for(int k = 0; k< 4; k++){
            if(delta_tmp[k][0] == 1e20 || delta_tmp[k][1] == 1e20)
                continue;

            if(delta_tmp[k][0] + delta_tmp[k][1] < delta){
                u_delta = delta_tmp[k][0];
                i_delta = delta_tmp[k][1];
                delta = u_delta + i_delta;
                new_time_tbl = time_tbl[k];
                cur = k;
            }
        }

        switch (cur) {
            case -1:
            case 0:
                u = u;
                i = i;
                break;
            case 1:
                u = u;
                i = i_tilde;
                break;
            case 2:
                u = u_tilde;
                i = i;
                break;
            case 3:
                u = u_tilde;
                i = i_tilde;
                break;
            default:
                My_Assert(false,"unknown cases!");
        }
    }

    if(delta == 1e20){
        move_result.reset();
        return false;
    }

    if(allow_infeasible){
        if(ns.policy.check_time_window(mcgrp,new_time_tbl)){
            move_result.reset();
            return false;
        }
    }

    move_result.choose_tasks(original_u, original_i);

    move_result.move_arguments.push_back(original_u);
    move_result.move_arguments.push_back(original_i);

    move_result.total_number_of_routes = ns.routes.activated_route_id.size();

    if (u_route == i_route) {

        move_result.num_affected_routes = 1;
        move_result.delta = delta;

        move_result.route_id.push_back(u_route);

        move_result.route_lens.push_back(ns.routes[u_route]->length + delta);

        move_result.route_loads.push_back(ns.routes[u_route]->load);

        move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers);

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;
        move_result.move_arguments.push_back(u);
        move_result.move_arguments.push_back(i);

        move_result.considerable = true;

        move_result.route_time_tbl.emplace_back(new_time_tbl[0]);


        if(ns.policy.has_rule(INFEASIBLE)){
            auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
            auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

            move_result.vio_load_delta = vio_load_delta;
            move_result.vio_time_delta = new_time_info_u.second - old_time_info_u.second;
            move_result.vio_time_custom_num_delta = new_time_info_u.first - old_time_info_u.first;
        }else{
            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = 0;
            move_result.vio_time_custom_num_delta = 0;
        }

        return true;
    }
    else {
        // Different routes
        move_result.num_affected_routes = 2;
        move_result.delta = delta;

        move_result.route_id.push_back(u_route);
        move_result.route_id.push_back(i_route);

        move_result.route_lens.push_back(ns.routes[u_route]->length + u_delta);
        move_result.route_lens.push_back(ns.routes[i_route]->length + i_delta);

        int u_route_load = ns.routes[u_route]->load - mcgrp.inst_tasks[u].demand + mcgrp.inst_tasks[i].demand;
        int i_route_load = ns.routes[i_route]->load - mcgrp.inst_tasks[i].demand + mcgrp.inst_tasks[u].demand;
        move_result.route_loads.push_back(u_route_load);
        move_result.route_loads.push_back(i_route_load);


        move_result.route_custs_num.push_back(ns.routes[u_route]->num_customers);
        move_result.route_custs_num.push_back(ns.routes[i_route]->num_customers);

        move_result.new_total_route_length = ns.cur_solution_cost + move_result.delta;

        move_result.move_arguments.push_back(u);
        move_result.move_arguments.push_back(i);

        move_result.vio_load_delta = vio_load_delta;

        move_result.considerable = true;

        move_result.route_time_tbl = new_time_tbl;


        if(ns.policy.has_rule(INFEASIBLE)){
            move_result.vio_load_delta = vio_load_delta;

            auto old_time_info_i = mcgrp.get_vio_time(ns.routes[i_route]->time_table);
            auto new_time_info_i = mcgrp.get_vio_time(move_result.route_time_tbl[0]);

            auto old_time_info_u = mcgrp.get_vio_time(ns.routes[u_route]->time_table);
            auto new_time_info_u = mcgrp.get_vio_time(move_result.route_time_tbl[1]);

            move_result.vio_time_delta =
                new_time_info_i.second + new_time_info_u.second - old_time_info_i.second - old_time_info_u.second;
            move_result.vio_time_custom_num_delta =
                new_time_info_i.first + new_time_info_u.first - old_time_info_i.first - old_time_info_u.first;
        }else{
            move_result.vio_load_delta = 0;
            move_result.vio_time_delta = 0;
            move_result.vio_time_custom_num_delta = 0;
        }

        return true;
    }

}

void NewSwap::move(HighSpeedNeighBorSearch &ns, const MCGRP &mcgrp)
{

#ifdef DEBUG
    hit_count++;
#endif

    DEBUG_PRINT("execute a swap move");

    My_Assert(move_result.considerable,"Invalid predictions");
    My_Assert(move_result.move_arguments.size() == 4, "Incorrect move arguments!");


    const int original_u = move_result.move_arguments[0];
    const int original_i = move_result.move_arguments[1];
    const int actual_u = move_result.move_arguments[2];
    const int actual_i = move_result.move_arguments[3];

    My_Assert(original_u >= 1 && original_u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(actual_u >= 1 && actual_u <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(original_i >= 1 && original_i <= mcgrp.actual_task_num,"Wrong arguments");
    My_Assert(actual_i >= 1 && actual_i <= mcgrp.actual_task_num,"Wrong arguments");

    if(move_result.num_affected_routes == 2){
        const int u_route = move_result.route_id[0];
        const int i_route = move_result.route_id[1];
        if(mcgrp.is_edge(original_u) && !mcgrp.is_edge(original_i)){
            ns.routes[u_route]->num_edges -= 1;
            ns.routes[i_route]->num_edges += 1;
        }else if(!mcgrp.is_edge(original_u) && mcgrp.is_edge(original_i)){
            ns.routes[u_route]->num_edges += 1;
            ns.routes[i_route]->num_edges -= 1;
        }
    }

    if(original_u != actual_u && original_i != actual_i) {
        //Edge Task has been moved
        //both Task switch
        My_Assert(mcgrp.inst_tasks[original_u].inverse == actual_u, "This is not the inverse");
        My_Assert(mcgrp.inst_tasks[original_i].inverse == actual_i,"Wrong tasks");
        My_Assert(ns.solution[original_u]->next != nullptr && ns.solution[actual_u]->next == nullptr,"Wrong arguments");
        My_Assert(ns.solution[original_i]->next != nullptr && ns.solution[actual_i]->next == nullptr,"Wrong arguments");

        if(move_result.num_affected_routes == 1) {
            //No need to change route id
            const int route_id = move_result.route_id[0];
            My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes[route_id]->num_customers > 1,"Invalid route");

            ns.routes[route_id]->length = move_result.route_lens[0];

            ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[route_id]->start = actual_i;
            }
            else if (ns.solution[original_u]->next->ID < 0){
                ns.routes[route_id]->end = actual_i;
            }

            if (ns.solution[original_i]->pre->ID < 0){
                ns.routes[route_id]->start = actual_u;
            }
            else if(ns.solution[original_i]->next->ID < 0){
                ns.routes[route_id]->end = actual_u;
            }


            ns.solution[actual_u]->route_id = route_id;
            ns.solution[actual_i]->route_id = route_id;
        }
        else{
            const int u_route = move_result.route_id[0];
            const int i_route = move_result.route_id[1];

            My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[u_route]->length = move_result.route_lens[0];
            ns.routes[i_route]->length = move_result.route_lens[1];

            ns.routes[u_route]->load = move_result.route_loads[0];
            ns.routes[i_route]->load = move_result.route_loads[1];

            ns.routes[u_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[i_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[u_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[i_route]->time_table = move_result.route_time_tbl[1];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[u_route]->start = actual_i;
            }

            if (ns.solution[original_u]->next->ID < 0){
                ns.routes[u_route]->end = actual_i;
            }

            if (ns.solution[original_i]->pre->ID < 0){
                ns.routes[i_route]->start = actual_u;
            }

            if(ns.solution[original_i]->next->ID < 0){
                ns.routes[i_route]->end = actual_u;
            }

            ns.solution[actual_u]->route_id = i_route;
            ns.solution[actual_i]->route_id = u_route;
        }

        //handle solution
        if(ns.solution[original_u]->next == ns.solution[original_i]){
            //...t-u-i-j...
            const auto original_i_next = ns.solution[original_i]->next;
            const auto original_u_pre = ns.solution[original_u]->pre;

            original_i_next->pre = ns.solution[actual_u];
            original_u_pre->next = ns.solution[actual_i];

            ns.solution[actual_i]->pre = original_u_pre;
            ns.solution[actual_i]->next = ns.solution[actual_u];
            ns.solution[actual_u]->pre = ns.solution[actual_i];
            ns.solution[actual_u]->next = original_i_next;
        }
        else if (ns.solution[original_i]->next == ns.solution[original_u]){
            //...h-i-u-v...
            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_u_next = ns.solution[original_u]->next;

            original_i_pre->next = ns.solution[actual_u];
            original_u_next->pre = ns.solution[actual_i];

            ns.solution[actual_i]->pre = ns.solution[actual_u];
            ns.solution[actual_i]->next = original_u_next;
            ns.solution[actual_u]->pre = original_i_pre;
            ns.solution[actual_u]->next = ns.solution[actual_i];
        }
        else{
            //...h-i-j...
            //...t-u-v...
            ns.solution[original_i]->pre->next = ns.solution[actual_u];
            ns.solution[original_i]->next->pre = ns.solution[actual_u];

            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_i_next = ns.solution[original_i]->next;
            ns.solution[actual_i]->pre = ns.solution[original_u]->pre;
            ns.solution[actual_i]->next = ns.solution[original_u]->next;

            ns.solution[original_u]->next->pre = ns.solution[actual_i];
            ns.solution[original_u]->pre->next = ns.solution[actual_i];

            ns.solution[actual_u]->pre = original_i_pre;
            ns.solution[actual_u]->next = original_i_next;
        }

        //clear original i.u Task info
        ns.solution[original_i]->clear();
        ns.solution[original_u]->clear();
    }
    else if (original_u != actual_u){
        //single Task switch
        My_Assert(mcgrp.inst_tasks[original_u].inverse == actual_u,"Wrong tasks");
        My_Assert(ns.solution[original_u]->next != nullptr && ns.solution[actual_u]->next == nullptr,
                  "Wrong arguments");

        if(move_result.num_affected_routes == 1) {
            //No need to change route id
            const int route_id = move_result.route_id[0];
            My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes[route_id]->num_customers > 1,"Invalid route");

            ns.routes[route_id]->length = move_result.route_lens[0];

            ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[route_id]->start = original_i;
            }
            else if (ns.solution[original_u]->next->ID < 0){
                ns.routes[route_id]->end = original_i;
            }

            if (ns.solution[original_i]->pre->ID < 0){
                ns.routes[route_id]->start = actual_u;
            }
            else if(ns.solution[original_i]->next->ID < 0){
                ns.routes[route_id]->end = actual_u;
            }

            ns.solution[actual_u]->route_id = route_id;
        }
        else{
            const int u_route = move_result.route_id[0];
            const int i_route = move_result.route_id[1];

            My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[u_route]->length = move_result.route_lens[0];
            ns.routes[i_route]->length = move_result.route_lens[1];

            ns.routes[u_route]->load = move_result.route_loads[0];
            ns.routes[i_route]->load = move_result.route_loads[1];

            ns.routes[u_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[i_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[u_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[i_route]->time_table = move_result.route_time_tbl[1];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[u_route]->start = original_i;
            }

            if (ns.solution[original_u]->next->ID < 0){
                ns.routes[u_route]->end = original_i;
            }

            if (ns.solution[original_i]->pre->ID < 0){
                ns.routes[i_route]->start = actual_u;
            }

            if(ns.solution[original_i]->next->ID < 0){
                ns.routes[i_route]->end = actual_u;
            }

            ns.solution[actual_u]->route_id = i_route;
            ns.solution[original_i]->route_id = u_route;
        }

        //handle solution
        if(ns.solution[original_u]->next == ns.solution[original_i]){
            //...t-u-i-j...
            const auto original_i_next = ns.solution[original_i]->next;
            const auto original_u_pre = ns.solution[original_u]->pre;

            original_i_next->pre = ns.solution[actual_u];
            original_u_pre->next = ns.solution[original_i];

            ns.solution[original_i]->pre = original_u_pre;
            ns.solution[original_i]->next = ns.solution[actual_u];
            ns.solution[actual_u]->pre = ns.solution[original_i];
            ns.solution[actual_u]->next = original_i_next;
        }
        else if (ns.solution[original_i]->next == ns.solution[original_u]){
            //...h-i-u-v...
            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_u_next = ns.solution[original_u]->next;

            original_i_pre->next = ns.solution[actual_u];
            original_u_next->pre = ns.solution[original_i];

            ns.solution[original_i]->pre = ns.solution[actual_u];
            ns.solution[original_i]->next = original_u_next;
            ns.solution[actual_u]->pre = original_i_pre;
            ns.solution[actual_u]->next = ns.solution[original_i];
        }
        else{
            //...h-i-j...
            //...t-u-v...
            ns.solution[original_i]->pre->next = ns.solution[actual_u];
            ns.solution[original_i]->next->pre = ns.solution[actual_u];

            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_i_next = ns.solution[original_i]->next;
            ns.solution[original_i]->pre = ns.solution[original_u]->pre;
            ns.solution[original_i]->next = ns.solution[original_u]->next;

            ns.solution[original_u]->next->pre = ns.solution[original_i];
            ns.solution[original_u]->pre->next = ns.solution[original_i];

            ns.solution[actual_u]->pre = original_i_pre;
            ns.solution[actual_u]->next = original_i_next;
        }

        //clear original u Task info
        ns.solution[original_u]->clear();
    }
    else if (original_i != actual_i){
        //single Task switch
        My_Assert(mcgrp.inst_tasks[original_i].inverse == actual_i,"Wrong tasks");
        My_Assert(ns.solution[original_i]->next != nullptr && ns.solution[actual_i]->next == nullptr,
                  "Wrong arguments");


        if(move_result.num_affected_routes == 1) {
            //No need to change route id
            const int route_id = move_result.route_id[0];
            My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes[route_id]->num_customers > 1,"Invalid route");

            ns.routes[route_id]->length = move_result.route_lens[0];

            ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[route_id]->start = actual_i;
            }
            else if (ns.solution[original_u]->next->ID < 0){
                ns.routes[route_id]->end = actual_i;
            }

            if (ns.solution[original_i]->pre->ID < 0){
                ns.routes[route_id]->start = original_u;
            }
            else if(ns.solution[original_i]->next->ID < 0){
                ns.routes[route_id]->end = original_u;
            }

            ns.solution[actual_i]->route_id = route_id;
        }
        else{
            const int u_route = move_result.route_id[0];
            const int i_route = move_result.route_id[1];

            My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[u_route]->length = move_result.route_lens[0];
            ns.routes[i_route]->length = move_result.route_lens[1];

            ns.routes[u_route]->load = move_result.route_loads[0];
            ns.routes[i_route]->load = move_result.route_loads[1];

            ns.routes[u_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[i_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[u_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[i_route]->time_table = move_result.route_time_tbl[1];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[u_route]->start = actual_i;
            }

            if (ns.solution[original_u]->next->ID < 0){
                ns.routes[u_route]->end = actual_i;
            }

            if (ns.solution[original_i]->pre->ID < 0){
                ns.routes[i_route]->start = original_u;
            }

            if(ns.solution[original_i]->next->ID < 0){
                ns.routes[i_route]->end = original_u;
            }

            ns.solution[original_u]->route_id = i_route;
            ns.solution[actual_i]->route_id = u_route;
        }

        //handle solution
        if(ns.solution[original_u]->next == ns.solution[original_i]){
            //...t-u-i-j...
            const auto original_i_next = ns.solution[original_i]->next;
            const auto original_u_pre = ns.solution[original_u]->pre;

            original_i_next->pre = ns.solution[original_u];
            original_u_pre->next = ns.solution[actual_i];

            ns.solution[actual_i]->pre = original_u_pre;
            ns.solution[actual_i]->next = ns.solution[original_u];
            ns.solution[original_u]->pre = ns.solution[actual_i];
            ns.solution[original_u]->next = original_i_next;
        }
        else if (ns.solution[original_i]->next == ns.solution[original_u]){
            //...h-i-u-v...
            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_u_next = ns.solution[original_u]->next;

            original_i_pre->next = ns.solution[original_u];
            original_u_next->pre = ns.solution[actual_i];

            ns.solution[actual_i]->pre = ns.solution[original_u];
            ns.solution[actual_i]->next = original_u_next;
            ns.solution[original_u]->pre = original_i_pre;
            ns.solution[original_u]->next = ns.solution[actual_i];
        }
        else{
            //...h-i-j...
            //...t-u-v...
            ns.solution[original_i]->pre->next = ns.solution[original_u];
            ns.solution[original_i]->next->pre = ns.solution[original_u];

            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_i_next = ns.solution[original_i]->next;
            ns.solution[actual_i]->pre = ns.solution[original_u]->pre;
            ns.solution[actual_i]->next = ns.solution[original_u]->next;

            ns.solution[original_u]->next->pre = ns.solution[actual_i];
            ns.solution[original_u]->pre->next = ns.solution[actual_i];

            ns.solution[original_u]->pre = original_i_pre;
            ns.solution[original_u]->next = original_i_next;
        }

        //clear original i Task info
        ns.solution[original_i]->clear();
    }
    else{
        //no Task switch
        //Modify routes info
        My_Assert(original_i == actual_i && original_u == actual_u,"Wrong arguments");
        if(move_result.num_affected_routes == 1) {
            //No need to change route id
            const int route_id = move_result.route_id[0];
            My_Assert(ns.routes.activated_route_id.find(route_id)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes[route_id]->num_customers > 1,"Invalid route");

            ns.routes[route_id]->length = move_result.route_lens[0];

            ns.routes[route_id]->time_table = move_result.route_time_tbl[0];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[route_id]->start = original_i;
            }
            else if(ns.solution[original_u]->next->ID < 0){
                ns.routes[route_id]->end = original_i;
            }

            if(ns.solution[original_i]->pre->ID < 0){
                ns.routes[route_id]->start = original_u;
            }
            else if(ns.solution[original_i]->next->ID < 0){
                ns.routes[route_id]->end = original_u;
            }

        }
        else{
            const int u_route = move_result.route_id[0];
            const int i_route = move_result.route_id[1];

            My_Assert(ns.routes.activated_route_id.find(u_route)!=ns.routes.activated_route_id.end(),"Invalid route");
            My_Assert(ns.routes.activated_route_id.find(i_route)!=ns.routes.activated_route_id.end(),"Invalid route");

            ns.routes[u_route]->length = move_result.route_lens[0];
            ns.routes[i_route]->length = move_result.route_lens[1];

            ns.routes[u_route]->load = move_result.route_loads[0];
            ns.routes[i_route]->load = move_result.route_loads[1];

            ns.routes[u_route]->num_customers = move_result.route_custs_num[0];
            ns.routes[i_route]->num_customers = move_result.route_custs_num[1];

            ns.routes[u_route]->time_table = move_result.route_time_tbl[0];
            ns.routes[i_route]->time_table = move_result.route_time_tbl[1];

            if(ns.solution[original_u]->pre->ID < 0){
                ns.routes[u_route]->start = original_i;
            }

            if(ns.solution[original_u]->next->ID < 0){
                ns.routes[u_route]->end = original_i;
            }

            if(ns.solution[original_i]->pre->ID < 0){
                ns.routes[i_route]->start = original_u;
            }

            if(ns.solution[original_i]->next->ID < 0){
                ns.routes[i_route]->end = original_u;
            }

            ns.solution[original_u]->route_id = i_route;
            ns.solution[original_i]->route_id = u_route;
        }

        //handle solution
        if(ns.solution[original_u]->next == ns.solution[original_i]){
            //...t-u-i-j...
            const auto original_i_next = ns.solution[original_i]->next;
            const auto original_u_pre = ns.solution[original_u]->pre;

            original_i_next->pre = ns.solution[original_u];
            original_u_pre->next = ns.solution[original_i];

            ns.solution[original_i]->pre = original_u_pre;
            ns.solution[original_i]->next = ns.solution[original_u];
            ns.solution[original_u]->pre = ns.solution[original_i];
            ns.solution[original_u]->next = original_i_next;
        }
        else if (ns.solution[original_i]->next == ns.solution[original_u]){
            //...h-i-u-v...
            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_u_next = ns.solution[original_u]->next;

            original_i_pre->next = ns.solution[original_u];
            original_u_next->pre = ns.solution[original_i];

            ns.solution[original_i]->pre = ns.solution[original_u];
            ns.solution[original_i]->next = original_u_next;
            ns.solution[original_u]->pre = original_i_pre;
            ns.solution[original_u]->next = ns.solution[original_i];
        }
        else{
            //...h-i-j...
            //...t-u-v...
            ns.solution[original_i]->pre->next = ns.solution[original_u];
            ns.solution[original_i]->next->pre = ns.solution[original_u];

            const auto original_i_pre = ns.solution[original_i]->pre;
            const auto original_i_next = ns.solution[original_i]->next;
            ns.solution[original_i]->pre = ns.solution[original_u]->pre;
            ns.solution[original_i]->next = ns.solution[original_u]->next;

            ns.solution[original_u]->next->pre = ns.solution[original_i];
            ns.solution[original_u]->pre->next = ns.solution[original_i];

            ns.solution[original_u]->pre = original_i_pre;
            ns.solution[original_u]->next = original_i_next;
        }
    }


    //modify global info
    ns.cur_solution_cost += move_result.delta;
    ns.total_vio_load += move_result.vio_load_delta;
    ns.total_vio_time += move_result.vio_time_delta;
    My_Assert(ns.valid_sol(mcgrp),"Prediction wrong!");


    ns.trace(mcgrp);

    move_result.reset();
}

vector<vector<RouteInfo::TimeTable>> NewSwap::expected_time_table(HighSpeedNeighBorSearch &ns,
                                                                  const MCGRP &mcgrp,
                                                                  int u,
                                                                  int u_tilde,
                                                                  int i,
                                                                  int i_tilde,
                                                                  bool allow_infeasible)
{
    vector<vector<RouteInfo::TimeTable>>
        res(2,vector<RouteInfo::TimeTable>({{-1, -1}}));

    const int i_route = ns.solution[i]->route_id;
    const int u_route = ns.solution[u]->route_id;

    if(u_route == i_route){
        vector<int> route_task;
        for(const auto& node:ns.routes[u_route]->time_table){
            if(node.task == u){
                route_task.push_back(i_tilde);
            }else if(node.task == i){
                route_task.push_back(u_tilde);
            }else{
                route_task.push_back(node.task);
            }
        }
        vector<int> time_tbl = mcgrp.cal_arrive_time(route_task);
        vector<RouteInfo::TimeTable> intermediate;

        My_Assert(time_tbl.size() == route_task.size(), "wrong time table");
        for(int k = 0; k<time_tbl.size(); k++){
            if(!allow_infeasible && time_tbl[k] > mcgrp.inst_tasks[route_task[k]].time_window.second)
                return res;
            intermediate.push_back({route_task[k],time_tbl[k]});
        }

        res[0] = intermediate;
        res[1] = intermediate;

        return res;
    }
    else{
        vector<int> route_task_u;
        vector<int> route_task_i;
        for(const auto& node:ns.routes[u_route]->time_table){
            if(node.task == u){
                route_task_u.push_back(i_tilde);
            }else{
                route_task_u.push_back(node.task);
            }
        }

        for(const auto& node:ns.routes[i_route]->time_table){
            if(node.task == i){
                route_task_i.push_back(u_tilde);
            }else{
                route_task_i.push_back(node.task);
            }
        }

        vector<int> time_tbl_u = mcgrp.cal_arrive_time(route_task_u);
        vector<int> time_tbl_i = mcgrp.cal_arrive_time(route_task_i);

        My_Assert(time_tbl_u.size() == route_task_u.size(), "wrong time table");
        My_Assert(time_tbl_i.size() == route_task_i.size(), "wrong time table");

        vector<RouteInfo::TimeTable> intermediate_u;
        vector<RouteInfo::TimeTable> intermediate_i;
        for(int k = 0; k<time_tbl_u.size(); k++){
            if(!allow_infeasible && time_tbl_u[k] > mcgrp.inst_tasks[route_task_u[k]].time_window.second)
                return res;
            intermediate_u.push_back({route_task_u[k],time_tbl_u[k]});
        }

        for(int k = 0; k<time_tbl_i.size(); k++){
            if(!allow_infeasible && time_tbl_i[k] > mcgrp.inst_tasks[route_task_i[k]].time_window.second)
                return res;
            intermediate_i.push_back({route_task_i[k],time_tbl_i[k]});
        }

        res[0] = intermediate_u;
        res[1] = intermediate_i;

        return res;
    }


}


bool NewSwap::update_score(HighSpeedNeighBorSearch &ns)
{
    if(move_result.delta > 0) return false;

    const int actual_u = move_result.move_arguments[2];
    const int actual_i = move_result.move_arguments[3];

    double penalty = (ns.best_solution_cost / ns.cur_solution_cost) * (-move_result.delta);

    if(ns.solution[actual_u]->next->ID == actual_i){
        // two tasks are neighbor
        ns.score_matrix[actual_u][actual_i] += penalty;
        ns.score_matrix[max(0,ns.solution[actual_u]->pre->ID)][actual_u] += penalty;
        ns.score_matrix[actual_i][max(0,ns.solution[actual_i]->next->ID)] += penalty;
    }else if(ns.solution[actual_i]->next->ID == actual_u){
        ns.score_matrix[actual_i][actual_u] += penalty;
        ns.score_matrix[max(0,ns.solution[actual_i]->pre->ID)][actual_i] += penalty;
        ns.score_matrix[actual_u][max(0,ns.solution[actual_u]->next->ID)] += penalty;
    }else{
        ns.score_matrix[max(0,ns.solution[actual_u]->pre->ID)][actual_u] += penalty;
        ns.score_matrix[actual_u][max(0,ns.solution[actual_u]->next->ID)] += penalty;

        ns.score_matrix[actual_i][max(0,ns.solution[actual_i]->next->ID)] += penalty;
        ns.score_matrix[max(0,ns.solution[actual_i]->pre->ID)][actual_i] += penalty;
    }

    return true;
}
