//
// Created by luke on 2019/11/29.
//


#ifndef INITILIZER_H
#define INITILIZER_H

#include "instance.h"
#include "utils.h"


class LocalSearch;
class Policy;

class PathConstructor{
protected:
    string name;
    const MCGRPTW& mcgrp;

public:
    PathConstructor(const MCGRPTW& mcgrp_, const string& name_);

    // mode : [allow_infeasible, feasible]
    virtual Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") = 0;
};

class CWScanner : public PathConstructor{
public:
    // Terminology used by Clarke & Wright
    struct Cell{
        int pair1;
        int pair2;
        int pair1mode;
        int pair2mode;
        int value;

        Cell(int pair1, int pair2, int pair1Mode, int pair2Mode, int value)
            : pair1(pair1), pair2(pair2), pair1mode(pair1Mode), pair2mode(pair2Mode), value(value)
        {}

        static bool cmp(const Cell& c1, const Cell& c2){
            return c1.value > c2.value;
        }
    };

    struct Route{
        int cost;
        int load;
        vector<int> sequence;
        vector<int> arrive_time;

        void clear(){
            cost = 0;
            load = 0;
            sequence.clear();
            arrive_time.clear();
        }

        Route(int cost, int load, const vector<int> &sequence, const vector<int> &arriveTime)
            : cost(cost), load(load), sequence(sequence), arrive_time(arriveTime)
        {}

    };

    vector<int> task_list;
    vector<Cell> saving_list;
    vector<Route> routes;
    vector<int> task_pos;
    vector<vector<int>> d;
    vector<vector<int>> t;

    static bool isinterior_(int task, Route route){
        return task != route.sequence.front() and task != route.sequence.back();
    }

    static bool isfirst_(int task, Route route){
        return task == route.sequence.front();
    }

    static bool islast_(int task, Route route){
        return task == route.sequence.back();
    }

    vector<int> dump_(){
        vector<int> solution {0};
        for (const auto& route : routes){
            if(route.sequence.size()) {
                solution.insert(solution.end(),route.sequence.begin(),route.sequence.end());
                solution.push_back(0);
            }
        }

        return solution;
    }


    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;


    void update_savings_(int new_route_id){
        const auto & route = routes[new_route_id];

        if (route.sequence.size() <= 1){
            cerr<<"invalid update";
            exit(-1);
        }

        bool update = false;

        for (auto& saving : saving_list){
            bool pair1_is_last = false;
            bool pair2_is_first = false;

            int pair1_route_id = -1;
            vector<int> pair1_possible = mcgrp.exclude(saving.pair1);
            pair1_possible.push_back(saving.pair1);
            for( auto i : pair1_possible){
                pair1_route_id = max(task_pos[i], pair1_route_id);
            }
            const auto& pair1_route = routes[pair1_route_id];

            int pair2_route_id = -1;
            vector<int> pair2_possible = mcgrp.exclude(saving.pair2);
            pair2_possible.push_back(saving.pair2);
            for(auto i : pair2_possible){
                pair2_route_id = max(task_pos[i], pair2_route_id);
            }
            const auto& pair2_route = routes[pair2_route_id];

            if (pair1_route_id != new_route_id && pair2_route_id != new_route_id){
                // this pair doesn't need to update
                continue;
            }

            if (pair1_route.sequence.back() == saving.pair1){
                pair1_is_last = true;
            }

            for(auto tmp : mcgrp.exclude(pair1_route.sequence.back())){
                if (saving.pair1 == tmp){
                    pair1_is_last = true;
                }
            }

            if( pair2_route.sequence.back() == saving.pair2){
                pair2_is_first = true;
            }

            for(auto tmp : mcgrp.exclude(pair2_route.sequence.back())){
                if(saving.pair2 == tmp){
                    pair2_is_first = true;
                }
            }

            if( not (pair1_is_last ^ pair2_is_first)){
                if (pair2_is_first and pair1_is_last){
                    // this pair in the same route
                    saving.value = DBL_MIN;
                }
                continue;
            }

            if(not (pair1_is_last and pair2_is_first)){
                // This pair cannot imporve solution in the later
                saving.value = DBL_MIN;
                continue;
            }

            auto y = saving.pair1;
            auto pre_y = 0;
            if (pair1_route.sequence.size() > 1){
                pre_y = route.sequence[route.sequence.size()-2];
            }

            auto z = saving.pair2;
            auto succ_z = route.sequence[1];
            if (pair2_route.sequence.size() > 1){
                succ_z = route.sequence[1];
            }

            int best_modey = -1;
            int best_modez = -1;
            double best_saving = DBL_MIN;

            vector<int> y_possible = mcgrp.exclude(y);
            y_possible.push_back(y);
            vector<int> z_possible = mcgrp.exclude(z);
            z_possible.push_back(z);


            for(auto modey : y_possible){
                for (auto modez : z_possible){
                    if(mcgrp.inst_tasks[modey].time_window.first + mcgrp.inst_tasks[modey].serve_time + \
                            mcgrp.min_time[mcgrp.inst_tasks[modey].tail_node][mcgrp.inst_tasks[modez].head_node] <= mcgrp.inst_tasks[modez].time_window.second){

                        double saving_cost = d[pre_y][y] + d[y][0] + d[0][z] + d[z][succ_z] - \
                                      d[pre_y][modey] - d[modez][succ_z] - d[modey][modez];

                        if(saving_cost > best_saving){
                            best_saving = saving_cost;
                            best_modey = modey;
                            best_modez = modez;
                        }

                    }


                }
            }

            if (best_saving > 0){
                saving.value = best_saving;
                saving.pair1mode = best_modey;
                saving.pair2mode = best_modez;
            }

            update = true;

        }

        if(update == true){
            sort(saving_list.begin(),saving_list.end(),Cell::cmp);
        }
    }

    void initialize_savings_(){
        saving_list.clear();

        for (auto y : task_list) {
            for (auto z : task_list) {

                if (y == z || mcgrp.inst_tasks[y].demand + mcgrp.inst_tasks[z].demand > mcgrp.capacity) {
                    continue;
                }

                bool duplicate_flag = false;
                for( auto tmp : mcgrp.exclude(z)){
                    if(y == tmp)
                        duplicate_flag = true;
                }

                for( auto tmp : mcgrp.exclude(y)){
                    if(z == tmp)
                        duplicate_flag = true;
                }

                if(duplicate_flag) continue;

                int best_modey = -1;
                int best_modez = -1;
                double best_saving = -DBL_MAX;

                vector<int> y_possible = mcgrp.exclude(y);
                y_possible.push_back(y);
                vector<int> z_possible = mcgrp.exclude(z);
                z_possible.push_back(z);

                for (int modey : y_possible) {
                    for (int modez : z_possible) {
                        if (mcgrp.inst_tasks[modey].time_window.first + mcgrp.inst_tasks[modey].serve_time + \
                                t[modey][modez] <= mcgrp.inst_tasks[modez].time_window.second){

                            double saving = d[0][y] + d[y][0] + d[0][z] + d[z][0] - \
                                     d[0][modey] - d[modez][0] - d[modey][modez];

                            if(saving > best_saving){
                                best_saving = saving;
                                best_modey = modey;
                                best_modez = modez;
                            }

                        }
                    }
                }

                if (best_saving > 0){
                    saving_list.emplace_back(Cell(y, z, best_modey, best_modez, best_saving));
                }

            }
        }

        sort(saving_list.begin(),saving_list.end(),Cell::cmp);

    }

    void print_saving_list(){
        unordered_map<int,int> tbl{
            {1,1},
            {2,3},
            {3,5},
            {4,7},
            {5,9},
            {6,11},
            {7,13},
            {8,15},
            {9,17},
            {10,19},
            {11,21},
            {12,2},
            {13,4},
            {14,6},
            {15,8},
            {16,10},
            {17,12},
            {18,14},
            {19,16},
            {20,18},
            {21,20},
            {22,22},
        };


        for(const auto & item : saving_list){
            cout << "["<< tbl[item.pair1] << "," << tbl[item.pair2] <<"]"<< " " << item.value << endl;
        }
    }

    CWScanner(const MCGRPTW &mcgrp);
};

class NearestScanner: public PathConstructor{
private:
    Distance& distance;
public:
    NearestScanner(const MCGRPTW &mcgrp, Distance& distance_);

    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};

class RTFScanner: public PathConstructor{
private:
    Distance& distance;
public:
    RTFScanner(const MCGRPTW &mcgrp, Distance &distance_);

    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};

class SampleScanner: public PathConstructor{
private:
    LearningDistance& distance;
    int sample_times;
public:
    SampleScanner(const MCGRPTW &mcgrp, LearningDistance &distance_, int sample_times_ = 50);

    Individual
    operator()(const vector<int> &taskList = vector<int>(), const string& mode="feasible") override;
};


/*!
 * execute merge and split to the current solution
 * @param ns
 * @param mcgrp
 */
void merge_split(LocalSearch &ns, const MCGRPTW &mcgrp,
                 const int merge_size);

/*!
 * @details use different policies to merge Task
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> split_task(const MCGRPTW &mcgrp, LocalSearch& ns, const vector<int> &tasks, Policy& policy);

/*!
 * @details use nearest L2 distance policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> nearest_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use nearest L2 distance to depot policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> nearest_depot_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use max yield  policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> maximum_yield_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use min yield  policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> minimum_yield_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

/*!
 * @details use mixture policy to merge tasks
 * @param mcgrp
 * @param tasks
 * @return
 */
vector<int> mixture_growing(const MCGRPTW &mcgrp, vector<int> tasks, Policy& policy);

vector<vector<int>> tour_splitting(const MCGRPTW &mcgrp, vector<int>& task_list);


struct BestServeSeq{
    vector<int> sequence;
    vector<int> arrive_time_tbl;
    int cost = INT32_MAX;
};

BestServeSeq viterbi_decoding(const MCGRPTW &mcgrp, const vector<int>& task_list, bool allow_infeasible);

#endif
