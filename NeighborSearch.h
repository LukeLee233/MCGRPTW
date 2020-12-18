#ifndef _NEIGHBORSEARCH_
#define _NEIGHBORSEARCH_

#include "RNG.h"
#include "utils.h"
#include "MCGRP.h"
#include "SearchPolicy.h"
#include "ConstructPolicy.h"
#include "config.h"
#include "operator.h"

extern vector<double> ratios;
extern vector<double> prob;


class HighSpeedNeighBorSearch
{
public:
    double best_solution_cost;
    vector<int> best_solution_neg;


    struct TASK_NODE{
        int ID;
        int route_id;   //tell the Task which route it belongs to
        TASK_NODE *pre;
        TASK_NODE *next;

        TASK_NODE(int id = 0){
            ID = id;
            route_id = -1;
            pre = nullptr;
            next = nullptr;
        }

        void clear(){
            route_id = -1;
            pre = nullptr;
            next = nullptr;
        }

    };

    struct DUMMYPOOL{
        //This class is used to manage dummy Task
        struct DUMMY_TASK{
            bool used;
            TASK_NODE task_info;

            DUMMY_TASK() : used(false){};
        };

        stack<int> unused;
        vector<DUMMY_TASK> pool;

        void extend();

        DUMMYPOOL(int initial_size = 10):pool(initial_size){
            for(int i = pool.size() - 1;i > 0; i--){
                unused.push(-i);
                pool[i].task_info.ID = -i;
            }
        }

        TASK_NODE* get_new_dummy(){
            if(unused.empty()){
//                extend();
                cerr<<"Bad new dummy Task!";
                abort();
            }

            int new_dummy_id = unused.top();
            unused.pop();
            My_Assert(new_dummy_id < 0 && !pool[-new_dummy_id].used,"Bad allocate for new dummy Task");
            pool[-new_dummy_id].used = true;

            return &pool[-new_dummy_id].task_info;
        }

        void free_dummy(int id){
            My_Assert(id < 0,"Bad free dummy Task");

            id = -id;

            if(pool[id].used == false){
                cerr<<"Double free!";
                abort();
            }

            pool[id].used = false;
            pool[id].task_info.clear();

            unused.push(pool[id].task_info.ID);
            return;
        }

        void clear(){
            for(int i = 0;i<pool.size();i++){
                if(pool[i].used){
                    pool[i].used = false;
                    unused.push(pool[i].task_info.ID);
                    pool[i].task_info.clear();
                }
            }

        }

    };

    struct SOLUTION{
        DUMMYPOOL dummypool;

        //for better memory management, use vector to implement double linked list
        vector<TASK_NODE> tasks;

        TASK_NODE *very_start;
        TASK_NODE *very_end;


        SOLUTION(unsigned int task_num): tasks(task_num + 1), dummypool(task_num+1){
            for(int id = 1;id<=task_num;id++){
                tasks[id].ID = id;
            }

            very_start = nullptr;
            very_end = nullptr;
        }

        TASK_NODE* operator[](int id)
        {
            //use ID to locate address
            My_Assert(id != 0,"Wrong ID");
            if(id > 0){
                return &tasks[id];
            }
            else{
                id = -id;
                if(dummypool.pool[id].used == false)
                {
                    cerr << "This dummy Task info is unavailable!" <<endl;
                    abort();
                }
                return &dummypool.pool[id].task_info;
            }

        }

        void print(){
            if(very_start == nullptr){
                return;
            }

            TASK_NODE * tmp = very_start;
            do{
                cout<<max (tmp->ID, 0) <<'-';
                tmp = tmp->next;
            }while(tmp != very_end);
            cout<<max (tmp->ID, 0)<<endl;
        }

        void clear(){
            dummypool.clear();
            very_start = nullptr;
            very_end = nullptr;
            for(int i = 0;i<tasks.size();i++){
                tasks[i].clear();
            }
        }

        pair<struct TASK_NODE*,struct TASK_NODE*>
        connect_tasks(const vector<int>& route_seq);

    };

    class ROUTESPOOL{
        struct ROUTE_NODE{
            bool used;
            RouteInfo route_info;

            ROUTE_NODE() : used(false){};
        };

        stack<int> unused;
        vector<ROUTE_NODE> pool;

    public:
        unordered_set<int> activated_route_id;
        ROUTESPOOL(int initial_size = 20):pool(initial_size){
            for(int i = pool.size() - 1; i>=0 ;i--){
                unused.push(i);
                pool[i].route_info.ID = i;
            }
        }

        RouteInfo* operator[](int id)
        {
            if(pool[id].used == false)
            {
                cerr << "This route info is unavailable!" <<endl;
                abort();
            }

            return &pool[id].route_info;
        }

        int allocate_route(){
            if(unused.empty()){
                cerr<<"Bad allocate!";
                abort();
            }

            int id = unused.top();
            unused.pop();
            pool[id].used = true;
            activated_route_id.insert(id);
            return id;
        }

        void free_route(int id){
            if(pool[id].used == false){
                cerr<<"Double free!";
                abort();
            }

            pool[id].used = false;
            pool[id].route_info.clear();

            unused.push(id);
            auto succeed = activated_route_id.erase(id);
            My_Assert(succeed == 1,"Bad Free");
            return;
        }

        void clear(){
            for(auto i : activated_route_id){
                My_Assert(pool[i].used == true,"Wrong routes!");
                pool[i].used = false;
                pool[i].route_info.clear();
                unused.push(i);
            }

            activated_route_id.clear();
        };

    };

    Policy policy;


    const int local_threshold = local_minimum_threshold;


    //dynamic info
    int neigh_size = 0;            //different policy use different neighbor size
    std::vector<int> search_space;
    /* current solution */
    SOLUTION solution;
    ROUTESPOOL routes;
    double cur_solution_cost;

    vector<int> task_set;

    double total_vio_load;
    int total_vio_time;

    //move operator
    unique_ptr<class XPreInsert> single_pre_insert;
    unique_ptr<class XPreInsert> double_pre_insert;
    unique_ptr<class XPreInsert> tribe_pre_insert;
    unique_ptr<class XPostInsert> single_post_insert;
    unique_ptr<class XPostInsert> double_post_insert;
    unique_ptr<class XPostInsert> tribe_post_insert;
    unique_ptr<class Invert> invert;
    unique_ptr<class NewSwap> swap;
    unique_ptr<class NewTwoOpt> two_opt;
    unique_ptr<class Extraction> extraction;
    unique_ptr<class Slice> slice;
    unique_ptr<class Attraction> attraction;

    // These variables are used to trace the search procedure
    string _stage;
    double _region_best;
    vector<int> _region_best_solution_neg;
    double _search_route_start;
    double _search_route_best;
    vector<int> _search_route_best_solution_neg;



public:

    HighSpeedNeighBorSearch(const MCGRP &mcgrp,int tabu_step_);
    ~HighSpeedNeighBorSearch();

    vector<vector<double>> score_matrix;
    vector<vector<double>> prob_matrix;
    void initialize_score_matrix(const MCGRP &mcgrp);
    void update_prob_matrix(vector<double>(*pf)(const vector<double>&));
    void print_prob_matrix(const string& filename="");
    void print_score_matrix(const string& filename="");

    vector<int> get_current_sol(string mode= "dummy");

    void delete_route(int route_id,const vector<int>& route_seq);
    int new_route(const MCGRP& mcgrp, const vector<int>& route_seq);

    void clear();

    inline vector<int> get_tasks_set(){return task_set;}

    inline double get_cur_cost(){
        return cur_solution_cost;
    }

    vector<int> get_sub_seq(int start, int end);

    /*!
     * sub procedure in decode_seq method

    /*!
     * @details 根据邻域表创建任务的邻域搜索空间
     * @param mcgrp
     * @param task
     */
    void create_search_neighborhood(const MCGRP &mcgrp, const vector<int>& chosen_seq,string mode = "", int offset = 0);

    /*!
     * @details pack current sol info to delimiter coding format
     * @param p
     */
    void dump_to_individual(const MCGRP &mcgrp, Individual &p);

    /*!
 * @details unpack sol info to negative coding format
 * @param dummy_seq
 */
    void unpack_seq(const std::vector<int> &dummy_seq, const MCGRP &mcgrp);

    /*!
     * @details neighbor search
     * @param mcgrp
     * @param max non improve_times
     * @param accept [First accept | Bext accept]
     */
    void neighbor_search(const MCGRP &mcgrp);

    /*!
     * @details downhill search
     * @param mcgrp
     */
    void descent_search(const MCGRP &mcgrp);

    /*!
     * @details random threshold search
     * @param mcgrp
     */
    void feasible_search(const MCGRP &mcgrp);

    /*!
     *  @details threshold exploration
     * @param mcgrp
     */
    void threshold_exploration(const MCGRP &mcgrp);

    /*!
     * @details descent exploration
     * @param mcgrp
     */
    void descent_exploration(const MCGRP &mcgrp);

    /*!
     * @details Infeasible search
     * @param mcgrp
     * @param start_point
     */
    void infeasible_search(const MCGRP &mcgrp);

    void small_step_infeasible_descent_exploration(const MCGRP &mcgrp);
    void small_step_infeasible_tabu_exploration(const MCGRP &mcgrp);

    void large_step_infeasible_exploration(const MCGRP &mcgrp);

    void update(const MCGRP &mcgrp, const vector<int>& best_buffer,const vector<int>& best_routes);

    /*!
     * @details repair overloaded solution
     * @param mcgrp
     */
    void repair_solution(const MCGRP &mcgrp);
    void _two_phase_repair(const MCGRP& mcgrp);
    void _repair_load(const MCGRP &mcgrp);
    void _repair_time_window(const MCGRP &mcgrp);
    void _tour_splitting_repair(const MCGRP &mcgrp);
    void _inplace_repair(const MCGRP &mcgrp);

    /*!
     * @details Check if any tasks missed.
     */
    bool missed(const MCGRP &mcgrp);

    /*!
     * @details check duplicate tasks
     * @param mcgrp
     * @return
     */
    bool check_duplicated(const MCGRP &mcgrp);

    bool valid_sol(const MCGRP& mcgrp);

    /*!
     * extract the searched solution to MCGRP instance
     * @param mcgrp: the implementation of a problem
     */
    void trace(const MCGRP &mcgrp);

    /*!
     * the implementation of the neighborhood search given a specific mode
     * @param mcgrp: the implementation of a problem
     * @param mode: search mode
     */
    void _neigh_search(const MCGRP &mcgrp, int mode);

    /*!
 * @details get the successor of the chosen task within the same route,exclude dummy task
 * @param ns
 * @param chosen_task
 * @return
 */
    vector<int> get_successor_tasks(int length, const int chosen_task);

    bool before(const int a,const int b);

    int get_time_window_violated_number(const MCGRP& mcgrp);
    double distance_to_feasible(const MCGRP& mcgrp);
    double _distance_to_feasible(const MCGRP& mcgrp, double vioc, double viot);
    double getFitness(const MCGRP& mcgrp,  Policy& policy);
    double getFitness(const MCGRP& mcgrp,  Policy& policy,const Individual& indi);
    double getFitnessDelta(const MCGRP& mcgrp, Policy& policy, const MoveResult &move_result);

    bool compress_empty_route(const int route_id);

    vector<RouteStable> get_route_stables(const MCGRP& mcgrp);
    void clear_stability(const MCGRP& mcgrp,const vector<int>& task_set);

#ifdef DEBUG
    int success_viterbi = 0;
#endif

    void viterbi_refine(const MCGRP &mcgrp);
};

#endif