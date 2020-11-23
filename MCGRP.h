#pragma once
#include "RNG.h"
#include "utils.h"
#include <bits/stdc++.h>
#include <mutex>

#define DUMMY 0
#define ARC_NO_INVERSE -1
#define NODE_NO_INVERSE -2


/*Global static info used by problem*/

class MCGRP
{
public:
    std::string instance_name;

    int node_num;
    int edge_num;
    int arc_num;

    int req_node_num;
    int req_edge_num;
    int req_arc_num;

    int nonreq_node_num;
    int nonreq_edge_num;
    int nonreq_arc_num;

    int actual_task_num;
    int total_arc_num;
    int vehicle_num;
    int capacity;
    int total_demand;

    int DEPOT;

    int neigh_size;

    std::vector<std::vector<int> > trav_cost;        //trave cost(without loading), trav_cost[节点数+1][节点数+1]

    std::vector<std::vector<int> > serve_cost;        //serve cost(with loading), serve_cost[节点数+1][节点数+1]

    std::vector<std::vector<int> > trave_time;
    std::vector<std::vector<int> > serve_time;


    std::vector<Arc> inst_arcs;

    //L2-distance (dijkstra)min_cost[node+1][node+1]
    std::vector<std::vector<int>> min_cost;
    std::vector<std::vector<int>> min_time;


    // shortest_path, 1st dimension: origin node; 2nd dimension: terminal node; 3rd dimension: sequence(the fisrt
    // entity represent the whole node number in the sequence
    std::vector<std::vector<std::vector<int> > > shortest_path;
    std::vector<Task> inst_tasks;
    int sentinel;


    vector<vector<double> > task_dist;
    vector<vector<NeighborInfo>> neighbor;

    double total_service_cost;
    RNG &_rng;

    mutable unordered_map<string, unique_ptr<class Distance>> distance_look_tbl;

    /* stage global sol information */
    mutable std::vector<int> best_sol_neg;
    mutable double best_total_route_length;
    mutable double best_sol_time;

    mutable mutex instance_mutex;


    MCGRP(const InstanceNumInfo &instance_info, RNG &rng);

    /*!
     * convert a delimiter format solution into individual style
     * @param seq
     * @return
     */
    Individual parse_delimiter_seq(const vector<int> &seq) const;

    inline double get_yield(int task_num) const
    {
        if (inst_tasks[task_num].inverse == NODE_NO_INVERSE) {
            return 1.0 * inst_tasks[task_num].demand / (inst_tasks[task_num].serv_cost + 0.001);
        }
        else {
            return 1.0 * inst_tasks[task_num].demand / inst_tasks[task_num].serv_cost;
        }
    }

    /*!
     * @details 检查当前解是否是一轮搜索中的全局最优解
     * @param total_route_length
     * @param sol_seq
     * @return
     */
    bool check_best_solution(const double total_route_length, const std::vector<int> &sol_seq) const;


    /*!
     * 重置算例信息
     * @param rng
     */
    void reset(RNG &rng);

    /*!
     * 读入文件信息，获得需求边信息数组，无需求边信息数组，遍历成本数组，服务成本数组，最短距离矩阵，总共待服务成本
     * @param input_file
     * @param instance_info
     */
    void load_file_info(std::string input_file, const InstanceNumInfo &instance_info);

    /*!
     * 使用dijkstra算法求解各点之间的最短路径及路径序列
     * This will capable for asymmetric graph
     */
    void dijkstra();
    void get_trave_matrix();

    void create_neighbor_lists(int neighbor_size);

private:
    void _build_neighbor_task(const Task& task, NeighborInfo& neighbor_info);
    void _build_neighbor_node(int start, int end, NeighborInfo& neighbor_info);
    unordered_map<int, vector<int>> end_task_lookup_tbl;
    unordered_map<int, vector<int>> start_task_lookup_tbl;
    const vector<int>& _same_end_task(int end_node);
    const vector<int>& _same_start_task(int start_node);

public:
    void register_distance(string name, unique_ptr<Distance> distance_);


    /*!
     * 获得任务序列的总成本
     * @param delimiter_seq
     * @return
     */
    int get_task_seq_total_cost(const std::vector<int> &delimiter_seq) const;

    /*!
     * 获得任务序列（一个解）的总超载量
     * @param route_seg_load
     * @return
     */
    int get_total_vio_load(const std::vector<int> &route_seg_load) const;

    /*!
     * check whether Task is edge Task
     * @param task_id
     * @return
     */
    inline bool is_edge(int task_id) const
    {
        return !(inst_tasks[task_id].inverse == ARC_NO_INVERSE ||
            inst_tasks[task_id].inverse == NODE_NO_INVERSE ||
            inst_tasks[task_id].inverse == DUMMY);
    }

    /*!
     * @details pack best sol info in mcgrp
     * @param p
     */
    void create_individual(Individual &p) const;

    /*!
     * test the correctness of the given solution
     * the solution is the negative coding style
     * which is a sequence of tasks
     * @param neg_seq
     * @param sol_cost
     * @return
     */
    bool valid_sol(const std::vector<int> &neg_seq, const double sol_cost) const;

    int get_travel_time(int source_task,int sink_task) const;

    vector<int> cal_arrive_time(const vector<int>& route) const;

    /*!
     *
     * @param source
     * @param sink
     * @param start
     * @param head_of_source true: need to sum the serve time of source Task, vice versa
     * @return
     */
    int cal_arrive_time(int source, int sink, int start, bool head_of_source) const;

    /*!
     * forecast the time table when modify the serve sequence
     * @param old_table
     * @param tasks
     * @param indicator_task
     * @param mode: 1. "insert_before": insert tasks in the route before a Task
     *              2. "insert_after": insert tasks in the route after a Task
     *              3. "remove": remove tasks in the route
     * @return:     {-1,-1} for infeasible time table
     */
    vector<RouteInfo::TimeTable> forecast_time_table(const vector<RouteInfo::TimeTable>& old_table,
                                                     const vector<int>& tasks,
                                                     string mode, int indicator_task = -1,
                                                     bool allow_infeasible = false) const;

    bool isTimeTableFeasible(const vector<RouteInfo::TimeTable>& tbl, bool meta = true) const;

    int get_vio_time(const vector<RouteInfo::TimeTable>& tbl) const;

    /*!
     * sort a set of tasks based on time
     * 1. the start of the time window
     * 2. the serve time
     * @param sequence
     */
    void time_window_sort(vector<int>& sequence) const;

    vector<vector<RouteInfo::TimeTable>> get_time_tbl(const vector<int>& sequence, string mode) const;

    int count_edges(const vector<int>& seq) const;

    void show_shortest_info(string output_file = "");
    void show_neighbor(int start, int end);
};

inline int MCGRP::get_travel_time(int source_task, int sink_task) const {
        return min_time[inst_tasks[source_task].tail_node][inst_tasks[sink_task].head_node];
}


// different Distance metric between tasks
class Distance{
protected:
    string name;

    vector<vector<double>> distance_matrix;
public:
    Distance(const MCGRP& mcgrp, const string &name);

    virtual double operator()(const MCGRP& mcgrp, const int task_a, const int task_b) = 0;
};

class CostDistance: public Distance{

public:
    CostDistance(const MCGRP& mcgrp);

    double operator()(const MCGRP& mcgrp, const int task_a, const int task_b) override;
};

class HybridDistance: public Distance{
private:
    double mean_cost;
    double deviation_cost;
    double mean_waiting_time;
    double deviation_waiting_time;

    vector<vector<double>> cost_matrix;
    vector<vector<double>> waiting_time_matrix;

    double beta;

public:
    HybridDistance(const MCGRP &mcgrp,double beta_);

    double operator()(const MCGRP& mcgrp, const int task_a, const int task_b) override;
};


class LearningDistance: public Distance{
private:
    vector<vector<double>> prob_matrix;

public:
    LearningDistance(const MCGRP &mcgrp, const vector<vector<double>> &probMatrix);

    double operator()(const MCGRP& mcgrp, const int task_a, const int task_b) override;

    const vector<double>& get_prob_vector(int task);
};