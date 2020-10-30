#pragma once
#include "RNG.h"
#include "utils.h"
#include <string>
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


    std::vector<arc> inst_arcs;

    //L2-distance (dijkstra)min_cost[node+1][node+1]
    std::vector<std::vector<int>> min_cost;
    std::vector<std::vector<int>> min_time;


    // shortest_path, 1st dimension: origin node; 2nd dimension: terminal node; 3rd dimension: sequence(the fisrt
    // entity represent the whole node number in the sequence
    std::vector<std::vector<std::vector<int> > > shortest_path;
    std::vector<task> inst_tasks;
    int sentinel;


    std::vector<std::vector<double> > task_dist;
    std::vector<std::vector<MCGRPNeighborInfo> > task_neigh_list;
    std::vector<std::vector<MCGRPNeighborInfo> > predecessor_task_neigh_list;
    std::vector<std::vector<MCGRPNeighborInfo> > successor_task_neigh_list;

    double total_service_cost;
    RNG &_rng;


    /* stage global sol information */
    mutable std::vector<int> best_sol_buff;     //negative coding format
    mutable double best_total_route_length;
    mutable double best_sol_time;

    /* global sol information */
    mutable std::vector<int> global_best_sol_buff;     //negative coding format
    mutable double global_best_total_route_length;
    mutable double global_best_sol_time;
    mutable mutex global_mut;


    MCGRP(const instance_num_information &instance_info, RNG &rng);

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
    void load_file_info(std::string input_file, const instance_num_information &instance_info);

    /*!
     * 使用dijkstra算法求解各点之间的最短路径及路径序列
     * This will capable for asymmetric graph
     */
    void dijkstra();
    void get_trave_matrix();

    /*!
     * 初始化各任务的邻域表（近该任务最近的几个任务）
     * @param neighbor_size
     */
    void create_neighbor_lists();

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
     * check whether task is edge task
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
     * @param head_of_source true: need to sum the serve time of source task, vice versa
     * @return
     */
    int cal_arrive_time(int source, int sink, int start, bool head_of_source) const;

    /*!
     * forecast the time table when modify the serve sequence
     * @param old_table
     * @param tasks
     * @param indicator_task
     * @param mode: 1. "insert_before": insert tasks in the route before a task
     *              2. "insert_after": insert tasks in the route after a task
     *              3. "remove": remove tasks in the route
     * @return:     {-1,-1} for infeasible time table
     */
    vector<MCGRPRoute::Timetable> forecast_time_table(const vector<MCGRPRoute::Timetable>& old_table,
                                                      const vector<int>& tasks,
                                                      string mode, int indicator_task = -1,
                                                      bool allow_infeasible = false) const;

    bool isTimetableFeasible(const vector<MCGRPRoute::Timetable>& tbl,bool meta = true) const;

    int get_vio_time(const vector<MCGRPRoute::Timetable>& tbl) const;

    /*!
     * sort a set of tasks based on time
     * 1. the start of the time window
     * 2. the serve time
     * @param sequence
     */
    void time_window_sort(vector<int>& sequence) const;

    vector<vector<MCGRPRoute::Timetable>> get_time_tbl(const vector<int>& sequence, string mode) const;

    int count_edges(const vector<int>& seq) const;
};

inline int MCGRP::get_travel_time(int source_task, int sink_task) const {
        return min_time[inst_tasks[source_task].tail_node][inst_tasks[sink_task].head_node];
}


