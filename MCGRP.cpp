#include "MCGRP.h"
#include "utils.h"
#include "json.hpp"


using namespace std;

using json = nlohmann::json;

extern timeb phase_start_time;

MCGRP::MCGRP(const InstanceNumInfo &instance_info, RNG &rng)
    : _rng(rng)
{
    //include dummy Task which id is 0
    int task_num = 2 * instance_info.req_edge_num + instance_info.req_arc_num + instance_info.req_node_num
        + 1;        //The number of the required tasks

    //the first node id is 1, add an offset.
    int offset_node_num = instance_info.node_num + 1;  //The number of the whole node

    best_sol_neg.clear();

    trav_cost.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; i++) {
        trav_cost[i].resize(offset_node_num);
    }

    serve_cost.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; ++i) {
        serve_cost[i].resize(offset_node_num);
    }

    trave_time.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; i++) {
        trave_time[i].resize(offset_node_num);
    }

    serve_time.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; ++i) {
        serve_time[i].resize(offset_node_num);
    }

    min_cost.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; ++i) {
        min_cost[i].resize(offset_node_num);
    }

    min_time.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; ++i) {
        min_time[i].resize(offset_node_num);
    }

    My_Assert(serve_cost[0].empty() && min_cost[0].empty() && trav_cost[0].empty(), "Dummy row should be null!");

    task_dist.resize(task_num);
    for (int i = 0; i < task_num; i++) {
        task_dist[i].resize(task_num);
    }

    shortest_path.resize(offset_node_num);
    for (int i = 1; i < offset_node_num; ++i) {
        shortest_path[i].resize(offset_node_num);
    }


    best_total_route_length = DBL_MAX;
    best_sol_time = DBL_MAX;
}

void MCGRP::create_individual(Individual &p) const
{
    My_Assert(!best_sol_neg.empty(), "Best solution cannot be empty!");

    p.sequence.clear();
    p.sequence = get_delimiter_coding(best_sol_neg);

    p.route_seg_load.clear();
    p.total_vio_load = 0;

    int load = 0;
    My_Assert(p.sequence.front() == DUMMY, "Incorrect delimiter coding!");

    double total_cost = 0;
    total_cost = inst_tasks[p.sequence.front()].serv_cost;
    for (auto cursor = 1; cursor < p.sequence.size(); cursor++) {
        total_cost += min_cost[inst_tasks[p.sequence[cursor - 1]].tail_node][inst_tasks[p.sequence[cursor]].head_node];
        total_cost += inst_tasks[p.sequence[cursor]].serv_cost;

        if (p.sequence[cursor] == DUMMY) {
            p.route_seg_load.push_back(load);
            p.total_vio_load += max((load - capacity), 0);
            load = 0;
        }
        else {
            load += inst_tasks[p.sequence[cursor]].demand;
        }
    }

    My_Assert(total_cost == best_total_route_length, "Wrong result of total route length!");

    p.total_cost = total_cost;
}

void MCGRP::reset(RNG &rng)
{
    _rng = rng;
    best_sol_neg.clear();
    best_total_route_length = DBL_MAX;
    best_sol_time = DBL_MAX;
}

void MCGRP::load_file_info(std::string input_file, const InstanceNumInfo &instance_info)
{
    std::ifstream fin(input_file);

    std::string dummy_string;

    node_num = instance_info.node_num;
    arc_num = instance_info.arc_num;
    edge_num = instance_info.edge_num;
    req_node_num = instance_info.req_node_num;
    req_arc_num = instance_info.req_arc_num;
    req_edge_num = instance_info.req_edge_num;
    nonreq_node_num = node_num - req_node_num;
    nonreq_arc_num = arc_num - req_arc_num;
    nonreq_edge_num = edge_num - req_edge_num;
    total_arc_num = 2 * req_edge_num + req_arc_num + 2 * nonreq_edge_num + nonreq_arc_num;


    actual_task_num = 2 * req_edge_num + req_arc_num + req_node_num;

    inst_tasks.clear();
    //one is dummy Task and the other is sentinel Task
    inst_tasks.resize(actual_task_num + 2);
    for(int ii = 0; ii < inst_tasks.size(); ii++){
        inst_tasks[ii].task_id = ii;
    }

    inst_arcs.clear();
    //one is dummy arcs
    inst_arcs.resize(total_arc_num + 1);

    total_demand = 0;
    while (fin >> dummy_string) {
        if (dummy_string == "Name:") {
            fin >> dummy_string;
            instance_name = dummy_string;
        }
        else if (dummy_string == "Optimal") {
            fin >> dummy_string;
            fin >> dummy_string;
            if (std::stod(dummy_string) != -1) {
                cout << "Best solution ever: " + dummy_string << endl;
            }
        }
        else if (dummy_string == "#Vehicles:") {
            fin >> dummy_string;
            vehicle_num = stoi(dummy_string);
        }
        else if (dummy_string == "Capacity:") {
            fin >> dummy_string;
            capacity = stoi(dummy_string);
        }
        else if (dummy_string == "Depot") {
            fin >> dummy_string;
            fin >> dummy_string;
            DEPOT = stoi(dummy_string);
        }
        else if (dummy_string == "ReE.") {
            getline(fin, dummy_string);
            for (int i = 1; i <= req_edge_num; i++) {
                fin >> dummy_string;
                inst_tasks[i].task_name = dummy_string;

                fin >> dummy_string;
                inst_tasks[i].head_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].tail_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].trave_cost = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].demand = stoi(dummy_string);
                total_demand += inst_tasks[i].demand;

                fin >> dummy_string;
                //inst_tasks[i].serv_cost = stoi(dummy_string);
                inst_tasks[i].serv_cost = inst_tasks[i].trave_cost;

                fin >> dummy_string;
                inst_tasks[i].trave_time = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].serve_time = stoi(dummy_string);

                Task::Window window;
                fin >> dummy_string;
                window.first = stoi(dummy_string);
                fin >> dummy_string;
                window.second = stoi(dummy_string);
                inst_tasks[i].time_window = window;

                inst_tasks[i].inverse = i + req_edge_num;

                inst_tasks[i + req_edge_num] = inst_tasks[i];
                inst_tasks[i + req_edge_num].task_id = i + req_edge_num;
                swap(inst_tasks[i + req_edge_num].head_node, inst_tasks[i + req_edge_num].tail_node);
                inst_tasks[i + req_edge_num].inverse = i;

                inst_arcs[i].head_node = inst_tasks[i].head_node;
                inst_arcs[i].tail_node = inst_tasks[i].tail_node;
                inst_arcs[i].trav_cost = inst_tasks[i].trave_cost;
                inst_arcs[i].trav_time = inst_tasks[i].trave_time;

                inst_arcs[i + req_edge_num] = inst_arcs[i];
                swap(inst_arcs[i + req_edge_num].head_node, inst_arcs[i + req_edge_num].tail_node);
            }
        }
        else if (dummy_string == "EDGE") {
            getline(fin, dummy_string);
            for (int i = 2 * req_edge_num + 1; i <= 2 * req_edge_num + nonreq_edge_num; i++) {
                fin >> dummy_string;

                fin >> dummy_string;
                inst_arcs[i].head_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_arcs[i].tail_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_arcs[i].trav_cost = stoi(dummy_string);

                fin >> dummy_string;
                inst_arcs[i].trav_time = stoi(dummy_string);

                inst_arcs[i + nonreq_edge_num] = inst_arcs[i];
                swap(inst_arcs[i + nonreq_edge_num].head_node, inst_arcs[i + nonreq_edge_num].tail_node);
            }
        }
        else if (dummy_string == "ReA.") {
            getline(fin, dummy_string);

            for (int i = 2 * req_edge_num + 1; i <= 2 * req_edge_num + req_arc_num; i++) {
                fin >> dummy_string;
                inst_tasks[i].task_name = dummy_string;

                fin >> dummy_string;
                inst_tasks[i].head_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].tail_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].trave_cost = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].demand = stoi(dummy_string);
                total_demand += inst_tasks[i].demand;

                fin >> dummy_string;
//                inst_tasks[i].serv_cost = stoi(dummy_string);
                inst_tasks[i].serv_cost = inst_tasks[i].trave_cost;

                fin >> dummy_string;
                inst_tasks[i].trave_time = stoi(dummy_string);

                fin >> dummy_string;
                inst_tasks[i].serve_time = stoi(dummy_string);

                Task::Window window;
                fin >> dummy_string;
                window.first = stoi(dummy_string);
                fin >> dummy_string;
                window.second = stoi(dummy_string);
                inst_tasks[i].time_window = window;

                inst_tasks[i].inverse = ARC_NO_INVERSE;
            }
            for (int i = 2 * req_edge_num + 2 * nonreq_edge_num + 1;
                 i <= 2 * req_edge_num + 2 * nonreq_edge_num + req_arc_num; i++) {
                inst_arcs[i].head_node = inst_tasks[i - 2 * nonreq_edge_num].head_node;
                inst_arcs[i].tail_node = inst_tasks[i - 2 * nonreq_edge_num].tail_node;
                inst_arcs[i].trav_cost = inst_tasks[i - 2 * nonreq_edge_num].trave_cost;
                inst_arcs[i].trav_time = inst_tasks[i - 2 * nonreq_edge_num].trave_time;
            }
        }
        else if (dummy_string == "ARC") {
            getline(fin, dummy_string);

            for (int i = 2 * req_edge_num + 2 * nonreq_edge_num + req_arc_num + 1;
                 i <= 2 * req_edge_num + 2 * nonreq_edge_num + req_arc_num + nonreq_arc_num; i++) {
                fin >> dummy_string;

                fin >> dummy_string;
                inst_arcs[i].head_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_arcs[i].tail_node = stoi(dummy_string);

                fin >> dummy_string;
                inst_arcs[i].trav_cost = stoi(dummy_string);

                fin >> dummy_string;
                inst_arcs[i].trav_time = stoi(dummy_string);
            }
        }
        else if (dummy_string == "ReN.") {
            getline(fin, dummy_string);

            for (int i = 2 * req_edge_num + req_arc_num + 1; i <= 2 * req_edge_num + req_arc_num + req_node_num; i++) {
                fin >> dummy_string;
                inst_tasks[i].task_name = dummy_string;

                dummy_string = dummy_string.substr(1);
                int node = stoi(dummy_string);

                inst_tasks[i].head_node = node;
                inst_tasks[i].tail_node = node;

                fin >> dummy_string;
                inst_tasks[i].demand = stoi(dummy_string);
                total_demand += inst_tasks[i].demand;

                fin >> dummy_string;
				inst_tasks[i].serv_cost = stoi(dummy_string);
                inst_tasks[i].trave_cost = 0;

                fin >> dummy_string;
                inst_tasks[i].serve_time = stoi(dummy_string);
                inst_tasks[i].trave_time = 0;

                Task::Window window;
                fin >> dummy_string;
                window.first = stoi(dummy_string);
                fin >> dummy_string;
                window.second = stoi(dummy_string);
                inst_tasks[i].time_window = window;

                inst_tasks[i].inverse = NODE_NO_INVERSE;
            }
        }
    }

    fin.close();

    /* dummy Task information */
    inst_tasks[DUMMY].task_name = "DUMMY";
    inst_tasks[DUMMY].tail_node = DEPOT;
    inst_tasks[DUMMY].head_node = DEPOT;
    inst_tasks[DUMMY].trave_cost = 0;
    inst_tasks[DUMMY].serv_cost = 0;
    inst_tasks[DUMMY].demand = 0;
    inst_tasks[DUMMY].inverse = DUMMY;
    inst_tasks[DUMMY].time_window = {0, INT32_MAX};

    /* sentinel Task info */
    sentinel = inst_tasks.size() - 1;
    inst_tasks[sentinel].task_name = "DUMMY";
    inst_tasks[sentinel].head_node = DEPOT;
    inst_tasks[sentinel].tail_node = DEPOT;
    inst_tasks[sentinel].trave_cost = 0;
    inst_tasks[sentinel].serv_cost = 0;
    inst_tasks[sentinel].demand = 0;
    inst_tasks[sentinel].inverse = DUMMY;
    inst_tasks[sentinel].time_window = {0, INT32_MAX};


    inst_arcs[DUMMY].head_node = DEPOT;
    inst_arcs[DUMMY].tail_node = DEPOT;
    inst_arcs[DUMMY].trav_cost = 0;


    for (int i = 1; i <= node_num; i++) {
        for (int j = 1; j <= node_num; j++) {
            trav_cost[i][j] = INT32_MAX;
            serve_cost[i][j] = 0;
            trave_time[i][j] = INT32_MAX;
            serve_time[i][j] = 0;
        }
    }

    My_Assert(total_arc_num == 2 * req_edge_num + 2 * nonreq_edge_num + req_arc_num + nonreq_arc_num,
              "Arc number is not right!");
    for (int i = 1; i <= total_arc_num; i++) {
        trav_cost[inst_arcs[i].head_node][inst_arcs[i].tail_node] = inst_arcs[i].trav_cost;
        trave_time[inst_arcs[i].head_node][inst_arcs[i].tail_node] = inst_arcs[i].trav_time;
    }

    My_Assert(actual_task_num == 2 * req_edge_num + req_arc_num + req_node_num,
              "tasks number is not right!");
    for (int i = 1; i <= actual_task_num; i++) {
        serve_cost[inst_tasks[i].head_node][inst_tasks[i].tail_node] = inst_tasks[i].serv_cost;
        serve_time[inst_tasks[i].head_node][inst_tasks[i].tail_node] = inst_tasks[i].serve_time;
    }

    dijkstra();
    get_trave_matrix();

    // pairwise distance matrix between tasks
    // dist(i,j) is the average distance between the vertices of the tasks
    for (int i = 0; i <= actual_task_num; i++) {
        for (int j = 0; j <= actual_task_num; j++) {
            if (i == j)
                task_dist[i][j] = 0.0;
            else
                task_dist[i][j] = double(min_cost[inst_tasks[i].tail_node][inst_tasks[j].head_node] +
                    min_cost[inst_tasks[j].tail_node][inst_tasks[i].head_node]) / 2.0;
        }
    }

    total_service_cost = 0;
    for (int i = 1; i <= req_edge_num; ++i)
        total_service_cost += inst_tasks[i].serv_cost;
    for (int i = 2 * req_edge_num + 1; i <= 2 * req_edge_num + req_arc_num + req_node_num; i++)
        total_service_cost += inst_tasks[i].serv_cost;

}

void MCGRP::dijkstra()
{
    for (int i = 1; i <= node_num; i++) {
        for (int j = 1; j <= node_num; j++) {

            shortest_path[i][j].push_back(i);    //start from i

            if (j == i) {
                min_cost[i][j] = 0;
            }
            else {
                min_cost[i][j] = INT32_MAX;
            }
        }
    }

//    show_shortest_info();

    std::vector<bool> mark(node_num + 1, false);
    std::vector<int> dist(node_num + 1, 0);
    std::vector<int> dist1(node_num + 1, 0);
    std::vector<int> nearest_neighbor;

    for (int i = 1; i <= node_num; i++) {// For each vertex
        mark[i] = true;
        //i: start node, j: terminal node, k: time step
        for (int j = 1; j <= node_num; j++) { //mark all nodes need to be reached
            if (j == i)
                continue;

            mark[j] = false;
            dist[j] = trav_cost[i][j];
            dist1[j] = dist[j];
        }

        for (int k = 1; k < node_num; k++) {    // there are vertex_num steps
            int minimum = INT32_MAX;
            nearest_neighbor.clear();

            for (int j = 1; j <= node_num; j++) {//check for the newly nearest node can be reached in this time step
                if (mark[j])    //already reach
                    continue;

                if (dist1[j] == INT32_MAX)    //cannot reach this step
                    continue;

                if (dist1[j] < minimum)
                    minimum = dist1[j];
            }

            if (minimum == INT32_MAX)    //node i in time step k cannot reach any node
                continue;

            for (int j = 1; j <= node_num; j++) {    //find out all the nearest neighbour nodes
                if (mark[j])
                    continue;

                if (dist1[j] == minimum) {
                    nearest_neighbor.push_back(j);
                }
            }

            int v = nearest_neighbor.front();    //choose the first neighbour j
            dist1[v] = INT32_MAX;
            mark[v] = true;

            if (shortest_path[i][v].back() != v) {
                shortest_path[i][v].push_back(v);
            }

            for (int j = 1; j <= node_num; j++) {
                if (mark[j])
                    continue;

                if (trav_cost[v][j] != INT32_MAX
                    && minimum + trav_cost[v][j] < dist[j]) {
                    //to check whether point v as a middle point can reach point j faster
                    dist[j] = minimum + trav_cost[v][j];
                    dist1[j] = minimum + trav_cost[v][j];
                    //	replace old route to j with new route which bypasses v
                    shortest_path[i][j] = shortest_path[i][v];
                }
            }

            for (int j = 1; j <= node_num; j++) {
                if (j == i)
                    continue;

                min_cost[i][j] = dist[j];
            }
        }
    }

    for (int i = 1; i <= node_num; i++) {    //find out the isolated point
        for (int j = 1; j <= node_num; j++) {
            if (shortest_path[i][j].size() == 1)
                shortest_path[i][j].clear();
        }
    }

//    show_shortest_info("shortest_path.json");
}

void MCGRP::create_neighbor_lists(int neighbor_size)
{
    neigh_size = neighbor_size;
    neighbor = vector<vector<NeighborInfo>>(node_num + 1, vector<NeighborInfo>(node_num + 1));

    // First build task neighbor
    for(const auto& task : inst_tasks){
        _build_neighbor_task(task, neighbor[task.head_node][task.tail_node]);
    }

    // Second build node neighbor
    for( int row = 1; row <= node_num; row++){
        for(int col = 1; col <= node_num; col++){
            _build_neighbor_node(row, col,neighbor[row][col]);
        }
    }
}

int MCGRP::get_task_seq_total_cost(const std::vector<int> &delimiter_seq) const
{
    int total_cost = 0;
    assert(delimiter_seq[delimiter_seq.size() - 1] == DUMMY);

    for (int i = 0; i < delimiter_seq.size() - 1; i++)
        total_cost += min_cost[inst_tasks[delimiter_seq[i]].tail_node][inst_tasks[delimiter_seq[i + 1]].head_node]
            + inst_tasks[delimiter_seq[i]].serv_cost;

    return total_cost;
}

int MCGRP::get_total_vio_load(const std::vector<int> &route_seg_load) const
{
    int total_vio_load = 0;
    for (auto load :route_seg_load) {
        if (load > capacity)
            total_vio_load += load - capacity;
    }

    return total_vio_load;
}

bool MCGRP::check_best_solution(const double total_route_length, const vector<int> &sol_seq) const
{
    lock_guard<mutex> lk(instance_mutex);
    // Determines if the current solution is the best found so far.
    if ((total_route_length < best_total_route_length)) {
        My_Assert(valid_sol(sol_seq, total_route_length), "Wrong solution!");

        best_total_route_length = total_route_length;
        best_sol_neg = sol_seq;
        struct timeb cur_time;
        ftime(&cur_time);
        best_sol_time =
            (cur_time.time - phase_start_time.time) + ((cur_time.millitm - phase_start_time.millitm) * 1.0 / 1000);

//        My_Assert(print(log_out, to_string((int) total_route_length)),"Wrong print");
#ifdef DEBUG
        print(log_out, to_string((int) total_route_length));
#endif
        return true;
    }
    else {

//        My_Assert(print(log_out, to_string((int) total_route_length)),"Wrong print");
#ifdef DEBUG
        print(log_out, to_string((int) total_route_length));
#endif
        return false;
    }
}

Individual MCGRP::parse_delimiter_seq(const vector<int> &seq) const
{
    My_Assert(!seq.empty() && seq.front() == DUMMY && seq.back() == DUMMY, "invalid seq");

    Individual buffer;

    buffer.sequence.clear();
    buffer.sequence = seq;

    buffer.route_seg_load.clear();
    buffer.total_vio_load = 0;

    int load = 0;
    double total_cost = 0;
    total_cost = inst_tasks[buffer.sequence.front()].serv_cost;

    for (auto cursor = 1; cursor < buffer.sequence.size(); cursor++) {
        total_cost +=
            min_cost[inst_tasks[buffer.sequence[cursor - 1]].tail_node][inst_tasks[buffer.sequence[cursor]].head_node];
        total_cost += inst_tasks[buffer.sequence[cursor]].serv_cost;

        if (buffer.sequence[cursor] == DUMMY) {
            buffer.route_seg_load.push_back(load);
            buffer.total_vio_load += max((load - capacity), 0);
            load = 0;
        }
        else {
            load += inst_tasks[buffer.sequence[cursor]].demand;
        }
    }

    buffer.total_cost = total_cost;

    buffer.time_tbl = get_time_tbl(seq,"delimiter");
    buffer.total_vio_time = 0;
    for(const auto & time_tbl : buffer.time_tbl)
        buffer.total_vio_time +=  get_vio_time(time_tbl);

    return buffer;
}

bool MCGRP::valid_sol(const vector<int> &neg_seq, const double sol_cost) const
{
    int valid_length = 0;
    int load = 0;
    int drive_time = 0;

    for (int j = 0; j < neg_seq.size() - 1; j++) {
        My_Assert(neg_seq[j] != 0, "solution can't has dummy Task!");
        if (neg_seq[j] < 0) {
            valid_length += min_cost[inst_tasks[DUMMY].tail_node][inst_tasks[-neg_seq[j]].head_node];
            load = inst_tasks[-neg_seq[j]].demand;

            drive_time = cal_arrive_time(DUMMY, -neg_seq[j],0, true);
            My_Assert(load <= capacity && drive_time <= inst_tasks[-neg_seq[j]].time_window.second,
                      "solution violate constraints!");
        }
        else {
            My_Assert(j != 0, "First Task must be negative!");
            valid_length += min_cost[inst_tasks[abs(neg_seq[j - 1])].tail_node][inst_tasks[neg_seq[j]].head_node];
            load += inst_tasks[neg_seq[j]].demand;
            drive_time = cal_arrive_time(abs(neg_seq[j - 1]),neg_seq[j],drive_time,true);

            My_Assert(load <= capacity && drive_time <= inst_tasks[neg_seq[j]].time_window.second,
                      "solution violate constraints!");
        }

        valid_length += inst_tasks[abs(neg_seq[j])].serv_cost;

        if (neg_seq[j + 1] < 0) {
            valid_length += min_cost[inst_tasks[abs(neg_seq[j])].tail_node][inst_tasks[DUMMY].head_node];
        }
    }

    int j = neg_seq.size() - 1;
    My_Assert(neg_seq[j] != 0, "solution can't has dummy Task!");

    if (neg_seq[j] < 0) {
        load = inst_tasks[abs(neg_seq[j])].demand;

        drive_time = cal_arrive_time(DUMMY,-neg_seq[j],0,true);

        valid_length += min_cost[inst_tasks[DUMMY].tail_node][inst_tasks[-neg_seq[j]].head_node];
        My_Assert(load <= capacity && drive_time <= inst_tasks[abs(neg_seq[j])].time_window.second,
                  "solution violate constraints!");
    }
    else {
        valid_length += min_cost[inst_tasks[abs(neg_seq[j - 1])].tail_node][inst_tasks[neg_seq[j]].head_node];

        load += inst_tasks[neg_seq[j]].demand;
        drive_time = cal_arrive_time(abs(neg_seq[j - 1]),neg_seq[j],drive_time,true);

        My_Assert(load <= capacity && drive_time <= inst_tasks[abs(neg_seq[j])].time_window.second,
                  "solution violate constraints!");
    }

    valid_length += inst_tasks[abs(neg_seq[j])].serv_cost;

    valid_length += min_cost[inst_tasks[abs(neg_seq[j])].tail_node][inst_tasks[DUMMY].head_node];

    return valid_length == sol_cost;
}

void MCGRP::get_trave_matrix()
{
    for (int i = 1; i <= node_num; i++) {
        for (int j = 1; j <= node_num; j++) {
            if (shortest_path[i][j].empty()) {
                min_time[i][j] = 0;
            }
            else {
                vector<int> start = shortest_path[i][j];
                vector<int> end = vector<int>(start.begin() + 1, start.end());
                start.pop_back();
                int sum = 0;
                for (int k = 0; k < start.size(); k++)
                    sum += trave_time[start[k]][end[k]];
                My_Assert(sum % min_cost[i][j] == 0,
                          "shortest path cost and shortest path time have wrong relationship");
                min_time[i][j] = sum;
            }
        }
    }

}

vector<int> MCGRP::cal_arrive_time(const vector<int> &route) const
{
    if (route.empty())
        return vector<int>();

    vector<int> arrive_time;
    My_Assert(route.front() != DUMMY && route.back() != DUMMY,
              "invalid route!");

    int drive_time = cal_arrive_time(DUMMY, route.front(),0,true);
    arrive_time.push_back(drive_time);
    for (int i = 1; i < route.size(); i++) {
        drive_time = cal_arrive_time(route[i-1],route[i],drive_time,true);
        arrive_time.push_back(drive_time);
    }

    return arrive_time;
}

int MCGRP::cal_arrive_time(int source, int sink, int start, bool head_of_source) const
{
    int t1 = 0;
    if(head_of_source){
        t1 = start + inst_tasks[source].serve_time
            + min_time[inst_tasks[source].tail_node][inst_tasks[sink].head_node];
    }else{
        t1 = start + min_time[inst_tasks[source].tail_node][inst_tasks[sink].head_node];
    }

    return max(inst_tasks[sink].time_window.first, t1);
}

vector<RouteInfo::TimeTable> MCGRP::forecast_time_table(const vector<RouteInfo::TimeTable> &old_table,
                                                        const vector<int> &tasks,
                                                        string mode,
                                                        int indicator_task,
                                                        bool allow_infeasible) const
{

    if (mode == "remove") {
        if (old_table.empty())
            return vector<RouteInfo::TimeTable>({{-1, -1}});

        vector<RouteInfo::TimeTable> buffer{{DUMMY, 0}};
        int cur = 0;
        for (; cur < old_table.size(); cur++) {
            if (old_table[cur].task != tasks.front())
                buffer.push_back(old_table[cur]);
            else
                break;
        }

        My_Assert(cur != old_table.size(), "removed Task doesn't exists!");
        cur += tasks.size();
        for (; cur < old_table.size(); cur++) {
            int drive_time = cal_arrive_time(buffer.back().task, old_table[cur].task,
                                             buffer.back().arrive_time, true);

            if (!allow_infeasible && drive_time > inst_tasks[old_table[cur].task].time_window.second)
                return vector<RouteInfo::TimeTable>({{-1, -1}});
            buffer.push_back(
                {old_table[cur].task, drive_time});
        }

        buffer.erase(buffer.begin());

        return buffer;
    }
    else if (mode.find("insert") != string::npos) {
        My_Assert(indicator_task != -1, "Wrong indicator!");
        vector<RouteInfo::TimeTable> buffer{{DUMMY, 0}};
        int cur = 0;
        for (; cur < old_table.size(); cur++) {
            if (old_table[cur].task != indicator_task)
                buffer.push_back(old_table[cur]);
            else
                break;
        }
        My_Assert(cur != old_table.size(), "indicator Task doesn't exists!");

        if (mode == "insert_before") {
            // pass
        }
        else if (mode == "insert_after") {
            buffer.push_back(old_table[cur++]);
        }
        else {
            My_Assert(false, "Unknown parameters!");
        }

        for (int i = 0; i < tasks.size(); i++) {
            int drive_time = cal_arrive_time(buffer.back().task,tasks[i],
                                             buffer.back().arrive_time,true);
            if (!allow_infeasible && drive_time > inst_tasks[tasks[i]].time_window.second)
                return vector<RouteInfo::TimeTable>({{-1, -1}});
            buffer.push_back({tasks[i], drive_time});
        }

        for (; cur < old_table.size(); cur++) {
            int drive_time = cal_arrive_time(buffer.back().task, old_table[cur].task, buffer.back().arrive_time,true);
            if (!allow_infeasible && drive_time > inst_tasks[old_table[cur].task].time_window.second)
                return vector<RouteInfo::TimeTable>({{-1, -1}});
            buffer.push_back(
                {old_table[cur].task, max(drive_time, inst_tasks[old_table[cur].task].time_window.first)});
        }

        buffer.erase(buffer.begin());

        return buffer;
    }
    else {
        My_Assert(false, "Unknown parameters!");
        return vector<RouteInfo::TimeTable>();
    }

}

int MCGRP::get_vio_time(const vector<RouteInfo::TimeTable> &tbl) const
{
    int ans = 0;
    for (const auto node : tbl)
        ans += max(node.arrive_time - inst_tasks[node.task].time_window.second, 0);
    return ans;
}

void MCGRP::time_window_sort(vector<int> &sequence) const
{
    function<bool(const int,const int)> comp = [this](const int a,const int b){
        if(this->inst_tasks[a].time_window.second == inst_tasks[b].time_window.second)
            return this->inst_tasks[a].serve_time < inst_tasks[b].serve_time;
        return this->inst_tasks[a].time_window.second < inst_tasks[b].time_window.second;
    };

    function<int(int,int)> partition = [&sequence,&comp](int low, int high){
        int pivot = sequence[high]; // pivot
        int i = (low - 1); // Index of smaller element

        for (int j = low; j <= high - 1; j++){
            // If current element is smaller than the pivot
            if (comp(j,high)){
                i++; // increment index of smaller element
                swap(sequence[i], sequence[j]);
            }
        }

        swap(sequence[i + 1], sequence[high]);
        return (i + 1);
    };

    function<void(int,int)> quickSort = [&quickSort, &partition](int low, int high){
        if (low < high)
        {
            /* pi is partitioning index, arr[p] is now at right place */
            int pi = partition(low, high);

            // Separately sort elements before
            // partition and after partition
            quickSort(low, pi - 1);
            quickSort(pi + 1, high);
        }
    };

    quickSort(0,(int)sequence.size() - 1);
}

bool MCGRP::isTimeTableFeasible(const vector<RouteInfo::TimeTable>& tbl, bool meta) const {
    if(meta)
        return tbl.empty() || tbl.front().task != -1;
    else{
        for(const auto & node : tbl){
            if(node.arrive_time > inst_tasks[node.task].time_window.second)
                return false;
        }
        return true;
    }
}

vector<vector<RouteInfo::TimeTable>> MCGRP::get_time_tbl(const vector<int> &sequence, string mode) const
{
    vector<vector<int>> route_tasks;
    if(mode == "delimiter"){
        for(int cursor = 0;cursor<sequence.size() - 1;cursor++){
            if(sequence[cursor] == DUMMY)
                route_tasks.push_back(vector<int>());
            else
                route_tasks.back().push_back(sequence[cursor]);
        }
    }else if(mode == "negative"){
        for(const auto task : sequence){
            if(task < 0)
                route_tasks.push_back({-task});
            else
                route_tasks.back().push_back(task);
        }
    }else{
        cerr << "Invalid mode\n";
        abort();
    }

    vector<vector<RouteInfo::TimeTable>> res;

    for(const auto & cur_route : route_tasks){
        res.push_back(vector<RouteInfo::TimeTable>());
        auto time_tbl = cal_arrive_time(cur_route);
        for(int i = 0;i<time_tbl.size();i++)
            res.back().push_back({cur_route[i],time_tbl[i]});
    }

    return res;
}

int MCGRP::count_edges(const vector<int> &seq) const
{
    int ans = 0;
    for(const auto & task : seq)
        if(is_edge(task)) ans++;
    return ans;
}

void MCGRP::show_shortest_info(string output_file)
{
    if(!output_file.empty()){
        json js;
        int columnIndex = 0;

        auto min_cost_bak = min_cost;
        min_cost_bak.erase(min_cost_bak.begin());
        for_each(min_cost_bak.begin(), min_cost_bak.end(),
                 [&](vector<int>& row){
                     row.erase(next(row.begin(),columnIndex));
                 });
        js["distance"] = min_cost_bak;


        auto shortest_path_bak = shortest_path;
        shortest_path_bak.erase(shortest_path_bak.begin());
        for_each(shortest_path_bak.begin(), shortest_path_bak.end(),
                 [&](vector<vector<int>>& row){
                     row.erase(next(row.begin(),columnIndex));
                 });
        js["path"] = shortest_path_bak;

        ofstream fout(output_file);
        fout << setw(4) << js << endl;
        fout.close();
    }

    cout << "MIN Cost Matrix:\n\n\t";
    for (int i = 1; i <= node_num; i++) cout << i << '\t';
    cout << endl << endl;
    for (int i = 1; i <= node_num; i++) {
        cout << i << '\t';
        for (int j = 1; j <= node_num; j++) {
            if (min_cost[i][j] == INT32_MAX)
                cout << "INF\t";
            else cout << min_cost[i][j] << '\t';
        }
        cout << endl << endl;
    }

    cout << "Shortest Path Information:\n\n";
    for (int i = 1; i <= node_num; i++) {
        for (int j = 1; j <= node_num; j++) {
            if (min_cost[i][j] == INT32_MAX) {
                cout << i << "->" << j << ": no path\n";
                continue;
            }

            cout << i << "->" << j << "(through " << shortest_path[i][j].size() << " nodes): ";
            for (int k = 0; k < shortest_path[i][j].size(); k++) {
                cout << shortest_path[i][j][k] << "->";
            }
            cout << "\b\b  " << endl;
        }
        cout << endl << endl;
    }
}

void MCGRP::_build_neighbor_task(const Task &task, NeighborInfo &neighbor_info)
{
    if(neighbor_info.start != -1) return;

    neighbor_info.start = task.head_node;
    neighbor_info.end = task.tail_node;


    std::vector<TaskNeighborInfo> predecessor_NList;   //tmp neighborhood list
    std::vector<TaskNeighborInfo> successor_NList;   //tmp neighborhood list

    for(int neighbor_id = 0; neighbor_id <= actual_task_num; neighbor_id++){
        if(neighbor_id == task.task_id)
            continue;

        if(is_edge(task.task_id) && task.inverse == neighbor_id)
            continue;

        neighbor_info.basic_neighbor.push_back(TaskNeighborInfo(
                neighbor_id, task_dist[task.task_id][neighbor_id]));
    }

    sort(neighbor_info.basic_neighbor.begin(), neighbor_info.basic_neighbor.end(), TaskNeighborInfo::cmp);

    // TODO(luke): maybe I should only keep one of the edge task?

    for(const auto& neighbor: neighbor_info.basic_neighbor){
        if(inst_tasks[task.task_id].time_window.first
            + inst_tasks[task.task_id].serve_time
            + get_travel_time(task.task_id, neighbor.task_id)
            <= inst_tasks[neighbor.task_id].time_window.second){
            successor_NList.push_back(neighbor);
        }

        if(inst_tasks[neighbor.task_id].time_window.first
            + inst_tasks[neighbor.task_id].serve_time
            + get_travel_time(neighbor.task_id, task.task_id)
            <= inst_tasks[task.task_id].time_window.second){
            predecessor_NList.push_back(neighbor);
        }
    }

    neighbor_info.predecessor_neighbor.insert({task.task_id,predecessor_NList});
    neighbor_info.successor_neighbor.insert({task.task_id,successor_NList});

//    show_neighbor(neighbor_info.start,neighbor_info.end);
}


void MCGRP::_build_neighbor_node(int start, int end, NeighborInfo &neighbor_info)
{
    if(neighbor_info.start == -1){
        neighbor_info.start = start;
        neighbor_info.end = end;

        for(int neighbor_id = 0; neighbor_id <= actual_task_num; neighbor_id++){
            double distance = double(min_cost[end][inst_tasks[neighbor_id].head_node]
                                         + min_cost[inst_tasks[neighbor_id].tail_node][start]) / 2.0;

            neighbor_info.basic_neighbor.push_back(TaskNeighborInfo(
                neighbor_id, distance));
        }

        sort(neighbor_info.basic_neighbor.begin(), neighbor_info.basic_neighbor.end(), TaskNeighborInfo::cmp);
    }


    const vector<int> & candidate_end_tasks = _same_end_task(end);
    for(const auto task : candidate_end_tasks){
        vector<TaskNeighborInfo> successor_NList;

        for(const auto& neighbor: neighbor_info.basic_neighbor){
            if(task == neighbor.task_id) continue;

            if(is_edge(task) && inst_tasks[task].inverse == neighbor.task_id) continue;

            if(inst_tasks[task].time_window.first
                + inst_tasks[task].serve_time
                + get_travel_time(task, neighbor.task_id)
                <= inst_tasks[neighbor.task_id].time_window.second){
                successor_NList.push_back(neighbor);
            }
        }

        neighbor_info.successor_neighbor.insert({task,successor_NList});
    }

    const vector<int> & candidate_start_tasks = _same_start_task(start);
    for(const auto task : candidate_start_tasks){
        vector<TaskNeighborInfo> predecessor_NList;

        for(const auto& neighbor: neighbor_info.basic_neighbor){
            if(task == neighbor.task_id) continue;

            if(is_edge(task) && inst_tasks[task].inverse == neighbor.task_id) continue;

            if(inst_tasks[neighbor.task_id].time_window.first
                + inst_tasks[neighbor.task_id].serve_time
                + get_travel_time(neighbor.task_id, task)
                <= inst_tasks[task].time_window.second){
                predecessor_NList.push_back(neighbor);
            }
        }

        neighbor_info.predecessor_neighbor.insert({task,predecessor_NList});
    }


//    show_neighbor(start,end);
}


void MCGRP::show_neighbor(int start, int end)
{
    const NeighborInfo& task = neighbor[start][end];

    if(task.start == -1){
        cout << "Neighbor information from " << start << " to " << end << " Undefined!\n";
    }

    cout << "Neighbor information from " << start << " to " << end<< ":"<<endl;
    cout << "Basic Neighbor:\n";

    for(const auto& neighbor_task : task.basic_neighbor){
        cout << "[" << neighbor_task.task_id << "," << neighbor_task.distance << "]->";
    }
    cout <<"\b\b\n";


    for(const auto& pair: task.successor_neighbor){
        int task_id  = pair.first;
        cout << "task " << task_id << " successor task:"<<endl;
        for(const auto& neighbor_task : pair.second){
            cout << "[" << neighbor_task.task_id << "," << neighbor_task.distance << "]->";
        }
        cout <<"\b\b\n\n"<<endl;
    }

    for(const auto& pair: task.predecessor_neighbor){
        int task_id  = pair.first;
        cout << "task " << task_id << " predecessor task:"<<endl;
        for(const auto& neighbor_task : task.predecessor_neighbor.at(task_id)){
            cout << "[" << neighbor_task.task_id << "," << neighbor_task.distance << "]->";
        }
        cout <<"\b\b\n\n"<<endl;
    }

}

const vector<int>& MCGRP::_same_end_task(int end_node)
{
    if(!end_task_lookup_tbl[end_node].empty()){
        return end_task_lookup_tbl[end_node];
    }

    for(int ii = 0; ii<=actual_task_num; ii++){
        if(inst_tasks[ii].tail_node == end_node)
            end_task_lookup_tbl[end_node].push_back(ii);
    }
    return end_task_lookup_tbl[end_node];
}

const vector<int> &MCGRP::_same_start_task(int start_node)
{
    if(!start_task_lookup_tbl[start_node].empty()){
        return start_task_lookup_tbl[start_node];
    }

    for(int ii = 0; ii<=actual_task_num; ii++){
        if(inst_tasks[ii].head_node == start_node)
            start_task_lookup_tbl[start_node].push_back(ii);
    }
    return start_task_lookup_tbl[start_node];
}
