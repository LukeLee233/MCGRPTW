#include "NeighborSearch.h"
#include <iostream>
#include "Swap.h"
#include "Invert.h"
#include "TwoOpt.h"
#include "Slice.h"
#include "Extraction.h"
#include <sys/timeb.h>
#include <vector>
#include <numeric>
#include <bitset>
#include <cmath>
#include <algorithm>
#include "config.h"
#include "utils.h"
#include "ConstructPolicy.h"

using namespace std;

vector<double> ratios = {0.07, 0.08, 0.09, 0.10};
//vector<double> ratios = {0.003,0.004,0.005,0.006};


vector<double> prob = {0.25, 0.25, 0.25, 0.25};


void HighSpeedNeighBorSearch::DUMMYPOOL::extend()
{
    unsigned int new_size = pool.size() * 2;
    unsigned int old_size = pool.size();

    //2^n growth
//    vector<DUMMY_TASK> pool_buffer;
//
//    std::move(std::begin(pool), std::end(pool), std::back_inserter(pool_buffer));
//    My_Assert(pool_buffer.size() == new_size / 2,"Incorrect move!");
    pool.resize(new_size);

//    pool = pool_buffer;

    //refill unused
    for (auto i = old_size; i < new_size; i++) {
        unused.push(i);
    }

    return;
}

HighSpeedNeighBorSearch::HighSpeedNeighBorSearch(const MCGRP &mcgrp)
    : solution(mcgrp.actual_task_num), routes(mcgrp.actual_task_num),
    task_set(mcgrp.actual_task_num), single_pre_insert(new XPreInsert(1)),
    single_post_insert(new XPostInsert(1)), double_pre_insert(new XPreInsert(2)),
    double_post_insert(new XPostInsert(2)),two_opt(new NewTwoOpt),
    invert(new Invert), swap(new NewSwap), extraction(new Extraction), slice(new Slice) {
    search_step = 0;
    equal_step = 0;
    cur_solution_cost = numeric_limits<decltype(best_solution_cost)>::max();
    total_vio_load = 0;
    total_vio_time = 0;
    best_solution_cost = numeric_limits<decltype(best_solution_cost)>::max();
    std::generate(task_set.begin(), task_set.end(), Generator());
    solution.print();
}

HighSpeedNeighBorSearch::~HighSpeedNeighBorSearch() = default;

void HighSpeedNeighBorSearch::neighbor_search(const MCGRP &mcgrp)
{
    int mode;
    if (neighbor_search_mode == "IDPRTR")
        mode = IDPRTR;
    else if (neighbor_search_mode == "RTRIDP")
        mode = RTRIDP;
    else
        mode = (int) mcgrp._rng.Randint(0, 1);

    _neigh_search(mcgrp, mode);

    trace(mcgrp);
}

void HighSpeedNeighBorSearch::_neigh_search(const MCGRP &mcgrp, int mode)
{
    if (mode == RTRIDP) {
        cout << "RTTR -> IDP\n";
        DEBUG_PRINT("RTTP Procedure started...");

        RTR_search(mcgrp);

        DEBUG_PRINT("RTTP Procedure done!");

        //Second Infeasible solve
        DEBUG_PRINT("IDP Procedure started...");

        infeasible_exploration(mcgrp);

        DEBUG_PRINT("IDP Procedure done!");

        My_Assert(total_vio_load == 0, "IDP procedure produce an infeasible solution");

        DEBUG_PRINT("RTTP Procedure started...");

        RTR_search(mcgrp);

        DEBUG_PRINT("RTTP Procedure done!");

    }
    else if (mode == IDPRTR) {
        cout << "IDP -> RTTP\n";
        DEBUG_PRINT("IDP Procedure started...");

        infeasible_exploration(mcgrp);

        DEBUG_PRINT("IDP Procedure done!");

        My_Assert(total_vio_load == 0, "IDP procedure produce an infeasible solution");

        DEBUG_PRINT("RTTP Procedure started...");

        RTR_search(mcgrp);

        DEBUG_PRINT("RTTP Procedure done!");
    }
    else {
        My_Assert(false, "Unknown search policy");
    }
}

void HighSpeedNeighBorSearch::unpack_seq(const std::vector<int> &dummy_seq, const MCGRP &mcgrp)
{

    My_Assert(dummy_seq.size() > 2, "You can't unpack an empty sequence!");

    //parse solution
    TASK_NODE *dummy = nullptr;
    dummy = solution.dummypool.get_new_dummy();
    My_Assert(solution.very_start == nullptr, "Incorrect state of solution");
    solution.very_start = dummy;

    //link first dummy
    My_Assert(dummy_seq[0] == DUMMY && dummy_seq[1] != DUMMY, "sequence should start from dummy Task!");
    solution.very_start->next = &solution.tasks[dummy_seq[1]];
    solution.tasks[dummy_seq[1]].pre = solution.very_start;

    int cursor = 1;
    while (true) {
        if (dummy_seq[cursor] == DUMMY) {
            if (cursor == dummy_seq.size() - 1) {
                //...0
                My_Assert(dummy->next == nullptr, "This dummy Task has been used!");
                solution.very_end = dummy;
                break;
            }

            My_Assert(dummy_seq[cursor + 1] != DUMMY, "too much dummy seq!");

            //...0-a...
            My_Assert(dummy != nullptr, "null dummy Task!");
            My_Assert(dummy->next == nullptr, "this dummy Task has been used!");
            dummy->next = &solution.tasks[dummy_seq[cursor + 1]];
            solution.tasks[dummy_seq[cursor + 1]].pre = dummy;

            cursor++;
        }
        else {
            if (dummy_seq[cursor + 1] == DUMMY) {
                //...a-0...
                // a new dummy Task is going to be created
                dummy = solution.dummypool.get_new_dummy();

                //connect dummy Task;
                solution.tasks[dummy_seq[cursor]].next = dummy;
                dummy->pre = &solution.tasks[dummy_seq[cursor]];
                cursor++;
            }
            else {
                //...a-b...
                solution.tasks[dummy_seq[cursor]].next = &solution.tasks[dummy_seq[cursor + 1]];
                solution.tasks[dummy_seq[cursor + 1]].pre = &solution.tasks[dummy_seq[cursor]];
                cursor++;
            }
        }
    }


    vector<vector<int>> seg;
    vector<int> buffer;
    for (cursor = 1; cursor < dummy_seq.size(); cursor++) {
        if (dummy_seq[cursor] == DUMMY) {
            seg.push_back(buffer);
            buffer.clear();
        }
        else {
            buffer.push_back(dummy_seq[cursor]);
        }
    }

    total_vio_load = 0;
    total_vio_time = 0;
    cur_solution_cost = 0;
    for (int i = 0; i < seg.size(); i++) {
        RouteInfo *new_route = routes[routes.allocate_route()];
        new_route->start = seg[i].front();
        new_route->end = seg[i].back();
        new_route->length = 0;
        new_route->load = 0;
        new_route->num_customers = seg[i].size();
        new_route->num_edges = 0;
        vector<int> arrive_time = move(mcgrp.cal_arrive_time(seg[i]));
        for (int j = 0; j < seg[i].size(); j++){
            new_route->time_table.push_back({seg[i][j], arrive_time[j]});
            total_vio_time +=
                max(0,
                    new_route->time_table.back().arrive_time - mcgrp.inst_tasks[seg[i][j]].time_window.second);
            if(mcgrp.is_edge(seg[i][j])) new_route->num_edges++;
        }

        //0-a...
        new_route->length +=
            mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[new_route->start].head_node];
        new_route->length += mcgrp.inst_tasks[DUMMY].serv_cost;

        for (int j = 0; j < seg[i].size() - 1; j++) {
            //...b-c...
            new_route->length += mcgrp.inst_tasks[seg[i][j]].serv_cost;
            new_route->length +=
                mcgrp.min_cost[mcgrp.inst_tasks[seg[i][j]].tail_node][mcgrp.inst_tasks[seg[i][j + 1]].head_node];
            new_route->load += mcgrp.inst_tasks[seg[i][j]].demand;
            My_Assert(solution[seg[i][j]]->route_id == -1, "The Task has been parsed!");
            solution[seg[i][j]]->route_id = new_route->ID;
        }

        //...d-0
        new_route->length += mcgrp.inst_tasks[new_route->end].serv_cost;
        new_route->length +=
            mcgrp.min_cost[mcgrp.inst_tasks[new_route->end].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
        new_route->length += mcgrp.inst_tasks[DUMMY].serv_cost;

        new_route->load += mcgrp.inst_tasks[new_route->end].demand;
        My_Assert(solution[new_route->end]->route_id == -1, "The Task has been parsed!");
        solution[new_route->end]->route_id = new_route->ID;

        total_vio_load += max((new_route->load - mcgrp.capacity), 0);
        cur_solution_cost += new_route->length;
    }

    My_Assert(valid_sol(mcgrp), "Wrong state");
}

vector<int> HighSpeedNeighBorSearch::get_solution(string mode)
{
    My_Assert(solution.very_start != nullptr, "Empty solution");
    vector<int> buffer;

    if (mode == "dummy") {
        TASK_NODE *tmp = solution.very_start;
        do {
            buffer.push_back(max(tmp->ID, 0));
            tmp = tmp->next;
        }
        while (tmp != solution.very_end);

        buffer.push_back(max(tmp->ID, 0));
    }
    else if (mode == "negative") {
        int sign = 1;
        TASK_NODE *tmp = solution.very_start;
        do {
            if (tmp->ID < 0) {
                sign = -1;
                tmp = tmp->next;
                continue;
            }
            buffer.push_back(sign * tmp->ID);
            tmp = tmp->next;
            sign = 1;
        }
        while (tmp != solution.very_end);

        if (tmp->ID > 0)
            buffer.push_back(tmp->ID);
    }
    else {
        cerr << "Unknown output style!\n";
        abort();
    }

    return buffer;
}

void HighSpeedNeighBorSearch::clear()
{
    total_vio_load = 0;
    best_solution_cost = numeric_limits<decltype(best_solution_cost)>::max();
    cur_solution_cost = 0;
    solution.clear();
    routes.clear();
}

void HighSpeedNeighBorSearch::create_individual(const MCGRP &mcgrp, Individual &p)
{
    p.sequence = get_solution();

    My_Assert(!routes.activated_route_id.empty(), "Routes info is empty!");
    p.route_seg_load.clear();
    p.total_vio_load = 0;
    for (auto id : routes.activated_route_id) {
        p.route_seg_load.push_back(routes[id]->load);
        p.total_vio_load += max((routes[id]->load - mcgrp.capacity), 0);
    }

    My_Assert(p.total_vio_load == total_vio_load, "incorrect total vio load!");
    p.total_cost = cur_solution_cost;
}

void HighSpeedNeighBorSearch::RTR_search(const MCGRP &mcgrp)
{
    double orig_val_for_uphill;
    double orig_val_for_downhill;
    local_minimum_likelihood = 1;
    int cnt = 1;
    do {
        struct timeb start_time;
        ftime(&start_time);

        threshold_exploration_version_0(mcgrp);
        orig_val_for_uphill = cur_solution_cost;

        do {
            orig_val_for_downhill = cur_solution_cost;
            descent_exploration_version_0(mcgrp);
        }
        while (cur_solution_cost < orig_val_for_downhill);

        struct timeb end_time;
        ftime(&end_time);

        cout << "Finish a feasible search process, spent: "
             << fixed << get_time_difference(start_time, end_time) << 's' << endl;

        if (cur_solution_cost < orig_val_for_uphill) {
            DEBUG_PRINT("total_route_length " + to_string(orig_val_for_uphill));
            DEBUG_PRINT("reset local likelihood");
            local_minimum_likelihood = 1;
            cnt++;
            if (cnt == max_RTR_search_cycle) {
                DEBUG_PRINT("maybe a local optimal");
                local_minimum_likelihood = local_threshold;
            }
        }
        else {
            local_minimum_likelihood++;
        }
    }
    while (local_minimum_likelihood < local_threshold);
}

void HighSpeedNeighBorSearch::descent_search(const MCGRP &mcgrp)
{
    descent_exploration_version_0(mcgrp);
}

void HighSpeedNeighBorSearch::create_search_neighborhood(const MCGRP &mcgrp, const vector<int>& chosen_seq,string mode, int offset)
{
    int original_neigh_size = neigh_size;
    search_space.clear();
    if(mode.empty() || chosen_seq.empty()) return;

    int start_node = mcgrp.inst_tasks[chosen_seq.front()].head_node;
    int end_node = mcgrp.inst_tasks[chosen_seq.back()].tail_node;

    const NeighborInfo& neighbor_info = mcgrp.neighbor[start_node][end_node];

    int loc = offset;
    if(mode == "basic"){

        neigh_size = min(neigh_size, (int)neighbor_info.basic_neighbor.size());
        while(search_space.size() < neigh_size){
            loc %= neighbor_info.basic_neighbor.size();
            search_space.push_back(neighbor_info.basic_neighbor[loc].task_id);
            loc++;
        }

    }else if(mode == "predecessor"){
        const auto& predecessor_seq = neighbor_info.predecessor_neighbor.at(chosen_seq.front());

        neigh_size = min(neigh_size, (int)predecessor_seq.size());
        while(search_space.size() < neigh_size){
            loc %= predecessor_seq.size();
            search_space.push_back(predecessor_seq[loc].task_id);
            loc++;
        }

    }else if(mode == "successor"){
        const auto& successor_seq = neighbor_info.successor_neighbor.at(chosen_seq.back());

        neigh_size = min(neigh_size,(int)successor_seq.size());
        while(search_space.size() < neigh_size){
            loc %= successor_seq.size();
            search_space.push_back(successor_seq[loc].task_id);
            loc++;
        }

    }else {
        My_Assert(false, "Unknown arguments");
    }


    mcgrp._rng.RandPerm(search_space);
    neigh_size = original_neigh_size;
}

bool HighSpeedNeighBorSearch::valid_sol(const MCGRP &mcgrp)
{
    if (solution.very_start->ID >= 0 || solution.very_end->ID >= 0) {
        return false;
    }

    auto cur_task = solution.very_start;

    double cost = 0;
    while (cur_task != solution.very_end) {
        cost += mcgrp.inst_tasks[max(cur_task->ID, 0)].serv_cost;
        cost += mcgrp.min_cost[mcgrp.inst_tasks[max(cur_task->ID, 0)].tail_node]
        [mcgrp.inst_tasks[max(cur_task->next->ID, 0)].head_node];
        cur_task = cur_task->next;
    }
    cost += mcgrp.inst_tasks[max(cur_task->ID, 0)].serv_cost;

    double vio_load = 0;
    double routes_cost_sum = 0;
    int vio_time = 0;
    for (auto id : routes.activated_route_id) {
        if(routes[id]->num_customers != routes[id]->time_table.size()
            || routes[id]->num_edges < 0
            || routes[id]->num_edges > routes[id]->num_customers)
            return false;
        if (routes[id]->load > mcgrp.capacity) {
            vio_load += routes[id]->load - mcgrp.capacity;
        }
        routes_cost_sum += routes[id]->length;
        for (auto plan : routes[id]->time_table) {
            vio_time += max(0, plan.arrive_time - mcgrp.inst_tasks[plan.task].time_window.second);
        }
    }

    return routes_cost_sum == cost && cost == cur_solution_cost
        && vio_load == total_vio_load && vio_time == total_vio_time;
}

void HighSpeedNeighBorSearch::threshold_exploration_version_0(const MCGRP &mcgrp)
{
    DEBUG_PRINT("Uphill and downhill...");

    // Based on the best solution find so far, experiment shows this is better than the policy,especially on big instance
    // which start from current solution.
    policy.benchmark = this->best_solution_cost;


    policy.tolerance = sel_ratio(prob, ratios, mcgrp._rng);
    const auto original_policy = policy.get();
    policy.set(BEST_ACCEPT | TOLERANCE | DELTA_ONLY);

    //you need to decide the dynamic neighbor size when you search based on different policy
//    neigh_size = mcgrp.neigh_size;
//    neigh_size = min(10,mcgrp.neigh_size);


    vector<NeighborOperator> neighbor_operator{
        NeighborOperator::SINGLE_INSERT,
        NeighborOperator::DOUBLE_INSERT,
        NeighborOperator::SWAP,
        NeighborOperator::INVERT,
        NeighborOperator::SLICE,
        NeighborOperator::EXTRACTION,
        NeighborOperator::TWO_OPT
    };

    int L = mcgrp._rng.Randint(28, 33);

    for (int k = 1; k < L; k++) {
        mcgrp._rng.RandPerm(neighbor_operator);

        int chosen_task = -1;

        for (auto cur_operator:neighbor_operator) {
            mcgrp._rng.RandPerm(task_set);

            for (int i = 0; i < mcgrp.actual_task_num; i++) {
                chosen_task = task_set[i];
                if (solution[chosen_task]->next == nullptr) {
                    if (mcgrp.is_edge(chosen_task)) {
                        chosen_task = mcgrp.inst_tasks[chosen_task].inverse;
                        My_Assert(solution[chosen_task]->next != nullptr, "An edge Task has been missed");
                    }
                    else {
                        My_Assert(false, "A non edge Task has been missed!");
                    }
                }

                switch (cur_operator) {
                    case SINGLE_INSERT:
                        single_pre_insert->search(*this, mcgrp, chosen_task);
                        single_post_insert->search(*this, mcgrp,chosen_task);
                        break;
                    case DOUBLE_INSERT:
                        double_pre_insert->search(*this, mcgrp, chosen_task);
                        double_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case SWAP:swap->search(*this, mcgrp, chosen_task);
                        break;
                    case INVERT:
                        if (mcgrp.req_edge_num != 0) {
                            invert->search(*this, mcgrp, chosen_task);
                        }
                        break;
                    case TWO_OPT:two_opt->search(*this, mcgrp, chosen_task);
                        break;
                    case SLICE:slice->search(*this, mcgrp, chosen_task);
                        break;
                    case EXTRACTION:extraction->search(*this, mcgrp, chosen_task);
                        break;
                    default:My_Assert(false, "unknown operator!");
                }

            }
        }
    }

    policy.benchmark = 0;
    policy.set(original_policy);
    policy.tolerance = 0;
    neigh_size = 0;
}

void HighSpeedNeighBorSearch::descent_exploration_version_0(const MCGRP &mcgrp)
{
    DEBUG_PRINT("Downhill...");

    const auto original_policy = policy.get();
    policy.set(BEST_ACCEPT | DOWNHILL | DELTA_ONLY);

//    neigh_size = mcgrp.neigh_size;

    vector<NeighborOperator> neighbor_operator{
        NeighborOperator::SINGLE_INSERT,
        NeighborOperator::DOUBLE_INSERT,
        NeighborOperator::INVERT,
        NeighborOperator::SWAP,
        NeighborOperator::TWO_OPT,
    };

    int chosen_task = -1;
    mcgrp._rng.RandPerm(neighbor_operator);
    for (auto cur_operator : neighbor_operator) {
        double start_val = 0;
        do {
            start_val = cur_solution_cost;
            mcgrp._rng.RandPerm(task_set);

            for (int i = 0; i < mcgrp.actual_task_num; i++) {
                chosen_task = task_set[i];
                if (solution[chosen_task]->next == nullptr) {
                    if (mcgrp.is_edge(chosen_task)) {
                        chosen_task = mcgrp.inst_tasks[chosen_task].inverse;
                        My_Assert(solution[chosen_task]->next != nullptr, "An edge Task has been missed");
                    }
                    else {
                        My_Assert(false, "A non edge Task has been missed!");
                    }
                }

                switch (cur_operator) {
                    case SINGLE_INSERT:
                        single_pre_insert->search(*this, mcgrp, chosen_task);
                        single_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case DOUBLE_INSERT:
                        double_pre_insert->search(*this, mcgrp, chosen_task);
                        double_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case SWAP:swap->search(*this, mcgrp, chosen_task);
                        break;
                    case INVERT:
                        if (mcgrp.req_edge_num != 0) {
                            invert->search(*this, mcgrp, chosen_task);
                        }
                        break;
                    case TWO_OPT:two_opt->search(*this, mcgrp, chosen_task);
                        break;
                    default:My_Assert(false, "unknown operator!");
                }
            }
        }
        while (cur_solution_cost < start_val);
    }

    policy.set(original_policy);
    neigh_size = 0;
}

void HighSpeedNeighBorSearch::infeasible_exploration(const MCGRP &mcgrp)
{
    // Based on the best solution find so far, experiment shows this is better than the policy
    // which start from current solution.
    this->clear();
    unpack_seq(get_delimiter_coding(mcgrp.best_sol_buff), mcgrp);

    My_Assert(total_vio_load == 0, "Cannot start from an infeasible point!");
    policy.nearest_feasible_cost = cur_solution_cost;
    policy.beta = cur_solution_cost / double(mcgrp.capacity * 15);

    small_step_infeasible_descent_search(mcgrp);

    DEBUG_PRINT("Trigger Infeasible Tabu Search");
    small_step_infeasible_tabu_search(mcgrp);

    small_step_infeasible_descent_search(mcgrp);


    //Here used to break the local minimum with merge-split operator
    large_step_infeasible_search(mcgrp);

    My_Assert(total_vio_load >= 0, "Wrong total violated load!");
    if (total_vio_load > 0 || total_vio_time > 0) {
        repair_solution(mcgrp);
    }

    policy.nearest_feasible_cost = 0;
    policy.beta = 0;
    this->best_solution_cost = cur_solution_cost;
}

void HighSpeedNeighBorSearch::trace(const MCGRP &mcgrp)
{
    for (int i : routes.activated_route_id) {
        My_Assert(routes[i]->time_table.size() == routes[i]->num_customers, "Wrong size");
    }

    if (total_vio_load == 0 && total_vio_time == 0) {
        if (cur_solution_cost < this->best_solution_cost) {
            this->best_solution_cost = cur_solution_cost;
        }

        vector<int> negative_coding_sol;
        if (cur_solution_cost < mcgrp.best_total_route_length) {
            negative_coding_sol = get_solution("negative");
        }

        mcgrp.check_best_solution(cur_solution_cost, negative_coding_sol);
    }
}

void HighSpeedNeighBorSearch::small_step_infeasible_descent_search(const MCGRP &mcgrp)
{
    auto original_policy = policy.get();
    policy.set(BEST_ACCEPT | DOWNHILL | FITNESS_ONLY);
    neigh_size = 10;

    vector<NeighborOperator> neighbor_operator{
        NeighborOperator::SINGLE_INSERT,
        NeighborOperator::DOUBLE_INSERT,
        NeighborOperator::SWAP,
        NeighborOperator::TWO_OPT
    };

    equal_step = 0;
    search_step = 0;

    int prior_search_step;
    int prior_equal_step;

    int search_step_delta;
    int equal_step_delta;

    double equal_ratio;

    static int too_far = 0;
    static int too_near = 0;
    static int count = 0;

    do {
        prior_equal_step = this->equal_step;
        prior_search_step = this->search_step;
        int chosen_task = -1;
        mcgrp._rng.RandPerm(neighbor_operator);

        if (total_vio_load > 0) {
            //start location should not be considered
            count++;

            double distance_to_feasible_zone = double(total_vio_load) / policy.nearest_feasible_cost;

            if (distance_to_feasible_zone > infeasible_distance_threshold)
                too_far++;
            else
                too_near++;

            if (count % 5 == 0) {

                if (too_near == 5) {
                    DEBUG_PRINT("Infeasible search too near");
                    policy.beta /= 2;
                }
                else if (too_far == 5) {
                    DEBUG_PRINT("Infeasible search too far");
                    policy.beta *= 2;
                }

                too_far = 0;
                too_near = 0;
                count = 0;
            }
        }


        for (auto cur_operator:neighbor_operator) {
            mcgrp._rng.RandPerm(task_set);

            for (int i = 0; i < mcgrp.actual_task_num; i++) {
                chosen_task = task_set[i];
                if (solution[chosen_task]->next == nullptr) {
                    if (mcgrp.is_edge(chosen_task)) {
                        chosen_task = mcgrp.inst_tasks[chosen_task].inverse;
                        My_Assert(
                            solution[chosen_task]->next != nullptr, "An edge Task has been missed");
                    }
                    else {
                        My_Assert(false, "A non edge Task has been missed!");
                    }
                }

                switch (cur_operator) {
                    case SINGLE_INSERT:
                        single_pre_insert->search(*this, mcgrp, chosen_task);
                        single_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case DOUBLE_INSERT:
                        double_pre_insert->search(*this, mcgrp, chosen_task);
                        double_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case SWAP:swap->search(*this, mcgrp, chosen_task);
                        break;
                    case TWO_OPT:two_opt->search(*this, mcgrp, chosen_task);
                        break;
                    default:My_Assert(false, "unknown operator!");
                }
            }
        }

        search_step_delta = search_step - prior_search_step;
        equal_step_delta = equal_step - prior_equal_step;

        if (search_step_delta == 0 && equal_step_delta == 0) {
            equal_ratio = 0;
        }
        else {
            equal_ratio = double(equal_step_delta) / double(search_step_delta);
        }

        My_Assert(equal_ratio >= 0 && equal_ratio <= 1, "Wrong ratio!");
    }
    while (
        search_step_delta > significant_search_delta
            && equal_ratio < local_ratio_threshold
        );


    policy.set(original_policy);
    neigh_size = 0;
}

void HighSpeedNeighBorSearch::small_step_infeasible_tabu_search(const MCGRP &mcgrp)
{
    policy.benchmark = this->best_solution_cost;
    policy.tolerance = sel_ratio(prob, ratios, mcgrp._rng);
    neigh_size = 10;

    auto original_policy = policy.get();
    policy.set(BEST_ACCEPT | TOLERANCE | FITNESS_ONLY);

    vector<NeighborOperator> neighbor_operator{
        NeighborOperator::SINGLE_INSERT,
        NeighborOperator::DOUBLE_INSERT,
        NeighborOperator::SWAP,
        NeighborOperator::TWO_OPT
    };

    int L = mcgrp._rng.Randint(28, 33);

    for (int k = 0; k < L; k++) {
        int chosen_task = -1;
        mcgrp._rng.RandPerm(neighbor_operator);

        for (auto cur_operator:neighbor_operator) {
            mcgrp._rng.RandPerm(task_set);

            for (int i = 0; i < mcgrp.actual_task_num; i++) {
                chosen_task = task_set[i];
                if (solution[chosen_task]->next == nullptr) {
                    if (mcgrp.is_edge(chosen_task)) {
                        chosen_task = mcgrp.inst_tasks[chosen_task].inverse;
                        My_Assert(
                            solution[chosen_task]->next != nullptr, "An edge Task has been missed");
                    }
                    else {
                        My_Assert(false, "A non edge Task has been missed!");
                    }
                }

                switch (cur_operator) {
                    case SINGLE_INSERT:
                        single_pre_insert->search(*this, mcgrp, chosen_task);
                        single_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case DOUBLE_INSERT:
                        double_pre_insert->search(*this, mcgrp, chosen_task);
                        double_post_insert->search(*this, mcgrp, chosen_task);
                        break;
                    case SWAP:swap->search(*this, mcgrp, chosen_task);
                        break;
                    case TWO_OPT:two_opt->search(*this, mcgrp, chosen_task);
                        break;
                    default:My_Assert(false, "unknown operator!");
                }
            }
        }

    }

    policy.benchmark = 0;
    policy.set(original_policy);
    policy.tolerance = 0;
    neigh_size = 0;
}

void HighSpeedNeighBorSearch::large_step_infeasible_search(const MCGRP &mcgrp)
{
    DEBUG_PRINT("Trigger large step break movement");
    // decide how many routes need to be merged, self-adaptive
//    int merge_size = routes.size() / 3;
    int merge_size = 10;

    // enlarge the capacity for infeasible search to decide whether use infeasible consideration
//    double scale = mcgrp.capacity * policy.beta;
//    int pseudo_capacity = int(mcgrp.capacity *(1+ (scale / (policy.beta*sqrt(1+scale*scale)))));
    int pseudo_capacity = mcgrp.capacity;
    merge_split(*this, mcgrp, merge_size, pseudo_capacity);

    trace(mcgrp);

//    mcgrp.check_best_infeasible_solution(cur_solution_cost,beta,total_vio_load,negative_coding_sol);
}

void HighSpeedNeighBorSearch::update(const MCGRP &mcgrp,
                                     const vector<int> &best_buffer,
                                     const vector<int> &best_routes)
{
    My_Assert(valid_sol(mcgrp), "Wrong state!");

    //remove old routes
    for (auto route_id : best_routes) {
        My_Assert(routes.activated_route_id.find(route_id) != routes.activated_route_id.end(), "Invalid route");
        cur_solution_cost -= routes[route_id]->length;
        if (routes[route_id]->load > mcgrp.capacity) {
            total_vio_load -= (routes[route_id]->load - mcgrp.capacity);
        }
        total_vio_time -= mcgrp.get_vio_time(routes[route_id]->time_table);

        //remove solution
        vector<int> seq;
        int cur_task = routes[route_id]->start;
        while (cur_task != routes[route_id]->end) {
            seq.push_back(cur_task);
            cur_task = solution[cur_task]->next->ID;
        }
        seq.push_back(cur_task);

        int dummy_marker = 0;
        HighSpeedNeighBorSearch::TASK_NODE *reserve_dummy;
        bool very_end_case = false;
        if (solution[seq.back()]->next == solution.very_end) {
            dummy_marker = solution[seq.front()]->pre->ID;
            reserve_dummy = solution[seq.back()]->next;
            very_end_case = true;
        }
        else {
            dummy_marker = solution[seq.back()]->next->ID;
            reserve_dummy = solution[seq.front()]->pre;
        }
        My_Assert(dummy_marker < 0, "Wrong tasks");

        if (very_end_case) {
            solution[dummy_marker]->pre->next = reserve_dummy;
            reserve_dummy->pre = solution[dummy_marker]->pre;
        }
        else {
            reserve_dummy->next = solution[dummy_marker]->next;
            solution[dummy_marker]->next->pre = reserve_dummy;
        }

        solution.dummypool.free_dummy(dummy_marker);

        for (auto task: seq) {
            solution[task]->clear();
        }

        //free old route
        routes.free_route(route_id);
    }

    My_Assert(valid_sol(mcgrp), "Wrong state!");

    My_Assert(best_buffer.front() == DUMMY && best_buffer.back() == DUMMY, "Wrong state");


    vector<vector<int >> seqs;
    vector<int> buffer;
    for (int i = 1; i < best_buffer.size(); i++) {
        if (best_buffer[i] == DUMMY) {
            seqs.push_back(buffer);
            buffer.clear();
        }
        else {
            buffer.push_back(best_buffer[i]);
        }
    }

    // insert new routes
    for (int row = 0; row < seqs.size(); row++) {
        int new_route = routes.allocate_route();
        double
            length = mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[seqs[row].front()].head_node];
        double load = 0;

        for (int col = 0; col < seqs[row].size() - 1; col++) {
            length += mcgrp.inst_tasks[seqs[row][col]].serv_cost
                + mcgrp.min_cost[mcgrp.inst_tasks[seqs[row][col]].tail_node][mcgrp.inst_tasks[seqs[row][col + 1]]
                    .head_node];

            load += mcgrp.inst_tasks[seqs[row][col]].demand;

            solution[seqs[row][col]]->route_id = new_route;
        }

        length += mcgrp.inst_tasks[seqs[row].back()].serv_cost
            + mcgrp.min_cost[mcgrp.inst_tasks[seqs[row].back()].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

        load += mcgrp.inst_tasks[seqs[row].back()].demand;
        solution[seqs[row].back()]->route_id = new_route;


        routes[new_route]->length = length;
        routes[new_route]->load = load;
        routes[new_route]->start = seqs[row].front();
        routes[new_route]->end = seqs[row].back();
        routes[new_route]->num_customers = seqs[row].size();
        vector<int> time_tbl = mcgrp.cal_arrive_time(seqs[row]);
        routes[new_route]->time_table.clear();
        for(int i = 0;i< time_tbl.size();i++){
            routes[new_route]->time_table.push_back({seqs[row][i],time_tbl[i]});
        }

        total_vio_time += mcgrp.get_vio_time(routes[new_route]->time_table);
        if (load > mcgrp.capacity) {
            total_vio_load += (load - mcgrp.capacity);
        }

        cur_solution_cost += length;

        //handle solution
        if (seqs[row].size() == 1) {
            auto task = seqs[row][0];
            auto new_dummy = solution.dummypool.get_new_dummy();
            new_dummy->pre = solution.very_end->pre;
            new_dummy->next = solution[task];
            new_dummy->pre->next = new_dummy;
            solution[task]->pre = new_dummy;
            solution[task]->next = solution.very_end;
            solution.very_end->pre = solution[task];
        }
        else {
            My_Assert(seqs[row].size() > 1, "Wrong state");

            auto new_dummy = solution.dummypool.get_new_dummy();
            new_dummy->pre = solution.very_end->pre;
            new_dummy->next = solution[seqs[row].front()];

            new_dummy->pre->next = new_dummy;

            solution[seqs[row].front()]->pre = new_dummy;
            solution[seqs[row].front()]->next = solution[seqs[row][1]];

            for (auto i = 1; i < seqs[row].size() - 1; i++) {
                auto pre_task = seqs[row][i - 1];
                auto task = seqs[row][i];
                auto next_task = seqs[row][i + 1];
                My_Assert(solution[task]->pre == nullptr && solution[task]->next == nullptr, "Wrong state!");
                solution[task]->next = solution[next_task];
                solution[task]->pre = solution[pre_task];
            }

            solution[seqs[row].back()]->pre = solution[seqs[row][seqs[row].size() - 2]];
            solution[seqs[row].back()]->next = solution.very_end;

            solution.very_end->pre = solution[seqs[row].back()];

        }

    }

    My_Assert(missed(mcgrp) && check_duplicated(mcgrp), "Wrong result!");
    My_Assert(valid_sol(mcgrp), "Wrong state!");
}

void HighSpeedNeighBorSearch::repair_solution(const MCGRP &mcgrp)
{
    // Repair infeasible solution
    // Only violated load and time window here
    // No missed Task, duplicated Task problem here

    My_Assert(missed(mcgrp), "Some Task missed!");
    My_Assert(check_duplicated(mcgrp), "Duplicated Task!");

    My_Assert(total_vio_load > 0 || total_vio_time > 0, "This is not a infeasible Task!");

/*
    if (total_vio_load > 0)
        _repair_load(mcgrp);

    My_Assert(total_vio_load == 0, "load reparation mistake!");

    if (total_vio_time > 0)
        _repair_time_window(mcgrp);

    My_Assert(total_vio_time == 0, "time window reparation mistake!");
*/

    _tour_splitting_repair(mcgrp);
    My_Assert(total_vio_load > 0 || total_vio_time > 0, "This is not a infeasible Task!");

}

void HighSpeedNeighBorSearch::_repair_load(const MCGRP &mcgrp)
{
    DEBUG_PRINT("Repair Load infeasible solution...");

    struct Route
    {
        int route_id;
        vector<int> task_seq;
        int load;
    };

    vector<Route> satisfied_routes;
    vector<Route> violated_routes;

    //generate tasks distribution in routes
    for (auto route_id : routes.activated_route_id) {
        Route tmp;
        tmp.route_id = route_id;
        tmp.load = routes[route_id]->load;
        int cur_task = routes[route_id]->start;
        while (cur_task != routes[route_id]->end) {
            tmp.task_seq.push_back(cur_task);
            cur_task = solution[cur_task]->next->ID;
        }
        tmp.task_seq.push_back(cur_task);

        if (tmp.load > mcgrp.capacity) {
            violated_routes.push_back(tmp);
        }
        else {
            satisfied_routes.push_back(tmp);
        }
    }

    My_Assert(!violated_routes.empty(), "Wrong state");
    //the Task set which need to be re-inserted
    vector<int> candidate_tasks;
    for (auto current_route : violated_routes) {
        //remove the violated routes info
        cur_solution_cost -= routes[current_route.route_id]->length;
        total_vio_load -= (routes[current_route.route_id]->load - mcgrp.capacity);
        total_vio_time -= mcgrp.get_vio_time(routes[current_route.route_id]->time_table);

        delete_route(current_route.route_id, current_route.task_seq);
        candidate_tasks.insert(candidate_tasks.end(),
                               current_route.task_seq.begin(),
                               current_route.task_seq.end());
    }

    My_Assert(valid_sol(mcgrp), "Wrong validation");

    //insert Task to repair
    HighSpeedNeighBorSearch::TASK_NODE *left_task;
    HighSpeedNeighBorSearch::TASK_NODE *right_task;
    int chosen_route = -1;
    int best_delta = numeric_limits<decltype(best_delta)>::max();
    int delta = numeric_limits<decltype(delta)>::max();
    for (auto task : candidate_tasks) {
        //reset chosen information before insert
        //for faster execution time, I omit inverse Task of edge Task!:)
        chosen_route = -1;
        best_delta = numeric_limits<decltype(best_delta)>::max();
        left_task = nullptr;
        right_task = nullptr;

        auto demand = mcgrp.inst_tasks[task].demand;

        // choose the best insert route
        for (int row = 0; row < satisfied_routes.size(); row++) {
            //load constraint check
            if (demand + satisfied_routes[row].load > mcgrp.capacity) {
                continue;
            }
            else {
                int a = 0;
                int b = 0;
                double ab = 0;
                double a_task = 0;
                double task_b = 0;

                //find a best position to insert
                a = solution[satisfied_routes[row].task_seq.front()]->pre->ID;
                My_Assert(a < 0, "Wrong Task");
                b = solution[satisfied_routes[row].task_seq.front()]->ID;
                My_Assert(b > 0, "Wrong Task");

                ab = mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[b].head_node];
                a_task = mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[task].head_node];
                task_b = mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[b].head_node];
                delta = -ab + a_task + mcgrp.inst_tasks[task].serv_cost + task_b;

                if (delta < best_delta) {
                    best_delta = delta;
                    chosen_route = row;
                    left_task = solution[a];
                    right_task = solution[b];
                }

                for (int col = 1; col < satisfied_routes[row].task_seq.size() - 1; col++) {
                    a = solution[satisfied_routes[row].task_seq[col]]->ID;
                    b = solution[satisfied_routes[row].task_seq[col + 1]]->ID;
                    My_Assert(a > 0, "Wrong Task");
                    My_Assert(b > 0, "Wrong Task");

                    ab = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[b].head_node];
                    a_task = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[task].head_node];
                    task_b = mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[b].head_node];
                    delta = -ab + a_task + mcgrp.inst_tasks[task].serv_cost + task_b;

                    if (delta < best_delta) {
                        best_delta = delta;
                        chosen_route = row;
                        left_task = solution[a];
                        right_task = solution[b];
                    }
                }

                a = solution[satisfied_routes[row].task_seq.back()]->ID;
                My_Assert(a > 0, "Wrong Task");
                b = solution[satisfied_routes[row].task_seq.back()]->next->ID;
                My_Assert(b < 0, "Wrong Task");

                ab = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
                a_task = mcgrp.min_cost[mcgrp.inst_tasks[a].tail_node][mcgrp.inst_tasks[task].head_node];
                task_b = mcgrp.min_cost[mcgrp.inst_tasks[task].tail_node][mcgrp.inst_tasks[DUMMY].head_node];
                delta = -ab + a_task + mcgrp.inst_tasks[task].serv_cost + task_b;

                if (delta < best_delta) {
                    best_delta = delta;
                    chosen_route = row;
                    left_task = solution[a];
                    right_task = solution[b];
                }
            }
        }

        //means a new route need to be created
        if (chosen_route == -1) {
            My_Assert(left_task == nullptr && right_task == nullptr, "Wrong tasks");

            satisfied_routes.push_back(Route());
            satisfied_routes.back().task_seq.push_back(task);
            satisfied_routes.back().load = mcgrp.inst_tasks[task].demand;
            satisfied_routes.back().route_id = new_route(mcgrp, {task});

        }
        else {
            int a = left_task->ID;
            int b = right_task->ID;
            if (a < 0) {
                satisfied_routes[chosen_route].task_seq.insert(satisfied_routes[chosen_route].task_seq.begin(), task);
            }
            else if (b < 0) {
                satisfied_routes[chosen_route].task_seq.insert(satisfied_routes[chosen_route].task_seq.end(), task);
            }
            else {
                auto ite = find(satisfied_routes[chosen_route].task_seq.begin(),
                                satisfied_routes[chosen_route].task_seq.end(),
                                b);
                My_Assert(ite != satisfied_routes[chosen_route].task_seq.end(), "Cannot find the Task!");
                satisfied_routes[chosen_route].task_seq.insert(ite, task);
            }

            satisfied_routes[chosen_route].load += mcgrp.inst_tasks[task].demand;


            const auto route_id = satisfied_routes[chosen_route].route_id;
            routes[route_id]->length += best_delta;
            routes[route_id]->num_customers += 1;
            routes[route_id]->load += demand;

            total_vio_time -= mcgrp.get_vio_time(routes[route_id]->time_table);
            vector<int> time_tbl = mcgrp.cal_arrive_time(satisfied_routes[chosen_route].task_seq);
            routes[route_id]->time_table.clear();
            for(int i = 0;i<time_tbl.size();i++)
                routes[route_id]->time_table.push_back({satisfied_routes[chosen_route].task_seq[i],time_tbl[i]});
            total_vio_time += mcgrp.get_vio_time(routes[route_id]->time_table);

            solution[task]->route_id = route_id;

            if (a < 0) {
                routes[route_id]->start = task;
            }
            else if (b < 0) {
                routes[route_id]->end = task;
            }

            cur_solution_cost += best_delta;
            //handle solution
            solution[a]->next = solution[task];
            solution[b]->pre = solution[task];
            solution[task]->pre = solution[a];
            solution[task]->next = solution[b];
        }

        My_Assert(valid_sol(mcgrp), "Repair method doesn't work properly!");
    }

    My_Assert(missed(mcgrp), "Some Task missed!");
    My_Assert(check_duplicated(mcgrp), "Duplicated Task!");
    My_Assert(total_vio_load == 0, "This is not a infeasible Task!");

    My_Assert(valid_sol(mcgrp), "Repair method doesn't work properly!");
}

void HighSpeedNeighBorSearch::_repair_time_window(const MCGRP &mcgrp)
{
    DEBUG_PRINT("Repair Time window infeasible solution...");

    struct Route
    {
        int route_id;
        vector<RouteInfo::TimeTable> time_tbl;
    };

    vector<Route> violated_routes;

    //generate tasks distribution in routes
    for (auto route_id : routes.activated_route_id) {
        Route tmp;
        tmp.route_id = route_id;
        tmp.time_tbl = routes[route_id]->time_table;
        if (!mcgrp.isTimeTableFeasible(tmp.time_tbl, false))
            violated_routes.push_back(tmp);
    }

    My_Assert(!violated_routes.empty(), "Wrong state");

    for (auto current_route : violated_routes) {
        cur_solution_cost -= routes[current_route.route_id]->length;
        total_vio_time -= mcgrp.get_vio_time(current_route.time_tbl);

        vector<int> candidate_tasks;
        for (const auto &node : current_route.time_tbl)
            candidate_tasks.push_back(node.task);
        delete_route(current_route.route_id, candidate_tasks);

        My_Assert(valid_sol(mcgrp), "Wrong validation");

        Individual individual =  nearest_scanning(mcgrp, candidate_tasks);

        for (const auto& cur_route : individual.time_tbl){
            vector<int> tasks_id;
            for(const auto& task : cur_route){
                tasks_id.push_back(task.task);
            }
            new_route(mcgrp, tasks_id);
        }
    }

    My_Assert(missed(mcgrp), "Some Task missed!");
    My_Assert(check_duplicated(mcgrp), "Duplicated Task!");
    My_Assert(total_vio_time == 0, "This is not a infeasible Task!");

    My_Assert(valid_sol(mcgrp), "Repair method doesn't work properly!");
}

bool HighSpeedNeighBorSearch::missed(const MCGRP &mcgrp)
{
    for (auto task_id : task_set) {
        if (solution.tasks[task_id].next == nullptr) {
            if (!mcgrp.is_edge(task_id)) {
                return false;
            }
            else {
                int inverse_task_id = mcgrp.inst_tasks[task_id].inverse;
                if (solution.tasks[inverse_task_id].next == nullptr) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool HighSpeedNeighBorSearch::check_duplicated(const MCGRP &mcgrp)
{
    //offset one
    vector<int> tasks(mcgrp.actual_task_num + 1, 0);

    vector<int> delimiter_coding_sol = get_solution();
    for (auto task : delimiter_coding_sol) {
        if (task == DUMMY) {
            continue;
        }

        if (tasks[task] != 0) {
            return false;
        }

        tasks[task]++;
        if (mcgrp.is_edge(task)) {
            tasks[mcgrp.inst_tasks[task].inverse]++;
        }
    }

    return true;
}

int HighSpeedNeighBorSearch::getSearch_step() const
{
    return search_step;
}

void HighSpeedNeighBorSearch::delete_route(int route_id, const vector<int> &route_seq)
{
    My_Assert(routes[route_id]->start == route_seq.front() && routes[route_id]->end == route_seq.back(),
              "Wrong route!");

    routes.free_route(route_id);

    //remove solution
    int dummy_marker = 0;
    HighSpeedNeighBorSearch::TASK_NODE *reserve_dummy;
    bool very_end_case = false;
    if (solution[route_seq.back()]->next == solution.very_end) {
        dummy_marker = solution[route_seq.front()]->pre->ID;
        reserve_dummy = solution[route_seq.back()]->next;
        very_end_case = true;
    }
    else {
        dummy_marker = solution[route_seq.back()]->next->ID;
        reserve_dummy = solution[route_seq.front()]->pre;
    }
    My_Assert(dummy_marker < 0, "Wrong tasks");

    if (very_end_case) {
        My_Assert(reserve_dummy == solution.very_end, "Wrong state");
        if (solution[dummy_marker] != solution.very_start) {
            solution[dummy_marker]->pre->next = reserve_dummy;
            reserve_dummy->pre = solution[dummy_marker]->pre;
            solution.dummypool.free_dummy(dummy_marker);
        }
        else {
            solution[dummy_marker]->next = reserve_dummy;
            reserve_dummy->pre = solution[dummy_marker];
        }

    }
    else {
        reserve_dummy->next = solution[dummy_marker]->next;
        solution[dummy_marker]->next->pre = reserve_dummy;
        solution.dummypool.free_dummy(dummy_marker);
    }


    for (auto task: route_seq)
        solution[task]->clear();
}

int HighSpeedNeighBorSearch::new_route(const MCGRP &mcgrp, const vector<int> &route_seq)
{
    if(route_seq.empty()) return -1;

    const auto new_route = routes.allocate_route();

    const auto
        dummy_task = mcgrp.min_cost[mcgrp.inst_tasks[DUMMY].tail_node][mcgrp.inst_tasks[route_seq.front()].head_node];
    const auto
        task_dummy = mcgrp.min_cost[mcgrp.inst_tasks[route_seq.back()].tail_node][mcgrp.inst_tasks[DUMMY].head_node];

    // update meta info
    routes[new_route]->num_customers = route_seq.size();
    routes[new_route]->start = route_seq.front();
    routes[new_route]->end = route_seq.back();

    // update length and load info
    routes[new_route]->length = dummy_task + task_dummy;
    routes[new_route]->load = 0;
    routes[new_route]->num_edges = 0;
    for (const int task : route_seq) {
        routes[new_route]->load += mcgrp.inst_tasks[task].demand;
        routes[new_route]->length += mcgrp.inst_tasks[task].serv_cost;
        solution[task]->route_id = new_route;
        if(mcgrp.is_edge(task)) routes[new_route]->num_edges++;
    }

    for (int i = 1; i < route_seq.size(); i++)
        routes[new_route]->length +=
            mcgrp.min_cost[mcgrp.inst_tasks[route_seq[i - 1]].tail_node][mcgrp.inst_tasks[route_seq[i]].head_node];

    // update time table
    routes[new_route]->time_table =
        RouteInfo::TimeTable::zip(route_seq, mcgrp.cal_arrive_time(route_seq));

    //handle solution
    auto new_dummy = solution.dummypool.get_new_dummy();
    auto ends = solution.connect_tasks(route_seq);

    new_dummy->pre = solution.very_end->pre;
    new_dummy->pre->next = new_dummy;
    new_dummy->next = ends.first;
    ends.first->pre = new_dummy;

    ends.second->next = solution.very_end;
    solution.very_end->pre = ends.second;

    return new_route;
}

void HighSpeedNeighBorSearch::_tour_splitting_repair(const MCGRP &mcgrp)
{
    DEBUG_PRINT("Tour splitting repair infeasible solution...");

    struct Route
    {
        int route_id;
        int load;
        vector<RouteInfo::TimeTable> time_tbl;
    };

    vector<Route> satisfied_routes;
    vector<Route> violated_routes;

    for (auto route_id : routes.activated_route_id) {
        Route tmp;
        tmp.route_id = route_id;
        tmp.load = routes[route_id]->load;
        tmp.time_tbl = routes[route_id]->time_table;

        if (tmp.load > mcgrp.capacity || !mcgrp.isTimeTableFeasible(tmp.time_tbl, false))
            violated_routes.push_back(tmp);
        else
            satisfied_routes.push_back(tmp);
    }

    My_Assert(!violated_routes.empty(), "Wrong state");

    vector<int> task_list;
    for (auto current_route : violated_routes) {
        vector<int> route_seq;
        cur_solution_cost -= routes[current_route.route_id]->length;
        total_vio_time -= mcgrp.get_vio_time(current_route.time_tbl);
        total_vio_load -= max(routes[current_route.route_id]->load - mcgrp.capacity,0);

        for (const auto &node : current_route.time_tbl)
            route_seq.push_back(node.task);

        delete_route(current_route.route_id, route_seq);
        task_list.insert(task_list.end(),route_seq.begin(),route_seq.end());
        My_Assert(valid_sol(mcgrp), "Wrong validation");
    }

//    Individual giant_tour = RTF(mcgrp,task_list,true);
//    auto new_routes = tour_splitting(mcgrp,giant_tour.sequence);
    auto new_routes = tour_splitting(mcgrp,task_list);

    for (const auto& route : new_routes){
        auto new_route_id = new_route(mcgrp, route);
        My_Assert(new_route_id != -1, "Wrong new route");
        cur_solution_cost += routes[new_route_id]->length;
    }

    My_Assert(missed(mcgrp), "Some Task missed!");
    My_Assert(check_duplicated(mcgrp), "Duplicated Task!");
    My_Assert(total_vio_time == 0, "This is not a infeasible Task!");

    My_Assert(valid_sol(mcgrp), "Repair method doesn't work properly!");
}

vector<int>
HighSpeedNeighBorSearch::get_successor_tasks(int length, const int chosen_task)
{
    My_Assert(chosen_task != DUMMY, "Chosen task can't be dummy");

    int current_node;
    vector<int> buffer;

    current_node = chosen_task;
    buffer.push_back(current_node);
    while (buffer.size() < length) {
        current_node = solution[current_node]->next->ID;

        if (current_node < 0) {
            return buffer;
        }
        else {
            buffer.push_back(current_node);
        }
    }

    return buffer;
}

pair<struct HighSpeedNeighBorSearch::TASK_NODE *, struct HighSpeedNeighBorSearch::TASK_NODE *>
HighSpeedNeighBorSearch::SOLUTION::connect_tasks(const vector<int> &route_seq)
{
    if (route_seq.empty())
        return {nullptr, nullptr};

    if (route_seq.size() == 1)
        return {&tasks[route_seq[0]], &tasks[route_seq[0]]};

    for (int i = 0; i < route_seq.size(); i++) {
        if (i == 0) {
            tasks[route_seq[i]].next = &tasks[route_seq[i + 1]];
        }
        else if (i == route_seq.size() - 1) {
            tasks[route_seq[i]].pre = &tasks[route_seq[i - 1]];
        }
        else {
            tasks[route_seq[i]].next = &tasks[route_seq[i + 1]];
            tasks[route_seq[i]].pre = &tasks[route_seq[i - 1]];
        }
    }

    return {&tasks[route_seq.front()], &tasks[route_seq.back()]};
}
