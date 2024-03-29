#include <array>
#include <iostream>
#include <boost/filesystem.hpp>
#include "RNG.h"
#include "utils.h"
#include "local_search.h"

using namespace std;

ofstream result_out;
ofstream log_out;

system_clock::time_point cur_time;
system_clock::time_point search_start_time;

void __My_Assert(const char *expr_str, bool expr, const char *file, int line, const char *msg)
{
    if (!expr) {
        std::cerr << "Assert failed:\t" << msg << "\n"
                  << "Expected:\t" << expr_str << "\n"
                  << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

vector<int> get_negative_coding(const vector<int> &sequence)
{
    int sign = 1;
    vector<int> negative_coding_sol;
    negative_coding_sol.clear();

    for (auto task : sequence) {
        if (task == DUMMY) {
            sign = -1;
        }
        else {
            negative_coding_sol.push_back(sign * task);
            sign = 1;
        }
    }

    return negative_coding_sol;
}

vector<int> get_delimiter_coding(const vector<int> &negative_coding)
{
    vector<int> delimiter_coding;
    delimiter_coding.clear();
    for (auto task:negative_coding) {
        if (task < 0) {
            delimiter_coding.push_back(DUMMY);
        }
        delimiter_coding.push_back(abs(task));
    }
    delimiter_coding.push_back(DUMMY);

    return delimiter_coding;
}

std::vector<std::string> read_directory(const std::string dir_path)
{
    cout<<"Reading directory...\n";
    using namespace boost::filesystem;
    std::vector<std::string> buffer;

    path p(dir_path);
    std::vector<directory_entry> v;
    if(is_directory(p)){
        copy(directory_iterator(p),directory_iterator(),back_inserter(v));
        for(std::vector<directory_entry>::const_iterator it = v.begin(); it != v.end(); ++it){
            if(is_regular_file(*it))
            {
                buffer.push_back((*it).path().filename().string());
            }
        }
    }

    return buffer;
}

bool print(std::ostream &os1, std::ostream &os2, const std::string &str)
{
    if(os1 << str << endl && os2 << str << endl){
        return true;
    }

    return false;
}

bool print(std::ostream &os1, const std::string &str)
{
    if(os1 << str << endl){
        return true;
    }

    return false;
}

//double sum(const vector<double> &vec)
//{
//    double tmp = 0;
//    for (auto &n: vec) {
//        tmp += n;
//    }
//    return tmp;
//}

vector<int> find_ele_positions(const std::vector<int> &seq, int target)
{
    vector<int> buffer;
    for (int i = 0; i < seq.size(); i++) {
        if (seq[i] == target) {
            buffer.push_back(i);
        }
    }

    return buffer;
}

double sel_ratio(const vector<double> &_prob, const vector<double> &ratio, RNG &rng)
{
    double rd = rng.Randfloat(0, 1);
    double fsum = 0;
    for (int i = 0; i < _prob.size(); ++i) {
        fsum += _prob[i];
        if (rd < fsum) {
            return ratio[i];
        }
    }

    My_Assert(false, "ERROR! Cannot select a ratio!");
}

vector<int> sort_solution(const vector<int>& negative_sol)
{
    map<int, vector<int>> routes;
    vector<int> route;
    for(const int task : negative_sol){
        if(task < 0){
            if(task != negative_sol.front()){
                routes[route[0]] = route;
                route.clear();
            }

            route.push_back(-task);
        }else{
            route.push_back(task);
        }
    }

    if(!route.empty())
        routes[route[0]] = route;

    vector<int> sorted;
    for(auto& item : routes){
        item.second[0] *= -1;
        sorted.insert(sorted.cend(), item.second.cbegin(), item.second.cend());
    }

    return sorted;
}

vector<double> stable_softmax(const vector<double> &x)
{
    double max_item = *max_element(x.begin(),x.end());
    vector<double> numerator = x;
    for(auto &item: numerator) {
        item -= max_item;
        item = exp(item);
    }

    double denominator = accumulate(numerator.begin(),numerator.end(),0);

    for(auto &item: numerator) {
        item /= denominator;
    }

    return numerator;
}
vector<double> stable_uniform(const vector<double> &x)
{
    vector<double> numerator = x;

    double denominator = accumulate(numerator.begin(),numerator.end(),0.0);
    if(denominator == 0) return numerator;

    for(auto &item: numerator) {
        item /= denominator;
    }

    return numerator;
}

double get_time_difference(const system_clock::time_point &start, const system_clock::time_point &end){
    auto duration = duration_cast<microseconds>(end - start);
    return double(duration.count()) * microseconds::period::num / microseconds::period::den;
}

vector<RouteInfo::TimeTable>
RouteInfo::TimeTable::zip(const vector<int> &task_list, const vector<int> &ArriveTime)
{
    My_Assert(task_list.size() == ArriveTime.size(), "sequence length does not same");
    vector<RouteInfo::TimeTable> ans;
    for(int i = 0; i< task_list.size();i++){
        ans.push_back({task_list[i],ArriveTime[i]});
    }

    return ans;
}

vector<vector<int>>
RouteInfo::TimeTable::unzip(const vector<TimeTable> &time_tbl)
{
    if(time_tbl.empty()) {
        return {};
    }

    vector<vector<int>> ans(2,vector<int>());
    for(const auto& tb : time_tbl){
        ans[0].push_back(tb.task);
        ans[1].push_back(tb.arrive_time);
    }

    return ans;
}

RouteScores::RouteScores(int routeId, int stableScore)
    : route_id(routeId), scores(stableScore)
{}
