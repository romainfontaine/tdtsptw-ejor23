#include <iomanip>
#include <set>
#include <map>
#include <unistd.h> // for getpid()
#include <thread>

#include "solving_algorithms.h"
#include "local_search.h"
#include "scc.hpp"

using namespace std;

bool check_path(const Path &path, const TDTSPTW &t, const uint &best) {
    State s = t.initial_state();
    for (uint i = 1; i < path.size(); i++) {
        if (!t.A(s).contains(path[i]))
            return false;
        s = t.tau(s, path[i]);
    }
    return s.t == best;
}

Path greedy_path(const TDTSPTW &tsp, bool &success) {
    Path path;
    path.push_back(tsp.start);
    success = true;
    State s = tsp.initial_state();
    while (!tsp.is_terminal(s)) {
        if (tsp.A(s).is_empty()) {
            success = false;
            return path;
        }
        path.push_back(tsp.greedy_action(s));
        s = tsp.tau(s, path.back());
    }
    return path;
}

uint compute_path_cost(const Path &path, const TDTSPTW &tsp) {
    State s = tsp.initial_state();
    for (uint i = 1; i < path.size(); i++) {
        s = tsp.tau(s, path[i]);
    }
    return s.t;
}

void print_path(const Path &p) {
    cout << "# Path : [";
    for (const Action &a: p)
        cout << a << ",";
    cout << "]" << endl;
}

struct ACSCompareState {
    bool operator()(const pair<uint, State> &lhs, const pair<uint, State> &rhs) {
        return lhs.first > rhs.first;
    }
};

tuple<uint, optional<Path>, bool> // obj, path, optimality_proof
greedy_upper_bound(TDTSPTW &t, const bool &infeasible, const bool &compute_greedy_ub, const bool &local_search,
                   const bool &process_tws, LocalSearch &ls, const real_time_point &starting_time) {
    const uint ub = t.upper_bound_from_tw();
    if (!compute_greedy_ub || infeasible)
        return {ub, {}, infeasible};

    uint obj = min(t.greedy_upper_bound(t.initial_state()), ub);
    bool found_solution(obj != ub);

    anytime_log(found_solution ? "" : "#", true, obj, 0, 0, starting_time);
    if (!found_solution)
        return {ub, {}, false};

    Path path;
    bool success;
    path = greedy_path(t, success);
    if (!success) {
        cerr << "Failed to compute greedy path" << endl;
        exit(1);
    }

    if (local_search) {
        Path new_path;
        obj = ls.search(path, new_path, obj, starting_time);
        path = new_path;
    }

    const bool optimality_proof = process_tws && t.process_tws_and_compute_lbs(obj);
    print_path(path);
    return {obj, path, optimality_proof};
}

tuple<bool, bool, uint, Path> // optimality proof, found_solution, objective, path
ACS(const Heuristic &IH, TDTSPTW &t, const uint &w,
    const bool &greedy_ub, const bool &gub_ls_twp, const bool &local_search, const bool &preprocess_tws,
    const bool &process_tws, const uint &tl, const bool &fewer_tl_checks, const bool &free_mem,
    const real_time_point &starting_time) {
    using PriorityQueues = vector<priority_queue<pair<uint, State>, vector<pair<uint, State >>, ACSCompareState> >;
    using DominanceSets = vector<ankerl::unordered_dense::set<State>>;

    const bool infeasible(t.process_tws_and_compute_lbs(preprocess_tws, true));  // DO compute cost LBs

    LocalSearch ls(t);
    uint N = t.N - 1;
    // Trick to avoid losing time freeing the memory of these (numerous) objects
    //  - If free_mem is false, memory is freed instantly by the OS when the program exits
    //  - otherwise, every destructor is called (which may take ~20sec for 10Gb of RAM usage)
    auto *queues_ = new PriorityQueues(N);
    auto *nd_ = new DominanceSets(N + 1);
    auto &queues = *queues_;
    auto &nd = *nd_;

    const auto tl_end = starting_time + chrono::seconds(tl);

    auto [best, p, optimality_proof] = greedy_upper_bound(t, infeasible, greedy_ub, local_search | gub_ls_twp,
                                                          process_tws | gub_ls_twp, ls, starting_time);
    bool found_solution(p.has_value());

    Path path, new_path;
    {
        const State &s_0 = t.initial_state();
        queues[0].emplace(0, s_0);
        // invariant : a state present in queue[x] is also present in nd[x]
        // (although it may have a better t value in nd[x])
        nd[0].insert(s_0);
    }
    static_array<State, MAX_SIZE> children;
    static_array<uint, MAX_SIZE> h_values;
    uint iteration = 0;
    uint opened = 0;
    uint first_non_empty = 0; // level of the first non-empty queue

    static_array<uint, MAX_SIZE> lb_per_level(N + 1, ~0u);
    uint best_lb = t.tws[t.end].early;
    
    while (!optimality_proof && std::chrono::steady_clock::now() < tl_end) { // stop on exhaustive search OR time limit
        ++iteration;
        for (uint level = first_non_empty;
             level < N && (fewer_tl_checks || std::chrono::steady_clock::now() < tl_end); level++) {
            const uint next_level = level + 1;
            for (uint j = 0; j < w && !queues[level].empty();) {
                {
                    const auto &[f, s] = queues[level].top();
                    // Discard s if : - it is bounded (i.e. it has no improving completion)
                    //                - or it is dominated (i.e. we know a better t for this pair (i,S))
                    if (f >= best || (s.t > nd[level].find(s)->t)) { // s must exist in nd[level]
                        queues[level].pop();
                        continue; // (without incrementing j)
                    }

                    lb_per_level[level] = min(f, lb_per_level[level]);

                    ++opened;
                    ++j;

                    children.clear();
                    h_values.clear();
                    t.A(s, children);
                    const auto &it = remove_if(children.begin(), children.end(),
                                               [&dom = as_const(nd[next_level]),
                                                       &best = as_const(best)]
                                                       (const State &s_) -> bool {
                                                   if (s_.t >= best)
                                                       return true;
                                                   const auto &it(dom.find(s_));
                                                   return it != dom.end() && s_.t >= it->t;
                                               });
                    children.resize(distance(children.begin(), it));
                    (t.*IH)(s, children, h_values); // Compute h-values
                    queues[level].pop(); // s is not needed anymore.
                }

                uint i = 0;
                for (const auto &s_: children) { // For each non-dominated child
                    const uint &h_ = h_values[i++];
                    const uint f_ = s_.t + h_;
                    // s_=(i,S,t) is not dominated ; there exist no state s'=(i,S,t') in queue such that t'<=t.
                    // therefore, there is *at most* one state in queue with the best-known t value for (i,S)
                    if (f_ < best) { // bound
                        // Update or insert
                        auto x = nd[next_level].find(s_);
                        if (x != nd[next_level].end())
                            x->t = s_.t;
                        else
                            nd[next_level].insert(s_);

                        if (next_level != N)
                            queues[next_level].emplace(f_, s_);
                        else {
                            best = s_.t;
                            if (!compute_path_from_nd(nd, t, path)) {
                                cerr << "Error while computing ACS's path" << endl;
                                exit(1);
                            }
                            const uint new_cost = compute_path_cost(path, t);
                            if (new_cost < best) {
                                // ACS may have partially found an improvement, but did not have enough time
                                // to propagate it to the final level
                                best = new_cost;
                            }
                            if (!check_path(path, t, best)) {
                                cerr << "Error while checking ACS's path" << endl;
                                exit(1);
                            }
                            anytime_log("", true, best, iteration, opened, starting_time);
                            found_solution = true;
                            if (local_search) {
                                best = ls.search(path, new_path, best, starting_time);
                                path = new_path;
                            }
                            const bool &opt_proven = process_tws && t.process_tws_and_compute_lbs(best);
                            print_path(path);
                            if (opt_proven) {
                                // set conditions to consider the instance as solved and to stop exploration:
                                first_non_empty = N; // i.e. optimality_proof = true
                                j = w; // terminate w loop
                                break;
                            }
                        }
                    }
                }
            }
        }
        while (first_non_empty < N && queues[first_non_empty].empty()) {
            shrinkToFit(queues[first_non_empty]); // first non-empty is empty, clear it.
            ++first_non_empty;
        }

        if (best_lb < lb_per_level[first_non_empty - 1]) {
            best_lb = lb_per_level[first_non_empty - 1];
            if (best_lb == best)
                first_non_empty = N; // i.e. optimality_proof = true
        }

        optimality_proof = first_non_empty == N; // if queue[N-1] is empty, queue[N] is necessarily empty.
    }

    anytime_log(optimality_proof ? "Complete"s : "tl"s, found_solution, best, iteration, opened, starting_time);

    {
        uint acc(0);
        cout << "# Remaining states in queues : " << endl << "# ";
        for (const auto &q: queues) {
            cout << setw(8) << q.size() << " ";
            acc += q.size();
        }
        cout << endl << "# " << acc << " total." << endl;
    }

    {
        uint acc(0);
        cout << "# Number of non-dominated states : " << endl << "# ";
        for (const auto &m: nd) {
            cout << setw(8) << m.size() << " ";
            acc += m.size();
        }
        cout << endl << "# " << acc << " total." << endl;
    }

    cout << "# " << peak_memory() << endl;

    if (free_mem) {
        delete queues_;
        delete nd_;
    }
    return make_tuple(optimality_proof, found_solution, best, path);
}
