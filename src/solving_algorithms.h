#ifndef TDTSPTW_SOLVING_ALGORITHMS_H
#define TDTSPTW_SOLVING_ALGORITHMS_H

#include "constants.h"
#include "util.hpp"
#include "unordered_dense.h"
#include "model.h"
#include <tuple>
#include <optional>

template<typename T>
void anytime_log(const string &prefix, const bool &found_solution, const uint &best,
                 const uint &iteration, const T &opened, const real_time_point &starting_time) {
    auto elapsed = chrono::duration_cast<chrono::duration<double >>(
            std::chrono::steady_clock::now() - starting_time).count();
    // If no solution has been found, the objective is set to None
    const string str_best = found_solution ? to_string(best) : "None";
    const auto t = CPUTime();
    cout << prefix << "(" << str_best << ", " << iteration << ", " << opened << ", (" << elapsed << ", "
         << t.first << ", " << t.second << "))" << endl;
}

void print_path(const Path &p);

bool check_path(const Path &path, const TDTSPTW &t, const uint &best);

template<typename T>
bool compute_path_from_nd(const vector<T> &nd, const TDTSPTW &tsp, Path &path) {
    return compute_path_from_nd(nd, tsp, path, [](const State &s) -> const State & { return s; });
}

template<typename T, typename F>
bool compute_path_from_nd(const vector<T> &nd, const TDTSPTW &tsp, Path &path, const F &f) {
    path.resize(tsp.N);
    path[0] = tsp.start;
    uint level = nd.size() - 1;

    State s = f(*nd[level].begin());
    // while not initial, try each (reverse) action (do not check if it is feasible or not - if infeasible, not found!)
    while (level != 0) {
        Bitset possible_actions = s.S;
        bool found(false);
        for (const auto &a: possible_actions) {
            State candidate(a, Bitset(s.S).remove(a));
            const auto &it = nd[level - 1].find(candidate);
            if (it != nd[level - 1].end()) {
                // state s_ may lead to state s
                const State &s_ = f(*it);

                // check the fact that state s_ actually leads to state s (by applying action s.i)
                // Using the original TWs (as "preprocess_during" may increase tw[x].early and prevent
                //  from retrieving the path)
                if (tsp.A<true>(s_).contains(s.i) && tsp.tau<true>(s_, s.i).t <= s.t) {
                    found = true;
                    path[level] = s.i;
                    // Go on to the previous level, looking for s_'s predecessor
                    level--;
                    s = s_;
                    break;
                }
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

using GrowthFn = function<uint(const uint &)>;

std::tuple<bool, bool, uint, Path> // optimality proof, found_solution, objective, path
ACS(const Heuristic &IH, TDTSPTW &t, const uint &w,
    const bool &greedy_ub, const bool &gub_ls_twp, const bool &local_search, const bool &preprocess_tws,
    const bool &process_tws, const uint &tl, const bool &fewer_tl_checks, const bool &free_mem,
    const real_time_point &starting_time);

template<bool use_fea_heuristic = false>
tuple<bool, bool, uint, Path> bottom_up(TDTSPTW &t, const bool &preprocess_tws, const int initial_ub,
                                        const real_time_point &starting_time, const bool compute_path = true) {
    const uint N = t.N;
    vector<ankerl::unordered_dense::set<State>> nd(N);
    static_array<State, MAX_SIZE> children;
    static_array<uint, MAX_SIZE> n_nd(N, 0);

    nd[0].insert(t.initial_state());
    n_nd[0] = 1;

    bool infeasible = false;
    if (initial_ub == -1)
        infeasible = t.process_tws_and_compute_lbs(preprocess_tws, use_fea_heuristic);
    else
        t.process_tws_and_compute_lbs(initial_ub);

    const uint ub = t.tws[t.end].late + 1;

    if (!infeasible)
        for (uint l = 0; l < N - 1; l++) {
            const auto &cur_nd = nd[l];
            auto &next_nd = nd[l + 1];
            for (const State &s: cur_nd) {
                children.clear();
                t.A(s, children);
                for (const auto &s_: children) {
                    const auto &best_it = next_nd.find(s_);
                    if (best_it != next_nd.end()) {
                        if (best_it->t > s_.t)
                            best_it->t = s_.t;
                    } else if (!use_fea_heuristic ||
                               t.h_0_out_feasibility_checks<false, LBCostType::medium>(s_, 0) != ub) {
                        next_nd.insert(s_);
                    }
                }
            }
            n_nd[l + 1] = nd[l + 1].size();
            if (!compute_path)
                nd[l] = ankerl::unordered_dense::set<State>(); // clear & free memory

            const auto elapsed = chrono::duration_cast<chrono::duration<double>>(
                    std::chrono::steady_clock::now() - starting_time).count();
            cout << "# (" << l << ", " << n_nd[l] << ", " << elapsed << ")" << endl;
        }
    const auto &it(nd[N - 1].begin());
    const bool found_solution(it != nd[N - 1].end());
    const uint best = found_solution ? it->t : ~0u;

    anytime_log("Complete", found_solution, best, 0, 0, starting_time);
    Path path;
    if (found_solution) {
        if (compute_path && compute_path_from_nd(nd, t, path)) {
            print_path(path);
            if (!check_path(path, t, best)) {
                cout << "Path is NOT valid!" << endl;
                exit(1);
            }
        } else {
            cout << "# Path : None" << endl;
        }
    }
    uint acc(0);
    cout << "# Number of non-dominated states : " << endl << "# ";
    for (const auto &n: n_nd) {
        cout << setw(8) << n << " ";
        acc += n;
    }
    cout << endl << "# " << acc << " total." << endl;

    cout << "# " << peak_memory() << endl;

    return make_tuple(true, found_solution, best, path);
}

#endif //TDTSPTW_SOLVING_ALGORITHMS_H
