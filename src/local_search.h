#ifndef TDTSPTW_LOCAL_SEARCH_H
#define TDTSPTW_LOCAL_SEARCH_H

#include <chrono>
#include "model.h"

class LocalSearch {
    const TDTSPTW &tsp;
    const bool FILTERING = true;
    uint n_neighbors_considered = 0;
    uint n_neighbors_evaluated = 0;

    uint two_opt_i = 1;
    uint two_opt_k = 3;

    uint one_shift_f = 1;
    uint one_shift_t = 1;

public:
    LocalSearch(const TDTSPTW &tdtsptw) : tsp(tdtsptw) {};

    uint search(const static_array<uint, MAX_SIZE> &base_solution, static_array<uint, MAX_SIZE> &path,
                const uint &cost, const real_time_point &tp);

private:
    void log_improvement(const uint &obj, const uint &it, const real_time_point &tp, const string &neighborhood);

    bool ls_check_path(const static_array<uint, MAX_SIZE> &path, uint &obj_out);

    bool ls_1shift(const static_array<uint, MAX_SIZE> &path, const uint &obj,
                   static_array<uint, MAX_SIZE> &new_path, uint &return_obj);

    bool ls_2opt(const static_array<uint, MAX_SIZE> &path, const uint &obj,
                 static_array<uint, MAX_SIZE> &new_path, uint &return_obj);

    bool ls_2opt_move(const static_array<uint, MAX_SIZE> &path, const uint &obj, static_array<uint, MAX_SIZE> &new_path,
                      uint &return_obj, const uint &i, const uint &k);

    bool
    ls_1shift_move(const static_array<uint, MAX_SIZE> &path, const uint &obj, static_array<uint, MAX_SIZE> &new_path,
                   uint &return_obj, const uint &f, const uint &t);
};


#endif //TDTSPTW_LOCAL_SEARCH_H
