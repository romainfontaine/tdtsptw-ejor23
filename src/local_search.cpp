#include "local_search.h"


uint LocalSearch::search(const static_array<uint, MAX_SIZE> &base_solution, static_array<uint, MAX_SIZE> &path,
                         const uint &cost, const real_time_point &tp) {
    n_neighbors_considered = 0;
    n_neighbors_evaluated = 0;

    auto t_begin = std::chrono::steady_clock::now();
    uint obj = cost;
    path = base_solution;

    static_array<uint, MAX_SIZE> new_path;
    uint new_obj;
    bool changes(true);
    uint it(0);
    while (changes) {
        it += 1;
        changes = false;
        if (ls_1shift(path, obj, new_path, new_obj)) {
            obj = new_obj;
            path = new_path;
            changes = true;
            log_improvement(obj, it, tp, "ls_1shift"s);
        }
        if (ls_2opt(path, obj, new_path, new_obj)) {
            obj = new_obj;
            path = new_path;
            changes = true;
            log_improvement(obj, it, tp, "ls_2opt"s);
        }
    }
    auto elapsed = chrono::duration_cast<chrono::duration<double >>(
            std::chrono::steady_clock::now() - (t_begin)).count();
    cout << "# Local search time: " << elapsed << "s. ";
    cout << "(" << n_neighbors_evaluated << "/" << n_neighbors_considered << ", " << it << " IT)" << endl;
    return obj;
}


bool LocalSearch::ls_check_path(const static_array<uint, MAX_SIZE> &path, uint &obj_out) {
    State s = tsp.initial_state();
    for (uint i = 1; i < path.size(); i++) {
        const uint &a = path[i];
        if (!tsp.A(s).contains(a)) {
            return false;
        }
        s = tsp.tau(s, a);
    }
    if (!tsp.is_terminal(s)) {
        cerr << "Invalid path..." << endl;
        exit(1);
    }
    obj_out = s.t;
    return true;
}

bool LocalSearch::ls_1shift_move(const static_array<uint, MAX_SIZE> &path, const uint &obj,
                                 static_array<uint, MAX_SIZE> &new_path, uint &return_obj,
                                 const uint &f, const uint &t) {

    if (f == t)
        return false;
    if (t == f - 1) // 1 shift of two adjacent nodes in the original path is symmetric
        //               (keep t==f+1, skip t==f-1 - i.e. only shift forward when |f-t|=1)
        return false;
    new_path.clear();
    n_neighbors_considered++;
    /*
     * Edge cases
     *   F             T        _             _     //     f-1|f  f|f+1   t|t+1       f-1|f+1 t|f     f|t+1
     * 0 1 2 3 4 5 6 7 8 0    0 2 3 4 5 6 7 8 1 0   // cut 01     12      80      add 02      81      10
     *
     *   F T                    _ _                 //     f-1|f  f|f+1   t|t+1       f-1|f+1 t|f     f|t+1
     * 0 1 2 3 4 5 6 7 8 0    0 2 1 3 4 5 6 7 8 0   // cut 01     12      23      add 02      21      13
     *
     *                                              //                    *****               ***     *****
     *   T             F        _             _     //     f-1|f  f|f+1   t-1|t       f-1|f+1 f|t     t-1|f
     * 0 1 2 3 4 5 6 7 8 0    0 8 1 2 3 4 5 6 7 0   // cut 78     80      01      add 70      81      08
     *
     *               T F                    _ _     //     f-1|f  f|f+1   t-1|t       f-1|f+1 f|t     t-1|f
     * 0 1 2 3 4 5 6 7 8 0    0 1 2 3 4 5 6 8 7 0   // cut 78     80      67      add 70      87      68
     *
     */
    if (f < t) {
        /*
         * Move the node later in the path
         *
         * Before
         *     F       T
         * 0 1 2 3 4 5 6 7 8 9
         *
         * After
         *     _       _
         * 0 1 3 4 5 6 2 7 8 9
         *
         * Valid heuristic for the TSP but not for the TDTSP :
         *   - The sum of the costs of the removed arcs must be higher than the
         *     sum of the costs of the added arcs
         *   - Not valid for the TDTSP because the cost of the path between the nodes may change
         *
         * Valid heuristics :
         *   - Arc (T,F) must be usable
         *   - Arc (F-1,F+1) must be usable
         *   - Arc (F,T+1) must be usable
         *   - F must not precede any node in the interval [F+1; T]
         *
         * Then, if all conditions are met, compute the cost and check the validity of the new solution
         */
        if (LocalSearch::FILTERING) {
            if (!tsp.E[path[t]].contains(path[f]))            // Arc (T,F) must be usable
                return false;
            if (!tsp.E[path[f - 1]].contains(path[f + 1]))  // Arc (F-1,F+1) must be usable
                return false;
            if (!tsp.E[path[f]].contains(path[t + 1]))      // Arc (F,T+1) must be usable
                return false;
            Bitset new_preceding_nodes;
            for (uint i = f + 1; i <= t; i++) {
                new_preceding_nodes.add(path[i]);
            }
            if (!(tsp.R[path[f]] & new_preceding_nodes).is_empty())
                return false;
        }


        for (uint i = 0; i < f; i++)
            new_path.push_back(path[i]);
        for (uint i = f + 1; i <= t; i++)
            new_path.push_back(path[i]);
        new_path.push_back(path[f]);
        for (uint i = t + 1; i < path.size(); i++)
            new_path.push_back(path[i]);
    } else {
        /*
         * Move the node earlier in the path
         *
         * Before
         *     T       F
         * 0 1 2 3 4 5 6 7 8 9
         *
         * After
         *     _       _
         * 0 1 6 2 3 4 5 7 8 9
         *
         * Valid heuristics :
         *   - Arc (F,T) must be usable
         *   - Arc (F-1,F+1) must be usable
         *   - Arc (T-1,F) must be usable
         *   - No node in the interval [T;F-1] must precede F
         */
        if (LocalSearch::FILTERING) {
            if (!tsp.E[path[f]].contains(path[t]))            // Arc (F,T) must be usable
                return false;
            if (!tsp.E[path[f - 1]].contains(path[f + 1]))  // Arc (F-1,F+1) must be usable
                return false;
            if (!tsp.E[path[t - 1]].contains(path[f]))      // Arc (T-1,F) must be usable
                return false;

            Bitset precedence_constraints;
            for (uint i = t; i < f; i++)
                precedence_constraints |= tsp.R[path[i]];

            if (precedence_constraints.contains(path[f]))
                return false;
        }


        for (uint i = 0; i < t; i++)
            new_path.push_back(path[i]);
        new_path.push_back(path[f]);
        for (uint i = t; i < f; i++)
            new_path.push_back(path[i]);
        for (uint i = f + 1; i < path.size(); i++)
            new_path.push_back(path[i]);
    }

    uint new_obj;
    n_neighbors_evaluated++;
    if (!ls_check_path(new_path, new_obj))
        return false;

    return_obj = min(obj, new_obj);
    return new_obj < obj;
}


bool LocalSearch::ls_1shift(const static_array<uint, MAX_SIZE> &path, const uint &obj,
                            static_array<uint, MAX_SIZE> &new_path, uint &return_obj) {

    /*
     * X = previous ending point (included in step 1)
        3 3 3 3 3 3 3 3 3 3
        3 3 3 3 3 3 3 3 3 3
        3 3 3 3 3 3 3 3 3 3
        4 4 4 4 4 4 4 X 1 1
        2 2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2 2
     */

    uint orig_f(one_shift_f), orig_t(one_shift_t);

    for (; one_shift_t < path.size() - 1; one_shift_t++) // Step 1
        if (ls_1shift_move(path, obj, new_path, return_obj, one_shift_f, one_shift_t))
            return true;

    for (one_shift_f = orig_f + 1; one_shift_f < path.size() - 1; one_shift_f++) // Step 2
        for (one_shift_t = 1; one_shift_t < path.size() - 1; one_shift_t++)
            if (ls_1shift_move(path, obj, new_path, return_obj, one_shift_f, one_shift_t))
                return true;

    for (one_shift_f = 1; one_shift_f < orig_f; one_shift_f++) // Step 3
        for (one_shift_t = 1; one_shift_t < path.size() - 1; one_shift_t++)
            if (ls_1shift_move(path, obj, new_path, return_obj, one_shift_f, one_shift_t))
                return true;

    one_shift_f = orig_f;
    for (one_shift_t = 1; one_shift_t < orig_t; one_shift_t++) // Step 4
        if (ls_1shift_move(path, obj, new_path, return_obj, one_shift_f, one_shift_t))
            return true;
    return false;
}

bool LocalSearch::ls_2opt_move(const static_array<uint, MAX_SIZE> &path,
                               const uint &obj, static_array<uint, MAX_SIZE> &new_path,
                               uint &return_obj, const uint &i, const uint &k) {
    new_path.clear();
    n_neighbors_considered++;
    /*
     * Examples :
     *   i             k        _             _
     * 0 1 2 3 4 5 6 7 8 0    0 8 7 6 5 4 3 2 1 0
     *   i     k                _     _
     * 0 1 2 3 4 5 6 7 8 0    0 4 3 2 1 5 6 7 8 0
     *   i   k                  _   _
     * 0 1 2 3 4 5 6 7 8 0    0 3 2 1 4 5 6 7 8 0
     *   i k                    _ _
     * 0 1 2 3 4 5 6 7 8 0    0 2 1 3 4 5 6 7 8 0
     *
     * Necessary conditions for a move to be valid (not sufficient) :
     * - Arc (i-1, k) can be used
     * - Arc (i, k+1) can be used
     * - For all j in [i, k[ :
     *   - arc (j+1, j) can be used
     *
     */
    if (LocalSearch::FILTERING) {
        if (!tsp.E[path[i - 1]].contains(path[k]))
            return false;
        if (!tsp.E[path[i]].contains(path[k + 1]))
            return false;
        for (uint j = i; j < k; j++)
            if (!tsp.E[path[j + 1]].contains(path[j]))
                return false;
        /*
         *   i             k        _             _
         * 0 1 2 3 4 5 6 7 8 0    0 8 7 6 5 4 3 2 1 0
         *
         *
         * 1 cannot precede any of [i;k]
         * 2 cannot precede any of [i+1;k]
         * 3 cannot precede any of [i+2;k]
         * ...
         * k-1 cannot precede k
         */
        Bitset interval;
        for (uint j = i; j <= k; j++)
            interval.add(path[j]);

        for (uint j = i; j < k; j++) {
            interval.remove(path[j]);
            if (!(tsp.R[path[j]] & interval).is_empty())
                return false;
        }
    }
    for (uint j = 0; j < i; j++)
        new_path.push_back(path[j]);
    for (uint j = k; j >= i; j--)
        new_path.push_back(path[j]);
    for (uint j = k + 1; j < path.size(); j++)
        new_path.push_back(path[j]);


    uint new_obj;
    n_neighbors_evaluated++;
    if (!ls_check_path(new_path, new_obj))
        return false;

    return_obj = min(obj, new_obj);
    return new_obj < obj; // found improvement
}

bool
LocalSearch::ls_2opt(const static_array<uint, MAX_SIZE> &path,
                     const uint &obj, static_array<uint, MAX_SIZE> &new_path,
                     uint &return_obj) {
    uint orig_i(two_opt_i), orig_k(two_opt_k);

    // Step 1
    for (; two_opt_k < path.size() - 1; two_opt_k++) {
        if (ls_2opt_move(path, obj, new_path, return_obj, two_opt_i, two_opt_k))
            return true;

    }

    // Step 2
    for (two_opt_i = orig_i + 1; two_opt_i < path.size() - 1; two_opt_i++) {
        for (two_opt_k = two_opt_i + 2; two_opt_k < path.size() - 1; two_opt_k++) {
            // + 2 because swap is included in 1shift
            if (ls_2opt_move(path, obj, new_path, return_obj, two_opt_i, two_opt_k))
                return true;

        }
    }

    // Step 3
    for (two_opt_i = 1; two_opt_i < orig_i; two_opt_i++) {
        for (two_opt_k = two_opt_i + 2;
             two_opt_k < path.size() - 1; two_opt_k++) { // + 2 because swap is included in 1shift
            if (ls_2opt_move(path, obj, new_path, return_obj, two_opt_i, two_opt_k))
                return true;

        }
    }

    // Step 4
    two_opt_i = orig_i;
    for (two_opt_k = two_opt_i + 2; two_opt_k < orig_k; two_opt_k++) { // + 2 because swap is included in 1shift
        if (ls_2opt_move(path, obj, new_path, return_obj, two_opt_i, two_opt_k))
            return true;

    }
    return false;
}

void
LocalSearch::log_improvement(const uint &obj, const uint &it, const real_time_point &tp, const string &neighborhood) {

    auto elapsed = chrono::duration_cast<chrono::duration<double >>(
            std::chrono::steady_clock::now() - tp).count();
    const auto t = CPUTime();
    cout << "(" << obj << ", " << it << ", " << n_neighbors_evaluated << ", (" << elapsed << ", " << t.first << ", "
         << t.second << "), '" << neighborhood << "')" << endl;
}

