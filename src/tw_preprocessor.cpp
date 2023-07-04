#include "tw_preprocessor.h"
#include "model.h"
#include <algorithm>
#include <cassert>

using namespace std;

TWPreprocessor::TWPreprocessor(TDTSPTW *tsp, const bool &triangle_inequality_assumed,
                               const bool &propagate_precedence_only,
                               const vector<pair<uint, uint>> &precedence_constraints)
        : tsp(tsp),
          N(tsp->N),
          E(tsp->E),
          E_(tsp->E_),
          R(tsp->R),
          R_(tsp->R_),
          e(N, 0),
          l(N, 0),
          s(N, 0),
          start(tsp->start),
          end(tsp->end),
          triangle_inequality_assumed(triangle_inequality_assumed),
          infeasible_triples_mem(N * N),
          propagate_precedence_only(propagate_precedence_only) {
    for (uint i = 0; i < N; i++) {
        e[i] = tsp->tws[i].early;
        l[i] = tsp->tws[i].late;
        s[i] = tsp->service_times[i];
    }
    for (uint i = 0; i < N; i++)
        removeArc(i, i);

    for (const auto &[i, j]: precedence_constraints)
        add_precedence(i, j);

    for (uint i = 0; i < N; i++) {
        if (i != start)
            add_precedence(start, i);
        if (i != end)
            add_precedence(i, end);
    }
    // Naive preprocessing in O(n^2):
    // - no inference of precedence constraints,
    // - no tw tightening.
    // Required to correctly compute Latest Departure Times.
    for (uint i = 0; i < N; i++)
        for (uint j = 0; j < N; j++)
            if (E[i][j])
                if (delta(i, j, e[i] + s[i]) > l[j] - s[j])
                    removeArc(i, j);
}

/*
 * Tightens TWs, infer precedence & filter arcs
 * Returns true if the instance is infeasible
 */
bool TWPreprocessor::preprocess() {
    if (propagate_precedence_only) {
        bool changes = true;
        for (it = 0; changes; it++) {
            changes = false;
            changes = transitiveClosurePrecedence() || changes;
            changes = removeArcsFromPrecedencePairs() || changes;
        }
        return false;
    }
    bool changes = true;
    for (it = 0; changes; it++) {
        changes = false;
        bool tw_changes;
        do { // Relatively cheap to compute, so repeat until convergence
            tw_changes = false;
            if (TW_R1(tw_changes) ||
                TW_R2(tw_changes) ||
                TW_R3(tw_changes) ||
                TW_R4((tw_changes))) {
                return true; // instance is infeasible
            }
            changes = tw_changes || changes;
        } while (tw_changes);

        changes = feasibilityPairs() || changes;
        changes = feasibilityTriples() || changes;

        changes = transitiveClosurePrecedence() || changes;
        changes = removeArcsFromPrecedencePairs() || changes;
    }
    return false; // instance may be feasible
}

/*
 * The first rule recognizes that the vehicle cannot depart from k earlier than it
 * can arrive there from another location
 */
bool TWPreprocessor::TW_R1(bool &changes) {
    for (uint k = 0; k < N; k++) {
        if (k == start)
            continue;
        static_array<uint, MAX_SIZE> earliest_arrival_times;
        for (uint i = 0; i < N; i++) {
            if (!E[i][k])
                continue;
            earliest_arrival_times.push_back(delta(i, k, e[i] + s[i]));
            if (earliest_arrival_times.back() <= e[k])
                break; // if any leq e[k] - stop here. It will not be possible to increase e[k]
        }

        if (earliest_arrival_times.empty()) {
            return true; // infeasible!
        }

        const uint new_val = clamp(*min_element(earliest_arrival_times.begin(), earliest_arrival_times.end()),
                                   e[k], l[k]);
        if (new_val > e[k]) {
            e[k] = new_val;
            changes = true;
            tightened_tws_e.add(k);
        }
    }
    return false;
}

/*
 * The second rule recognizes that the vehicle can only depart a location
 * at times that enable it to reach another location during its time window
 */
bool TWPreprocessor::TW_R2(bool &changes) {
    for (uint k = 0; k < N; k++) {
        if (k == end)
            continue;
        static_array<uint, MAX_SIZE> latest_departure_times;
        for (const uint &i: E[k]) {
            // find out the latest time we can leave k and still reach i before the end of its TW
            //   (while keeping some time to service it)
            const int latest_time = delta_inv(k, i, l[i] - s[i], e[k] + s[k], l[k]);
            if (latest_time >= 0) {
                latest_departure_times.push_back(static_cast<uint>(latest_time));
                if (latest_departure_times.back() >= l[k])
                    break; // if any geq l[k] - stop here. it will not be possible to decrease l[k]
            }
        }

        if (latest_departure_times.empty()) {
            return true; // infeasible!
        }

        const uint new_val = clamp(*max_element(latest_departure_times.begin(), latest_departure_times.end()),
                                   e[k], l[k]);
        if (new_val < l[k]) {
            l[k] = new_val;
            changes = true;
            tightened_tws_l.add(k);
        }
    }
    return false;
}

/*
 * Rule 3
 * recognizes that it is not necessary for the vehicle to depart from location k
 * to another location only to wait for that location’s time window to open
 */
bool TWPreprocessor::TW_R3(bool &changes) {
    for (uint k = 0; k < N; k++) {
        if (k == end)
            continue;

        static_array<uint, MAX_SIZE> latest_departure_times;
        for (const uint &i: E[k]) {
            // find out the latest time we can leave k and still reach i's opening
            // default value: leave as soon as possible when it is not possible to arrive early
            const int latest_time = delta_inv(k, i, e[i], e[k] + s[k], l[k]) - static_cast<int>(s[k]);
            latest_departure_times.push_back(latest_time >= 0 ? static_cast<uint>(latest_time) : e[k]);
            if (latest_departure_times.back() <= e[k])
                break; // if any lte e[k] - stop here. it will not be possible to increase e[k]
        }

        if (latest_departure_times.empty()) {
            return true; // infeasible!
        }

        const uint new_val = clamp(*min_element(latest_departure_times.begin(), latest_departure_times.end()),
                                   e[k], l[k]);
        if (new_val > e[k]) {
            e[k] = new_val;
            changes = true;
            tightened_tws_e.add(k);
        }
    }
    return false;
}

/*
 * Rule 4
 * recognizes that the vehicle need not depart from
 * location k later than the latest time at which it can arrive there
 */
bool TWPreprocessor::TW_R4(bool &changes) {
    for (uint k = 0; k < N; k++) {
        if (k == start)
            continue;
        static_array<uint, MAX_SIZE> latest_arrival_times;
        for (uint i = 0; i < N; i++) {
            if (!E[i][k])
                continue;
            latest_arrival_times.push_back(delta(i, k, l[i]) + s[k]);
            if (latest_arrival_times.back() >= l[k])
                break; // if any geq l[k] - stop here. it will not be possible to decrease l[k]
        }

        if (latest_arrival_times.empty()) {
            return true; // infeasible!
        }

        const uint new_val = clamp(*max_element(latest_arrival_times.begin(), latest_arrival_times.end()),
                                   e[k], l[k]);
        if (new_val < l[k]) {
            l[k] = new_val;
            changes = true;
            tightened_tws_l.add(k);
        }
    }
    return false;
}

/*
 * delta(i, j, t) returns the arrival time at j, leaving i at t (and waiting for the opening of j if necessary)
 * Note : it *does not* include service time at j.
 */
uint TWPreprocessor::delta(const uint &i, const uint &j, const uint &t) const {
    return max(e[j], t + this->tsp->cost(i, j, t));
}

/*
 * returns the latest departure time t' in [lb;ub] from i to j to reach j at time t
 * if t' < 0, infeasible
 */

int TWPreprocessor::delta_inv(const uint &i, const uint &j, const uint &t,
                              const uint &lb, const uint &ub) const {
    return this->tsp->cost_inv(i, j, t, lb, ub);
}

/*
 * delta_indirect(i, j, t) returns the arrival time at j, leaving i at t
 * using the shortest TW-aware path from i to j (i.e. the one arriving the earliest)
 * Note : as delta(i, j, t), *does not* include service time at j
 * (Dijkstra's algorithm)
 */
uint TWPreprocessor::delta_indirect(const uint &a, const uint &b, const uint &t) const {
    assert(!triangle_inequality_assumed); // should not call this fn if triangle inequality is assumed
    using pii = pair<uint, uint>;
    static priority_queue<pii, vector<pii>, greater<pii>> pq; // static to avoid allocating the vector at each call.
    clear(pq);

    static_array<uint, MAX_SIZE> sp(N, ~0u);
    sp[a] = t;
    pq.emplace(t, a);

    Bitset visited;
    while (!pq.empty()) {
        auto [t_i, i] = pq.top();
        pq.pop();
        if (i == b)
            return sp[b] - s[b]; // Arrival time, excluding service time
        if (t_i > sp[i])
            continue;
        visited.add(i);
        for (const uint &j: E[i] & ~visited) {
            const uint t_j = delta(i, j, sp[i]) + s[j]; // includes waiting for opening + service time
            if (t_j < sp[j] && t_j <= l[j]) {
                sp[j] = t_j;
                pq.emplace(t_j, j);
            }
        }
    }
    return l[b] + 1; // infeasible
}

/*
 * Adds the precedence constraint "a must be visited before b"
 * (Maintains both a precedence matrix and a precedence set)
 * Returns true if a change was made, false otherwise.
 */
bool TWPreprocessor::add_precedence(const uint &a, const uint &b) {
    if (R[a][b])
        return false;
    R[a].add(b);
    R_[b].add(a);
    removeArc(b, a);
    ++precedence_inferred;
    return true;
}

/*
 * If i->j is infeasible, j must precede i
 * Returns true if a change was made, false otherwise.
 */
bool TWPreprocessor::feasibilityPairs() {
    bool changes(false);
    for (uint i = 0; i < N; i++)
        for (uint j = 0; j < N; j++) {
            if (i == j || i == end || j == start || R[j][i])
                continue;
            if (infeasible_direct(i, j)) {
                // j cannot be a direct successor of i
                changes = removeArc(i, j) || changes;
                if (triangle_inequality_assumed || infeasible_indirect(i, j)) {
                    // j precedes i
                    add_precedence(j, i);
                    precedence_inferred_pairs++;
                    changes = true;
                }
            }
        }
    return changes;
}

/*
 * Tests feasibility between sequences of three nodes:
 *  for each pair (i, j) :
 *      if there exists a node k such that
 *        - ijk, kij and ikj are infeasible, j must precede i
 *        - ijk, kij are infeasible, arc (i, j) cannot be used in a valid solution
 * Returns true if a change was made, false otherwise.
 */
bool TWPreprocessor::feasibilityTriples() {
    bool changes(false);
    // First loops to compute feasibility for each (i,j,k) - the next loops will use three values each
    // (overlapping subproblems)
    static static_array<static_array<Bitset, MAX_SIZE>, MAX_SIZE> mem; // static to avoid stack overflow
    for (uint i = 0; i < N; i++) {
        mem[i].fill(N); // reset (as it is static)
        if (i == start || i == end)
            continue;
        for (uint j = 0; j < N; j++) {
            if (i == j || j == start || j == end)
                continue;
            for (uint k = 0; k < N; k++) {
                if (k == i || k == j || k == start || k == end)
                    continue;
                if (infeasible_mem(i, j, k))
                    mem[i][j].add(k);
            }
        }
    }

    for (uint i = 0; i < N; i++) {
        if (i == start || i == end)
            continue;
        for (uint j = 0; j < N; j++) {
            if (i == j || j == start || j == end)
                continue;
            for (uint k = 0; k < N && !R[j][i]; k++) { // stop if j precedes i (impossible to do better)
                if (k == i || k == j || k == start || k == end)
                    continue;
                const bool infeasible_ikj = mem[i][k][j]; // if !infeasible_ikj, the best we can do is remove arc (i,j)
                if (!infeasible_ikj && !E[i][j])
                    continue; // impossible to do better, arc (i,j) cannot be used already
                if (mem[i][j][k] && mem[k][i][j]) {
                    if (infeasible_ikj) {
                        if (add_precedence(j, i)) {
                            precedence_inferred_triples++;
                            changes = true;
                        }
                    } else
                        changes = removeArc(i, j) || changes;
                }
            }
        }
    }
    return changes;
}

/*
 * Determines if j can be a direct successor of i (in the best case scenario)
 */
bool TWPreprocessor::infeasible_direct(const uint &i, const uint &j) const {
    return delta(i, j, e[i] + s[i]) + s[j] > l[j];
}

/*
 * Determines if j can be an indirect successor of i (in the best case scenario)
 */
bool TWPreprocessor::infeasible_indirect(const uint &i, const uint &j) const {
    return delta_indirect(i, j, e[i] + s[i]) + s[j] > l[j];
}

/*
 * Determines if the sequence i -> j -> k (in the best case scenario) can be part of a feasible solution
 * (may need to compute shortest paths if triangle inequality is not verified)
 */
template<bool heuristic>
bool TWPreprocessor::infeasible(const uint &i, const uint &j, const uint &k) const {
    // Note : DO NOT check if arcs (i,j) and (j,k) are usable because we want to determine if
    //  i can be visited before j and j before k
    const uint end_i = e[i] + s[i];
    const uint end_j = s[j] + (heuristic ? delta(i, j, end_i) : delta_indirect(i, j, end_i));
    if (end_j > l[j])
        return true;
    const uint end_k = s[k] + (heuristic ? delta(j, k, end_j) : delta_indirect(j, k, end_j));
    return end_k > l[k];
}

bool TWPreprocessor::infeasible_mem(const uint &i, const uint &j, const uint &k) {
    // If it has once been proved infeasible, it will always be.
    if (infeasible_triples_mem[i * N + j].contains(k))
        return true;
    const bool ret = R[j][i] ||
                     R[k][i] ||
                     R[k][j] ||
                     infeasible(i, j, k);
    if (ret)
        infeasible_triples_mem[i * N + j].add(k);
    return ret;
}

bool TWPreprocessor::infeasible(const uint &i, const uint &j, const uint &k) const {
    const bool heuristic = infeasible<true>(i, j, k);
    if (triangle_inequality_assumed || !heuristic) {
        // a) if triangle inequality is verified, no need to compute shortest paths (i.e. "heuristic" is sufficient)
        // b) if it is feasible without computing the shortest paths, it is feasible
        return heuristic;
    }
    return infeasible<false>(i, j, k);
}

/*
 * Updates the precedence constraints so that if (i,k) and (k,j) belongs to it, (i,j) belongs to it as well.
 */
bool TWPreprocessor::transitiveClosurePrecedence() {
    bool changes(false);
    for (uint i = 0; i < N; i++)
        for (uint j = 0; j < N; j++) {
            if (i == j || R[i][j])
                continue;
            // if there exists at least a node k such that i precedes k and k precedes j
            if (!(R[i] & R_[j]).is_empty()) {
                add_precedence(i, j);
                changes = true;
                break;
            }
        }
    return changes;
}

/*
 * Remove arc (i, k) if there exists j such that (i, j) in P and (j, k) in P
 */
bool TWPreprocessor::removeArcsFromPrecedencePairs() {
    bool changes(false);
    for (uint i = 0; i < N; i++)
        for (const uint k: E[i]) { // for each usable arc (i, k) - k != i is implicit
            // if there exists at least a node j such that i precedes j and j precedes k
            if (!(R[i] & R_[k]).is_empty()) {
                removeArc(i, k); // always makes a change, since arc (i, k) is usable at the beginning
                changes = true;
            }
        }
    return changes;
}

const uint &TWPreprocessor::early(const uint &i) const {
    return e[i];
}

const uint &TWPreprocessor::late(const uint &i) const {
    return l[i];
}

void TWPreprocessor::late(const uint &i, const uint &v) {
    l[i] = v;
}

bool TWPreprocessor::removeArc(const uint &a, const uint &b) {
    if (E[a][b]) {
        E[a].remove(b);
        E_[b].remove(a);
        ++removed_arcs;
        return true;
    }
    return false;
}

void TWPreprocessor::reset_stats() {
    removed_arcs = 0;
    precedence_inferred = 0;
    precedence_inferred_pairs = 0;
    precedence_inferred_triples = 0;
    tightened_tws_e.clear();
    tightened_tws_l.clear();
}

void TWPreprocessor::print_stats() const {
    cout << removed_arcs << " RA, " << precedence_inferred << " PI (" << precedence_inferred_pairs << " PI2, "
         << precedence_inferred_triples << " PI3), (" << tightened_tws_e.count() << ", "
         << tightened_tws_l.count() << ") TWT, " << it << " IT, LB=" << e[end] << endl;
}

