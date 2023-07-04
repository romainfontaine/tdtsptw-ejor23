#ifndef TDTSPTW_TW_PREPROCESSOR_H
#define TDTSPTW_TW_PREPROCESSOR_H

#include <vector>
#include <iostream>
#include "bitsets.hpp"
#include "util.hpp"
#include "constants.h"
#include "static_array.hpp"

using namespace std;

class TDTSPTW;

class TWPreprocessor {
    TDTSPTW *tsp;
    uint N;
    vector<Bitset> &E;
    vector<Bitset> &E_;
    vector<Bitset> &R;
    vector<Bitset> &R_;
    vector<uint> e;
    vector<uint> l;
    vector<uint> s;
    const uint start;
    const uint end;
    const bool triangle_inequality_assumed;
    static_array<Bitset, MAX_SIZE * MAX_SIZE> infeasible_triples_mem;

    bool propagate_precedence_only;

    uint it{};
    uint removed_arcs{};
    uint precedence_inferred{};
    uint precedence_inferred_pairs{};
    uint precedence_inferred_triples{};
    Bitset tightened_tws_l{};
    Bitset tightened_tws_e{};

public:
    TWPreprocessor(TDTSPTW *tsp,
                   const bool &triangle_inequality_assumed,
                   const bool &propagate_precedence_only = false,
                   const vector<pair<uint, uint>> &precedence_constraints = vector<pair<uint, uint>>());

    bool preprocess();

    bool add_precedence(const uint &a, const uint &b);

    const uint &early(const uint &i) const;

    const uint &late(const uint &i) const;

    void late(const uint &i, const uint &v);

    void reset_stats();

    void print_stats() const;

private:
    bool TW_R1(bool &changes);

    bool TW_R2(bool &changes);

    bool TW_R3(bool &changes);

    bool TW_R4(bool &changes);

    uint delta(const uint &i, const uint &j, const uint &t) const;

    int delta_inv(const uint &i, const uint &j, const uint &t, const uint &lb, const uint &ub) const;

    uint delta_indirect(const uint &a, const uint &b, const uint &t) const;

    bool feasibilityPairs();

    bool feasibilityTriples();

    bool infeasible_direct(const uint &i, const uint &j) const;

    bool infeasible_indirect(const uint &i, const uint &j) const;

    template<bool heuristic>
    bool infeasible(const uint &i, const uint &j, const uint &k) const;

    bool infeasible(const uint &i, const uint &j, const uint &k) const;

    bool infeasible_mem(const uint &i, const uint &j, const uint &k);

    bool transitiveClosurePrecedence();

    bool removeArcsFromPrecedencePairs();

    bool removeArc(const uint &a, const uint &b);
};

#endif //TDTSPTW_TW_PREPROCESSOR_H
