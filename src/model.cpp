#include <vector>
#include <iostream>
#include <algorithm>
#include "cost_models.h"
#include "model.h"
#include "msa.hpp"
#include "tw_preprocessor.h"

using namespace std;

TDTSPTW::TDTSPTW(CostModel *cm, const uint &N_, const vector<TimeWindow> &tws, const vector<uint> &service_times,
                 const uint &start, const bool &triangle_inequality_verified, const bool &duplicate_depot,
                 const bool &propagate_precedence_only, const vector<pair<uint, uint>> &prec) :
        cm(cm),
        N(duplicate_depot ? N_ + 1 : N_),
        start(start),
        end(duplicate_depot ? N_ : N_ - 1),
        tws(tws),
        original_tws(tws),
        service_times(service_times),
        E(N, Bitset::full_set(N)),
        E_(N, Bitset::full_set(N)),
        R(N, Bitset()),
        R_(N, Bitset()),
        N_bs(Bitset::full_set(N)),
        N_0(Bitset(N_bs).remove(end)) {
    for (uint i = 0; i < N_; i++) {
        if (tws[i].early + service_times[i] > tws[i].late) {
            cout << "Infeasible instance (service times)" << endl;
            exit(1);
        }
    }
    if (duplicate_depot) {
        this->tws.push_back(tws[start]);
        this->original_tws.push_back(tws[start]);
        this->service_times.push_back(service_times[start]);
    }
    // can only leave start between [0; service_time[start]]
    this->tws[start].late = this->service_times[start];
    this->original_tws[start].late = this->service_times[start];

    p = new TWPreprocessor(this, triangle_inequality_verified, propagate_precedence_only, prec);
}

TDTSPTW::~TDTSPTW() {
    delete cm;
    delete p;
}

bool TDTSPTW::process_tws_and_compute_lbs(const uint &instance_ub) {
    // Precondition: there exist at least a solution of cost "instance_ub"
    // Returns true if the search is done (i.e. early[end] + service[end] == late[end])
    tws[end].late = instance_ub;
    process_tws_and_compute_lbs(true, true); // DO process TWS and DO compute LBs
    // Always returns true as there exist at least a solution of cost instance_ub
    return tws[end].early + service_times[end] == tws[end].late; // Cannot improve "instance_ub"
}

bool TDTSPTW::process_tws_and_compute_lbs(const bool &preprocess, const bool &compute_lbs) {
    // Returns true if the instance is proven to be infeasible
    if (preprocess) {
        p->late(end, tws[end].late);
        auto t_begin = std::chrono::steady_clock::now();
        bool infeasible = p->preprocess();

        // Get tighter time windows from preprocessor
        for (uint i = 0; i < N; i++) {
            this->tws[i].early = p->early(i);
            this->tws[i].late = p->late(i);
        }

        auto elapsed = chrono::duration_cast<chrono::duration<double >>(
                std::chrono::steady_clock::now() - (t_begin)).count();
        cout << "# Preprocessing time: " << elapsed << "s. - ";
        p->print_stats();
        p->reset_stats();
        if (infeasible) {
            cout << "# TW preprocessing proved infeasibility" << endl;
            return true;
        }
    }
    if (compute_lbs) {
        auto t_begin = std::chrono::steady_clock::now();

        this->cm->precompute_lbs(this->tws, this->service_times, this->E, this->E_);

        auto elapsed = chrono::duration_cast<chrono::duration<double >>(
                std::chrono::steady_clock::now() - (t_begin)).count();
        cout << "# Precompute LB time: " << elapsed << "s." << endl;
    }
    return false;
}

TDTSPTW::TDTSPTW(const PWConstantInstance &i, const uint &start) :
        TDTSPTW(i.force_fifo ? reinterpret_cast<CostModel *>(new PWConstantCostModel<true>(i.mat, i.N, i.n_ts,
                                                                                           i.w_ts, start))
                             : reinterpret_cast<CostModel *>(new PWConstantCostModel<false>(i.mat, i.N, i.n_ts,
                                                                                            i.w_ts, start)),
                i.N, i.tws, i.s, start) {};

TDTSPTW::TDTSPTW(const IGPInstance &inst,
                 const uint &start) :
        TDTSPTW(get_igp_cost_model(inst, start),
                inst.N, inst.tws, vector<uint>(inst.N, 0),
                start) {};

TDTSPTW::TDTSPTW(const ConstantInstance &inst,
                 const uint &start, const bool triangle_inequality_verified) :
        TDTSPTW(new ConstantCostModel(inst.mat, inst.N, start),
                inst.N, inst.tws, vector<uint>(inst.N, 0), start, triangle_inequality_verified) {};

void TDTSPTW::dump() {
    cout << "Time windows, service times  : " << endl;
    for (uint i = 0; i < N; i++)
        cout << i << " : [" << tws[i].early << ";" << tws[i].late << "], " << service_times[i] << " "
             << (i == start ? "(start)" : "")
             << (i == end ? "(end)" : "") << endl;
    cout << endl;

    cout << "R matrix : " << endl;
    for (uint i = 0; i < N; i++) {
        for (uint j = 0; j < N; j++) {
            cout << (R[i].contains(j) ? '1' : '_') << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "E matrix : " << endl;
    for (uint i = 0; i < N; i++) {
        for (uint j = 0; j < N; j++) {
            cout << (E[i].contains(j) ? '1' : '_') << " ";
        }
        cout << endl;
    }
    cout << endl;
}

uint TDTSPTW::cost(const uint &i, const uint &j, const uint &t) const {
    return cm->cost(i, j, t);
}

int TDTSPTW::cost_inv(const uint &i, const uint &j, const uint &t, const uint &lb, const uint &ub) const {
    return cm->cost_inv(i, j, t, lb, ub);
}

State TDTSPTW::initial_state() const {
    return State(start, Bitset(), 0);
}

Bitset TDTSPTW::A_superset(const State &s) const {
    // Superset of possible actions (does not consider arrival time / precedence constraints)
    Bitset leftToVisit = N_0 & ~s.get_visited();
    if (leftToVisit.is_empty())
        leftToVisit = Bitset::unit_set(end);
    return leftToVisit & E[s.i];
}

void TDTSPTW::A(const State &s, static_array<State, MAX_SIZE> &children) const {
    const Bitset not_visited = ~s.get_visited();
    for (const auto &a: A_superset(s))
        if ((R_[a] & not_visited).is_empty()) {
            // all predecessors of a have been visited
            const State &s_ = tau(s, a);
            if (s_.t <= tws[s_.i].late)
                children.push_back(s_);
        }
}

bool TDTSPTW::is_terminal(const State &s) const {
    return s.i == this->end;
}

uint TDTSPTW::upper_bound_from_tw() const {
    return tws[end].late + service_times[end] + 1;
}

Action TDTSPTW::greedy_action(State s) const {
    // Precondition : A(s) is not empty (i.e. at least one feasible action
    auto sort = [this, s](const uint &a, const uint &b) -> bool {
        // sort by (e_a, l_a, c(i, a, t)) - Lexicographic order
        return this->tws[a].late < this->tws[b].late ||

               (this->tws[a].late == this->tws[b].late &&
                this->tws[a].early < this->tws[b].early) ||

               (this->tws[a].late == this->tws[b].late &&
                this->tws[a].early == this->tws[b].early &&
                this->cost(s.i, a, s.t) < this->cost(s.i, b, s.t));
    };
    const Bitset e = A(s);
    return *std::min_element(e.begin(), e.end(), sort);
}

uint TDTSPTW::greedy_upper_bound(State s) const {
    // Tries to compute greedily an upper bound for completing the partial solution s.
    // Returns ~0u if no feasible completion was found.
    if (is_terminal(s))
        return s.t;
    if (A(s).is_empty())
        return ~0u;
    return greedy_upper_bound(tau(s, greedy_action(s)));
}

uint TDTSPTW::h_0(const State &s, const uint &thread_number) const {
    return 0;
}
