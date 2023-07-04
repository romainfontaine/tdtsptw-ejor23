#ifndef TDTSPTW_MODEL_H
#define TDTSPTW_MODEL_H

#include <string>
#include <vector>
#include <sstream>
#include "cost_models.h"
#include "constants.h"
#include "bitsets.hpp"
#include "omp.h"
#include "msa.hpp"
#include "util.hpp"

using namespace std;

struct State {
    uint i;
    uint t;
    Bitset S;

    explicit State(const uint &i = 0, const Bitset &S = Bitset(), const uint &t = 0) {
        this->i = i;
        this->S = S;
        this->t = t;
    }

    friend std::ostream &operator<<(std::ostream &stream, const State &s) {
        stream << "State(" << s.i << "," << s.S << "," << s.t << ")";
        return stream;
    }

    bool operator==(const State &other) const {
        return (i == other.i && S == other.S);
    }

    Bitset get_visited() const {
        return Bitset::unit_set(i) | S;
    }
};

namespace std {
    template<>
    struct hash<State> {
        size_t operator()(const State &k) const {
            // t is excluded for dominance
            return hash<uint>()(k.i) + hash<Bitset>()(k.S);
        }
    };
}

class CostModel;

class TDTSPTW;

using SeqHeuristic = uint (TDTSPTW::*)(const State &, const uint &) const;

using Heuristic = void (TDTSPTW::*)(const State &s,
                                    static_array<State, MAX_SIZE> &children,
                                    static_array<uint, MAX_SIZE> &h_values) const;

class TWPreprocessor;

class TDTSPTW {
    CostModel *cm;
public:
    const uint N, start, end;
    vector<TimeWindow> tws;
    vector<TimeWindow> original_tws;
    vector<uint> service_times;
    vector<Bitset> E; // Usable arcs : arc (i,j) may be used if E[i][j]
    vector<Bitset> E_; // Transposed
    vector<Bitset> R; // Precedence constraints : node i must precede j iff R[i][j]
    vector<Bitset> R_; // Transposed

    const Bitset N_bs;
    const Bitset N_0;
    TWPreprocessor *p;

    TDTSPTW(CostModel *cm, const uint &N_, const vector<TimeWindow> &tws, const vector<uint> &service_times,
            const uint &start, const bool &triangle_inequality_verified = false, const bool &duplicate_depot = true,
            const bool &propagate_precedence_only = false,
            const vector<pair<uint, uint>> &prec = vector<pair<uint, uint>>());

    ~TDTSPTW();

    bool process_tws_and_compute_lbs(const uint &instance_ub);

    bool process_tws_and_compute_lbs(const bool &preprocess, const bool &compute_lbs);

    TDTSPTW(const PWConstantInstance &i, const uint &start);

    TDTSPTW(const IGPInstance &inst,
            const uint &start);

    TDTSPTW(const ConstantInstance &inst,
            const uint &start, const bool triangle_inequality_verified);

    TDTSPTW(const TDTSPTW &) = delete;

    static CostModel * get_igp_cost_model(const IGPInstance &inst, const uint &start){
        if(inst.use_correct_rounding){
            if(inst.memorize_costs)
                return new IGPCostModel<true,true>(inst, start);
            else
                return new IGPCostModel<true,false>(inst, start);
        }else{
            if(inst.memorize_costs)
                return new IGPCostModel<false,true>(inst, start);
            else
                return new IGPCostModel<false,false>(inst, start);
        }
    }

    void dump();

    uint cost(const uint &i, const uint &j, const uint &t) const;

    int cost_inv(const uint &i, const uint &j, const uint &t, const uint &lb, const uint &ub) const;

    State initial_state() const;

    template<bool use_original_tws = false>
    uint arrival_time(const State &s, const Action &a) const {
        // Time after applying action A to state S (includes optional waiting at A and its service time)
        const uint arrival_time = s.t + cost(s.i, a, s.t);
        const uint begin_service = max(arrival_time, use_original_tws ? original_tws[a].early : tws[a].early);
        return begin_service + service_times[a];
    }

    template<bool use_original_tws = false>
    State tau(const State &s, const Action &a) const {
        return State(a, s.S | Bitset::unit_set(s.i), arrival_time<use_original_tws>(s, a));
    }

    Bitset A_superset(const State &s) const;

    void A(const State &s, static_array<State, MAX_SIZE> &children) const;

    template<bool use_original_tws = false>
    Bitset A(const State &s) const {
        const Bitset not_visited = ~s.get_visited();
        Bitset ret;
        for (const uint &dest: A_superset(s))
            if (arrival_time<use_original_tws>(s, dest) <= tws[dest].late &&
                (use_original_tws || (R_[dest] & not_visited).is_empty()))
                // moving from s.i to dest + servicing dest can be done before dest's deadline
                //   AND all of dest prerequisites have been visited (if !use_original_tws)
                ret.add(dest); // reachable
        return ret;
    }

    bool is_terminal(const State &s) const;

    uint upper_bound_from_tw() const;

    Action greedy_action(State s) const;

    uint greedy_upper_bound(State s) const;

    uint get_ldt(const uint &a, const uint &b) const {
        return cm->get_ldt(a, b);
    }

    uint h_0(const State &s, const uint &thread_number) const;

    template<bool h_out = false,
            LBCostType lb_type,
            bool use_exact_i_cost = false,
            bool use_precedence_constraints_i = true>
    uint h_0_out_feasibility_checks(const State &s, const uint &thread_number = 0) const {
        if ((N_0 & ~s.get_visited()).is_empty())
            // Two cases: a) s is a final state (cost end -> end will be 0), or
            //            b) all nodes have been visited except the end depot.
            return cost(s.i, end, s.t);

        CostLB lb = CostModel::lbs[static_cast<int>(lb_type)];

        const Bitset X = N_0 & ~s.get_visited();     // X - unvisited, without end depot
        const Bitset Xn = X | Bitset::unit_set(end); // X U {end} - unvisited, with end depot

        Bitset has_incoming_arc;
        uint min_time = ~0u;

        // Arcs leaving s.i & compute min. time of arrival on another node.
        for (const uint &y: X) { // Note that end is excluded (at least 1 node to visit before)
            if ((s.t > cm->get_ldt(s.i, y)) ||
                (use_precedence_constraints_i && (R_[y] & Xn).any()))
                continue;
            const uint c(use_exact_i_cost ?
                         (       // include waiting time if the LB does (i.e. != naive).
                                 lb_type == LBCostType::naive ?
                                 cost(s.i, y, s.t) :
                                 max(s.t + cost(s.i, y, s.t), tws[y].early) - s.t
                         ) :
                         (cm->*lb)(s.i, y, s.t)); // Use LB otherwise.
            min_time = min(min_time, s.t + c + service_times[y]);
            has_incoming_arc.add(y);
        }
        if (min_time == ~0u) // Infeasible - no arc leaving s.i
            return tws[end].late + 1;

        const auto &E_t = cm->get_E_t<false>(min_time);
        for (const uint &x: X)
            if ((Xn & E_t[x]).is_empty())
                return tws[end].late + 1; // At least one node has no *outgoing* arc.

        const auto &E_t_ = cm->get_E_t<true>(min_time);
        for (const uint &x: Xn)
            if (!has_incoming_arc[x] && (X & E_t_[x]).is_empty())
                return tws[end].late + 1; // At least one node has no *incoming* arc.

        if (h_out)
            return min_time - s.t;
        else
            return 0;
    }
	template<bool use_exact_costs=false, uint fea_first=0> // fea_first = 0 -> no fea | 1 -> fea | 2 -> fea+
	uint h_fea_spp(const State &s, const uint &thread_number = 0) const {
        static_assert(0<=fea_first && fea_first<=2);
        if(fea_first != 0){
            constexpr bool exact_costs = fea_first == 2;
            const auto fea = h_0_out_feasibility_checks<false, LBCostType::medium, exact_costs>(s, thread_number);
            if (fea != 0 || s.i == end)
                return fea;
        } else if ((N_0 & ~s.get_visited()).is_empty())
			// Two cases: a) s is a final state (cost end -> end will be 0), or
			//            b) all nodes have been visited except the end depot.
			return cost(s.i, end, s.t);

        const uint INF = tws[end].late + 1;
		CostLB lb = &CostModel::lb_medium;

        const Bitset X = N_0 & ~s.get_visited();     // X - unvisited, without end depot
        const Bitset Xn = X | Bitset::unit_set(end); // X U {end} - unvisited, with end depot

		static_array<uint, MAX_SIZE> sp(N, ~0u); // Useless to set "sp[s.i] = s.t;"
        Bitset visited; // Useless to add s.i to this set

		using pii = pair<uint, uint>;
		const auto comp = []( const pii &a, const pii &b) { return a.first > b.first; };
        using pq_t = priority_queue<pii, vector<pii>, decltype(comp)>;
		static vector<pq_t> pqs(MAX_THREAD_NUMBER, pq_t(comp));
        pq_t &pq(pqs[thread_number]);
		clear(pq);
        // Initialize Dijkstra with arcs leaving s.i
        const auto &E_ti = cm->get_E_t<false>(s.t)[s.i];
        for (const uint &j: E_ti & X) { // s.i to any unvisited node, except end.
            if ((R_[j] & Xn).any()) // j is preceded by an unvisited node
                continue;
            const auto &t_j = max(tws[j].early, s.t + cost(s.i, j, s.t));
            sp[j] = t_j;
            pq.emplace(t_j, j);
        }
        // Main Dijkstra loop - single source shortest paths using TWs + detecting nodes without outgoing edges.
        while (!pq.empty()) {
            const auto [t_i, i] = pq.top();
            pq.pop();
            if (visited[i]) // equivalent to "if(t_i > sp[i])"
                continue;
            visited.add(i);
            if(i == end)
                continue;
            const auto &e = cm->get_E_t(t_i)[i] & Xn; // s.i excluded
            if(e.is_empty())
                return INF; // i has no outgoing arc
            for (const uint &j: e & ~visited) {
                const uint t_j = max(tws[j].early, t_i + (use_exact_costs ?
                                                          cost(i, j, t_i) :
                                                          (cm->*lb)(i, j, t_i)));
                if (t_j < sp[j]) {
                    sp[j] = t_j;
                    pq.emplace(t_j, j);
                }
            }
        }
        if (visited != Xn)
            return INF; // At least one node in Xn has no incoming arc
		return 0;
	}

    template<LBCostType lb_type,
            typename Callback,
            typename CallbackEndRow,
            LDTFilteringMode ldt_mode,
            bool exact_i_cost,
            bool use_precedence_constraints_i,
            bool dfs_reachability_check>
    bool compute_subgraph(const State &s,
                          Callback callback, CallbackEndRow callback_end_row, uint &sum_services) const {
        // returns false if infeasibility was proven.
        CostLB lb = CostModel::lbs[static_cast<int>(lb_type)];

        const Bitset X = N_0 & ~s.get_visited();     // X - unvisited, without end depot
        const Bitset Xn = X | Bitset::unit_set(end); // X U {end} - unvisited, with end depot

        static_array<Bitset, MAX_SIZE> adj(N); // Only needed for DFS - should be optimized out.

        Bitset has_incoming_arc;
        uint min_time = ~0u;

        // Arcs leaving s.i & compute min. time of arrival on another node.
        for (const uint &y: X & E[s.i]) { // end excluded as at least one node must be visited before
            if ((ldt_mode != LDTFilteringMode::NONE && s.t > cm->get_ldt(s.i, y)) ||
                (use_precedence_constraints_i && (R_[y] & Xn).any()))
                continue;
            const uint c(exact_i_cost ?
                         (       // include waiting time if the LB does (i.e. != naive).
                                 lb_type == LBCostType::naive ?
                                 cost(s.i, y, s.t) :
                                 max(s.t + cost(s.i, y, s.t), tws[y].early) - s.t
                         ) :
                         (cm->*lb)(s.i, y, s.t)); // Use LB otherwise.
            min_time = min(min_time, s.t + c + service_times[y]);
            callback(s.i, y, c);
            has_incoming_arc.add(y);
            if (dfs_reachability_check)
                adj[s.i].add(y);
        }
        if (min_time == ~0u) // Infeasible - no arc leaving s.i
            return false;
        callback_end_row();

        const auto &E_t = cm->get_E_t(min_time);
        // Other arcs (necessarily used at or after min_time)
        for (const uint &x: X) {
            sum_services += service_times[x];
            bool no_outgoing_arc = true;
            for (const uint &y: Xn & (ldt_mode == LDTFilteringMode::PRECOMPUTED ?
                                      E_t[x] : E[x])) {
                if (ldt_mode == LDTFilteringMode::CONDITION && min_time > cm->get_ldt(x, y))
                    continue;
                callback(x, y, (cm->*lb)(x, y, min_time));
                has_incoming_arc.add(y);
                no_outgoing_arc = false;
                if (dfs_reachability_check)
                    adj[x].add(y);
            }
            if (no_outgoing_arc)
                return false;
            callback_end_row();
        }
        // Make sure all nodes in Xn have at least one incoming arc.
        if ((Xn & ~has_incoming_arc).any())
            return false;

        // Make sure all nodes are reachable from s.i
        if (dfs_reachability_check) {
            Bitset reachable = adj[s.i];
            Bitset reachable_unvisited(reachable);
            Bitset visited;
            while ((Xn & ~reachable).any() && !reachable_unvisited.is_empty()) {
                visited |= reachable_unvisited;
                for (const auto &x: reachable_unvisited)
                    reachable |= adj[x];
                reachable_unvisited = reachable & ~visited;
            }
            if ((Xn & ~reachable).any())
                return false;
        }
        return true;
    }

    template<LBCostType lb_type,
            typename Callback,
            typename CallbackEndCol,
            LDTFilteringMode ldt_mode,
            bool exact_i_cost,
            bool use_precedence_constraints_i,
            bool dfs_reachability_check>
    bool compute_subgraph_col_major(const State &s,
                                    Callback callback, CallbackEndCol callback_end_col,
                                    uint &sum_services) const {
        // returns false if infeasibility was proven.
        CostLB lb = CostModel::lbs[static_cast<int>(lb_type)];

        const Bitset X = N_0 & ~s.get_visited();     // X - unvisited, without end depot
        const Bitset Xn = X | Bitset::unit_set(end); // X U {end} - unvisited, with end depot

        static_array<Bitset, MAX_SIZE> adj(N); // Only needed for DFS - should be optimized out.

        static_array<uint, MAX_SIZE> leaving_si(N);
        Bitset leaving_si_bs;

        Bitset has_outgoing_arc;
        uint min_time = ~0u;
        // Arcs leaving s.i & compute min. time of arrival on another node.
        for (const uint &y: X & E[s.i]) { // Note that end is excluded (at least 1 node to visit before)
            if ((ldt_mode != LDTFilteringMode::NONE && s.t > cm->get_ldt(s.i, y)) ||
                (use_precedence_constraints_i && (R_[y] & Xn).any()))
                continue;
            const uint c(exact_i_cost ?
                         (       // include waiting time if the LB does (i.e. != naive).
                                 lb_type == LBCostType::naive ?
                                 cost(s.i, y, s.t) :
                                 max(s.t + cost(s.i, y, s.t), tws[y].early) - s.t
                         ) :
                         (cm->*lb)(s.i, y, s.t)); // Use LB otherwise.
            min_time = min(min_time, s.t + c + service_times[y]);
            leaving_si[y] = c;
            leaving_si_bs.add(y);
            if (dfs_reachability_check)
                adj[s.i].add(y);
        }
        if (min_time == ~0u) // Infeasible - no arc leaving s.i
            return false;

        const auto &E_t_ = cm->get_E_t<true>(min_time);

        for (const uint &y: Xn) {
            sum_services += service_times[y];
            bool y_has_incoming = leaving_si_bs[y];
            if (leaving_si_bs[y])
                callback(s.i, y, leaving_si[y]);

            for (const uint &x: X & (ldt_mode == LDTFilteringMode::PRECOMPUTED ?
                                     E_t_[y] : E_[y])) { // s.i excluded
                if (ldt_mode == LDTFilteringMode::CONDITION && min_time > cm->get_ldt(x, y))
                    continue;
                callback(x, y, (cm->*lb)(x, y, min_time));
                y_has_incoming = true;
                has_outgoing_arc.add(x);
                if (dfs_reachability_check)
                    adj[x].add(y);
            }
            if (!y_has_incoming)
                return false;
            callback_end_col();
        }
        if ((X & ~has_outgoing_arc).any())
            return false;

        // Make sure all nodes are reachable from s.i
        if (dfs_reachability_check) {
            Bitset reachable = adj[s.i];
            Bitset reachable_unvisited(reachable);
            Bitset visited;
            while ((Xn & ~reachable).any() && !reachable_unvisited.is_empty()) {
                visited |= reachable_unvisited;
                for (const auto &x: reachable_unvisited)
                    reachable |= adj[x];
                reachable_unvisited = reachable & ~visited;
            }
            if ((Xn & ~reachable).any())
                return false;
        }
        return true;
    }

    template<LBCostType lb_type,
            LDTFilteringMode ldt_mode = LDTFilteringMode::PRECOMPUTED,
            bool use_exact_i_cost = false,
            bool use_precedence_constraints_i = true,
            bool dfs_reachability_check = false>
    uint h_outgoing_arcs(const State &s, const uint &thread_number = 0) const {
        if ((N_0 & ~s.get_visited()).is_empty())
            // Two cases: a) s is a final state (cost end -> end will be 0), or
            //            b) all nodes have been visited except the end depot.
            return cost(s.i, end, s.t);

        uint sum = 0;
        uint min_c = ~0u;
        const auto &callback_arc = [&min_c](const uint &x, const uint &y, const uint &c) {
            min_c = min(min_c, c);
        };
        const auto &callback_end_row = [&sum, &min_c]() {
            sum += min_c;
            min_c = ~0u;
        };
        uint sum_services(0);
        if (!compute_subgraph<lb_type, decltype(callback_arc), decltype(callback_end_row),
                ldt_mode, use_exact_i_cost, use_precedence_constraints_i,
                dfs_reachability_check>(s, callback_arc, callback_end_row, sum_services))
            return tws[end].late + 1;
        return sum + sum_services;
    }

    template<LBCostType lb_type,
            LDTFilteringMode ldt_mode = LDTFilteringMode::PRECOMPUTED,
            bool use_exact_i_cost = false,
            bool use_precedence_constraints_i = true,
            bool dfs_reachability_check = false>
    uint h_outgoing_incoming_arcs(const State &s, const uint &thread_number = 0) const {
        if ((N_0 & ~s.get_visited()).is_empty())
            // Two cases: a) s is a final state (cost end -> end will be 0), or
            //            b) all nodes have been visited except the end depot.
            return cost(s.i, end, s.t);

        uint sum_oa = 0;
        uint min_oa = ~0u;
        static_array<uint, MAX_SIZE> min_ia(N, ~0u);

        const auto &callback_arc = [&min_oa, &min_ia](const uint &x, const uint &y, const uint &c) {
            min_oa = min(min_oa, c);
            min_ia[y] = min(min_ia[y], c);
        };
        const auto &callback_end_row = [&min_oa, &sum_oa]() {
            sum_oa += min_oa;
            min_oa = ~0u;
        };
        uint sum_services(0);
        if (!compute_subgraph<lb_type, decltype(callback_arc), decltype(callback_end_row),
                ldt_mode, use_exact_i_cost, use_precedence_constraints_i,
                dfs_reachability_check>(s, callback_arc, callback_end_row, sum_services))
            return tws[end].late + 1;

        const Bitset &Xn(N_bs &
                         ~s.get_visited());
        uint sum_ia(0);
        for (const auto &x: Xn)
            sum_ia += min_ia[x];

        return max(sum_oa, sum_ia) + sum_services;
    }

    template<LBCostType lb_type,
            LDTFilteringMode ldt_mode = LDTFilteringMode::PRECOMPUTED,
            bool use_exact_i_cost = false,
            bool use_precedence_constraints_i = true,
            bool dfs_reachability_check = false>
    uint h_msa(const State &s, const uint &thread_number = 0) const {
        if ((N_0 & ~s.get_visited()).is_empty())
            return cost(s.i, end, s.t);

        const auto Xn = N_bs & ~s.get_visited();

        auto &msa = *MSA<MAX_SIZE>::get_instance(thread_number);
        msa.init(N);

        static_array<pair<uint, uint>, MAX_SIZE> row;
        uint cur_y;

        const auto &callback_arc = [&row, &cur_y](const uint &x, const uint &y, const uint &c) {
            cur_y = y;
            row.emplace_back(x, c);
        };

        const auto &callback_end_col = [&msa, &cur_y, &row]() {
            using uu = pair<uint, uint>;

            make_heap(row.begin(), row.end(), [](const uu &x, const uu &y) {
                return x.second > y.second;
            }); // insert in heap order.
            msa.add_edges(cur_y, row);
            row.clear();
        };

        uint sum_services(0);
        if (!compute_subgraph_col_major<lb_type,
                decltype(callback_arc), decltype(callback_end_col), ldt_mode, use_exact_i_cost,
                use_precedence_constraints_i, dfs_reachability_check>(s, callback_arc, callback_end_col, sum_services))
            return tws[end].late + 1;

        return msa.solve(s.i, tws[end].late + 1, Xn) + sum_services;
    }

    template<SeqHeuristic h>
    void h_par(const State &s, static_array<State, MAX_SIZE> &children,
               static_array<uint, MAX_SIZE> &h_values) const {
#ifndef _OPENMP
        cerr<<"Cannot use parallel heuristics when compiled without OpenMP, exiting."<<endl;
    exit(1);
#endif
        const uint size = children.size();
        h_values.resize(size);
#pragma omp parallel for shared(children, h_values, s, size), schedule(dynamic, 1), default(none)
        for (uint i = 0; i < size; i++) {
            h_values[i] = (this->*h)(children[i], static_cast<uint>(omp_get_thread_num()));
        }
    }

    template<SeqHeuristic H>
    void h_seq(const State &s, static_array<State, MAX_SIZE> &children,
               static_array<uint, MAX_SIZE> &h_values) const {
        for (const auto &s_: children)
            h_values.push_back((this->*H)(s_, 0));
    }
};

#endif //TDTSPTW_MODEL_H
