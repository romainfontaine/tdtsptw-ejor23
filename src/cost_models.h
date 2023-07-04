#ifndef TDTSPTW_COST_MODELS_H
#define TDTSPTW_COST_MODELS_H

#include <vector>
#include <cassert>
#include <algorithm>
#include "constants.h"
#include "instances.h"
#include "bitsets.hpp"
#include "static_array.hpp"
#include "util.hpp"

class CostModel;

using CostLB = uint(CostModel::*)(const uint &, const uint &, const uint &) const;

typedef unsigned int uint;

using namespace std;

class CostModel {
protected:
    uint N;
    vector<vector<uint>> lb_naive_;
    vector<static_array<uint, MAX_SIZE>> lb_medium_;
    vector<vector<uint>> ldt;
    vector<static_array<Bitset, MAX_SIZE>> E_t; // Sets of usable arcs at or after time t
    vector<static_array<Bitset, MAX_SIZE>> E_t_; // Transposed
    vector<uint> E_t_time_map;

public:
    CostModel(const uint &N) : N(N),
                               lb_naive_(N, vector<uint>(N)),
                               lb_medium_(N, static_array<uint, MAX_SIZE>(N)),
                               ldt(N, vector<uint>(N)) {
    };

    virtual ~CostModel() {};

    virtual uint cost(const uint &i, const uint &j, const uint &t) const = 0;

    uint lb_naive(const uint &i, const uint &j, const uint &t) const {
        return lb_naive_[i][j];
    };

    uint lb_medium(const uint &i, const uint &j, const uint &t) const {
        return lb_medium_[i][j];
    };

    uint get_ldt(const uint &i, const uint &j) const {
        // Latest time to leave i and still reach j on time
        // precondition : E[i].contains(j)
        return ldt[i][j];
    }

    template<bool transposed = false>
    static_array<Bitset, MAX_SIZE> const &get_E_t(const uint &t) const {
        if (t >= E_t_time_map.size())
            return (transposed ? E_t_ : E_t).back();
        return (transposed ? E_t_ : E_t)[E_t_time_map[t]];
    }

    void precompute_LDTs(const vector<TimeWindow> &tws,
                         const vector<uint> &service_times,
                         const vector<Bitset> &E,
                         const vector<Bitset> &E_) {
        for (uint i = 0; i < N - 1; i++)
            for (const auto &j: E[i])
                ldt[i][j] = static_cast<uint>(cost_inv(i, j, tws[j].late - service_times[j],
                                                       tws[i].early + service_times[i], tws[i].late));

        precompute_E_t(E, E_);
    }

    void precompute_E_t(const vector<Bitset> &E,
                        const vector<Bitset> &E_) {
        uint n_usable_arcs(0);
        for (uint i = 1; i < N - 1; i++) // Arcs leaving 0 excluded as they can only be used at t=0
            n_usable_arcs += E[i].count();

        vector<tuple<uint, uint, uint>> sorted_LDTs; // will store (ldt, i, j) for each usable arc (i, j) in E.
        sorted_LDTs.reserve(n_usable_arcs);

        for (uint i = 1; i < N - 1; i++)
            for (const auto &j: E[i])
                sorted_LDTs.emplace_back(ldt[i][j], i, j);

        sort(sorted_LDTs.begin(), sorted_LDTs.end()); // Non decreasing order (by (ldt, i, j), lexico.)

        E_t.clear();
        // E_t[t] gives a list of bitset, for each (interesting) time.
        // Each list of bitsets represents the arcs that can be used at time t or after (0(n*n/64) memory),
        // There are at most O(#E) values of t,
        // => O(#E*n*n/64) => O(n^4) - with a small constant factor.

        E_t_.clear(); // E_t transposed.

        vector<uint> times;

        uint current_time(0);
        static_array<Bitset, MAX_SIZE> next_E(E); /* Contains arcs usable at or after *current_time*
                                                                 (All usable arcs, initially) */
        static_array<Bitset, MAX_SIZE> next_E_(E_); // Transposed

        for (const auto &[ldt_xy, x, y]: sorted_LDTs) {
            if (current_time != ldt_xy) {
                E_t.push_back(next_E);
                E_t_.push_back(next_E_);

                times.push_back(current_time);
                current_time = ldt_xy;
            }
            next_E[x].remove(y);
            next_E_[y].remove(x);
        }
        E_t.push_back(next_E);
        E_t_.push_back(next_E_);

        times.push_back(current_time);

        E_t_time_map.clear();
        E_t_time_map.push_back(0);
        for (uint i = 0; i < times.size() - 1; i++)
            for (uint j = times[i]; j < times[i + 1]; j++)
                E_t_time_map.push_back(i);
        E_t_time_map.push_back(times.size() - 1);
    }

    virtual void precompute_lbs(const vector<TimeWindow> &tws,
                                const vector<uint> &service_times,
                                const vector<Bitset> &E,
                                const vector<Bitset> &E_) {
        precompute_LDTs(tws, service_times, E, E_);
        // Naive : min cost for each t in [ei+si;li+si]
        for (uint i = 0; i < N - 1; i++) {
            const uint low = tws[i].early + service_times[i];
            const uint high = tws[i].late;
            for (const auto &j: E[i]) {
                lb_naive_[i][j] = cost(i, j, low);
                for (uint t = low + 1; t <= high; t++)
                    lb_naive_[i][j] = min(lb_naive_[i][j], cost(i, j, t));
            }
        }
        // Same as naive but takes into account 1) waiting times (tw's opening) 2) j's TW
        // should dominate the naive one as it is tighter & it is not slower to access
        for (uint i = 0; i < N - 1; i++) {
            const uint low = tws[i].early + service_times[i];
            for (const auto &j: E[i]) {
                const uint high = ldt[i][j];

                lb_medium_[i][j] = max(tws[j].early, low + cost(i, j, low)) - low;
                for (uint t = low + 1; t <= high; t++) {
                    lb_medium_[i][j] = min(lb_medium_[i][j],
                                           max(tws[j].early, t + cost(i, j, t)) - t);
                }
                assert(lb_medium_[i][j] >= lb_naive_[i][j]);
            }
        }
    }

    template<LBCostType lb_type>
    void dump_cost_bounds() {
        for (uint i = 0; i < N; i++) {
            for (uint j = 0; j < N; j++)
                cout << (this->*lbs[static_cast<uint>(lb_type)])(i, j, 0) << " ";
            cout << endl;
        }
    }

    virtual int cost_inv(const uint &i, const uint &j, const uint &t,
                         const uint &lb, const uint &ub) const = 0;

    int cost_inv_linear(const uint &i, const uint &j, const uint &t,
                        const uint &lb, const uint &ub) const {
        // To be used for non-fifo cost functions
        assert(lb <= ub);
        // for k in reverse([lb;ub]) - both bounds included
        for (uint k = ub; k != static_cast<uint>(lb - 1); k--) {
            if (k + cost(i, j, k) <= t)
                return static_cast<int>(k);
        }
        return -1; // infeasible
    }

    int cost_inv_binary_search_rightmost(const uint &i, const uint &j, const uint &t,
                                         const uint &lb, const uint &ub) const {
        // Only use when cost functions are fifo (i.e. arrival time function is non-decreasing)
        assert(lb <= ub);

        // https://en.wikipedia.org/wiki/Binary_search_algorithm#Procedure_for_finding_the_rightmost_element
        uint L(lb), R(ub + 1); // ub+1 for ub to be included

        while (L < R) {
            uint m = (L + R) / 2; // L+R is always positive, i.e. floor((L+R)/2)
            if (m + cost(i, j, m) > t)
                R = m;
            else
                L = m + 1;
        }
        int ret = static_cast<int>(R) - 1;
        // infeasible if ret is lower than the LB
        return ret < static_cast<int>(lb) ? -1 : ret;
    }

    static constexpr CostLB lbs[] = {&CostModel::lb_naive,
                                     &CostModel::lb_medium};
};

template<bool MAKE_FIFO>
class PWConstantCostModel : public CostModel {
    vector<vector<vector<uint>>> mat;
    uint n_timesteps, w_timesteps;
public:
    PWConstantCostModel(const vector<vector<vector<uint>>> &m, const uint &N_,
                        const uint &n_timesteps, const uint &w_timesteps, const uint &start) :
            CostModel(N_ + 1),
            mat(m),
            n_timesteps(n_timesteps),
            w_timesteps(w_timesteps) {

        // Duplicate end depot
        for (uint i = 0; i < N_; i++) {
            mat[i].push_back(mat[i][start]);
        }
        mat.push_back(mat[start]);

        // if costs are made fifo, we must ensure the precondition that the slope is always >= -1.
        if (MAKE_FIFO && !is_fifo())
            force_fifo(); // force the slope to be >= -1
    }

    int cost_inv(const uint &i, const uint &j, const uint &t,
                 const uint &lb, const uint &ub) const {
        if (!MAKE_FIFO)
            return CostModel::cost_inv_linear(i, j, t, lb, ub);

        assert(CostModel::cost_inv_linear(i, j, t, lb, ub) ==
               CostModel::cost_inv_binary_search_rightmost(i, j, t, lb, ub));
        // Equivalent (because instance is FIFO) and faster
        return CostModel::cost_inv_binary_search_rightmost(i, j, t, lb, ub);
    }

    void force_fifo() {
        for (uint i = 0; i < N; i++)
            for (uint j = 0; j < N; j++)
                for (int ts = static_cast<int>(n_timesteps) - 1; ts >= 0; ts--) // ok thanks to the guard
                    this->mat[i][j][ts] = min(this->mat[i][j][ts], this->mat[i][j][ts + 1] + w_timesteps);
    }

    bool is_fifo() {
        for (uint i = 0; i < N; i++)
            for (uint j = 0; j < N; j++)
                for (uint ts = 0; ts < n_timesteps; ts++)
                    if (this->mat[i][j][ts] > this->mat[i][j][ts + 1] + w_timesteps)
                        return false;
        return true;
    }

    uint cost_fifo(const uint &o, const uint &d, const uint &t) const {
        uint ts = t / w_timesteps; // rounds toward 0, i.e floor when positive (always the case).
        assert(ts <= n_timesteps + 1);
        uint x = w_timesteps - (t % w_timesteps);
        return min(mat[o][d][ts], mat[o][d][ts + 1] + x);
    }

    uint cost_raw(const uint &o, const uint &d, const uint &t) const {
        uint ts = t / w_timesteps; // rounds toward 0, i.e floor when positive (always the case).
        assert(ts <= n_timesteps);
        return mat[o][d][ts];
    }

    uint cost(const uint &i, const uint &j, const uint &t) const {
        return MAKE_FIFO ? cost_fifo(i, j, t) : cost_raw(i, j, t);
    }

    template<LBCostType lb>
    void dump_cost_bounds(const uint &i, const uint &j) {
        for (uint t = 0; t < n_timesteps; t++)
            cout << (this->*lb)(i, j, t * w_timesteps) << " ";
        cout << endl;
    }
};

class ConstantCostModel : public CostModel {
    vector<static_array<uint, MAX_SIZE>> mat;
public:
    ConstantCostModel(const vector<vector<uint>> &matrix, const uint &N_, const uint &start,
                      const bool &duplicate_depot = true) :
            CostModel(duplicate_depot ? N_ + 1 : N_),
            mat(N_) {
        for (uint i = 0; i < N_; i++) {
            mat[i].resize(N_);
            for (uint j = 0; j < N_; j++) {
                mat[i][j] = matrix[i][j];
            }
        }
        if (duplicate_depot) {
            // Duplicate end depot & make sure that the diagonal only contains zeros (precondition in heuristics)
            for (uint i = 0; i < N_; i++) {
                mat[i][i] = 0;
                mat[i].push_back(mat[i][start]);
            }
            mat.push_back(mat[start]);
        }
    }

    uint cost(const uint &o, const uint &d, const uint &t) const {
        return mat[o][d];
    }

    int cost_inv(const uint &o, const uint &d, const uint &t,
                 const uint &lb, const uint &ub) const {
        return static_cast<int>(t) - static_cast<int>(mat[o][d]);
        // infeasible if t < mat[o][d] (return value will be negative)
    }

    void precompute_lbs(const vector<TimeWindow> &tws,
                        const vector<uint> &service_times,
                        const vector<Bitset> &E,
                        const vector<Bitset> &E_) override {
        precompute_LDTs(tws, service_times, E, E_);
        for (uint i = 0; i < N; i++)
            for (uint j = 0; j < N; j++)
                lb_naive_[i][j] = mat[i][j];
        for (uint i = 0; i < N; i++)
            for (uint j = 0; j < N; j++)
                lb_medium_[i][j] = mat[i][j];
        for (uint i = 0; i < N - 1; i++) {
            for (const auto &j: E[i]) {
                const uint latest_arrival_time = cost(i, j, 0) + // t does not matter, constant case
                                                 tws[i].late;
                if (tws[j].early > latest_arrival_time) // it is only possible to arrive before e[j]
                    lb_medium_[i][j] += tws[j].early - latest_arrival_time; // incl. waiting time
                assert(lb_medium_[i][j] >= lb_naive_[i][j]);
            }
        }
    }
};

template<bool use_correct_rounding, bool memorize_costs=false>
class IGPCostModel : public CostModel {
    static const uint INFEASIBLE_COST = 1u<<16;
    vector<vector<double>> distances;
    vector<vector<uint>> C;
    vector<pair<uint, uint>> times;
    uint n_zones;
    uint n_timesteps;
    vector<vector<double>> speeds;
    vector<uint> timestep_mapping;
    uint max_t;
	vector<uint16_t> mat;
	uint N;
public:
    IGPCostModel(const IGPInstance &inst, const uint &start) :
            CostModel(inst.N + 1),
            distances(inst.distances),
            C(inst.C),
            times(inst.times),
            n_zones(inst.n_zones),
            n_timesteps(inst.n_timesteps),
            speeds(inst.speeds),
			N(inst.N+1) {
        // Duplicate end depot
        for (uint i = 0; i < inst.N; i++) {
            distances[i].push_back(distances[i][start]);
        }
        distances.push_back(distances[start]);

        uint current_ts = 0;
        for (const auto &p: times) {
            for (uint i = p.first; i < p.second; i++) {
                timestep_mapping.push_back(current_ts);
            }
            ++current_ts;
        }
        max_t = times.back().second;
		if(memorize_costs){
			mat.resize(N*N*max_t);
			for (uint i = 0; i <N; i++)
				for (uint j = 0; j <N; j++)
					for(uint t = 0; t<max_t;t++)
						mat[i*max_t*N+j*max_t+t]=cost_(i,j,t);
		}
    }

    uint cost(const uint &i, const uint &j, const uint &t_0) const override {
		if(memorize_costs)
			return mat[i*max_t*N+j*max_t+t_0];
		else
			return cost_(i, j, t_0);
	}

    uint cost_(const uint &i, const uint &j, const uint &t_0) const {
        const uint c = C[i][j];
        if (use_correct_rounding) {
            // see Penelope's thesis p. 51
            double t = t_0;
            if (t >= max_t)
                return INFEASIBLE_COST;

            uint k = timestep_mapping[t_0];
            double d = distances[i][j];
            double t_ = t + d / speeds[c][k];
            while (t_ > times[k].second) {
                d -= speeds[c][k] * (times[k].second - t);
                t = times[k].second;
                ++k;
                if (k >= n_timesteps)
                    return INFEASIBLE_COST;
                t_ = t + d / speeds[c][k];
            }
            return round(t_ - t_0);
        } else {
            // function adapted from Vu et al.
            uint t = t_0;
            if (t >= max_t)
                return INFEASIBLE_COST;

            uint k = timestep_mapping[t_0];
            double d = distances[i][j];
            uint t_ = t + d / speeds[c][k];
            while (t_ > times[k].second) {
                d -= speeds[c][k] * (times[k].second - t);
                t = times[k].second;
                ++k;
                if (k >= n_timesteps)
                    return INFEASIBLE_COST;
                t_ = t + d / speeds[c][k];
            }
            return t_ - t_0;
        }
    }

    int cost_inv(const uint &i, const uint &j, const uint &t,
                 const uint &lb, const uint &ub) const override {
        assert(CostModel::cost_inv_linear(i, j, t, lb, ub) ==
               CostModel::cost_inv_binary_search_rightmost(i, j, t, lb, ub));
        // Equivalent (because instance is FIFO) and faster
        return CostModel::cost_inv_binary_search_rightmost(i, j, t, lb, ub);
    }

};

#endif //TDTSPTW_COST_MODELS_H
