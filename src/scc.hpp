#ifndef TDTSPTW_SCC_HPP
#define TDTSPTW_SCC_HPP

class SCC {
public:
    using Vbs = static_array<Bitset, MAX_SIZE>;
    using Vuint = static_array<uint, MAX_SIZE>;

    static void tarjan(const Vbs &adj,
                       Vbs &partitions, Vuint &partitions_t) {
        uint cur_num = 0;
        Vuint P;

        Vuint num(adj.size(), ~0u);
        Vuint num_reachable(adj.size());

        Bitset in_P;

        for (uint v = 0; v < adj.size(); v++) {
            if (num[v] == ~0u)
                tarjan_dfs(v, adj, cur_num, num, num_reachable, in_P, P, partitions, partitions_t);
        }
    }

    static void topological_sort(const uint &i,
                             const Vbs &cfc_adj,
                             Bitset &marked,
                             Vuint &L) {
        dfs_topo(i, cfc_adj, marked, L);
        reverse(L.begin(), L.end());
    }

private:
    static void tarjan_dfs(const uint &v,
                           const Vbs &adj,
                           uint &cur_num,
                           Vuint &num,
                           Vuint &num_reachable,
                           Bitset &in_P,
                           Vuint &P,
                           Vbs &partitions,
                           Vuint &partitions_t) {
        num[v] = cur_num;
        num_reachable[v] = cur_num;
        ++cur_num;
        P.push_back(v);
        in_P.add(v);

        for (const auto &w: adj[v]) {
            if (num[w] == ~0u) {
                tarjan_dfs(w,
                           adj, cur_num, num, num_reachable, in_P, P, partitions, partitions_t);
                num_reachable[v] = min(num_reachable[v], num_reachable[w]);
            } else if (in_P[w]) {
                num_reachable[v] = min(num_reachable[v], num[w]);
            }
        }
        if (num_reachable[v] == num[v]) {
            partitions.push_back(Bitset());
            uint w;
            do {
                w = P.back();
                P.pop_back();
                in_P.remove(w);
                partitions.back().add(w);
                partitions_t[w] = partitions.size() - 1;
            } while (w != v);
        }
    }

    static void dfs_topo(const uint &i,
                         const Vbs &cfc_adj,
                         Bitset &marked,
                         Vuint &L) {
        if (marked[i])
            return;
        marked.add(i);
        for (const auto &j: cfc_adj[i])
            dfs_topo(j, cfc_adj, marked, L);
        L.push_back(i);
    }
};

#endif //TDTSPTW_SCC_HPP
