#ifndef TDTSPTW_MSA_HPP
#define TDTSPTW_MSA_HPP

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "util.hpp"
#include <cstring> // memcpy

using namespace std;


template<uint MAX_SIZE>
class MSA {

    struct UnionFind {
        static_array<int, MAX_SIZE> p;

        explicit UnionFind(const uint &n) : p(n, -1) {};

        uint find(const uint &u) {
            if (p[u] < 0)
                return u;
            p[u] = static_cast<int>(find(static_cast<uint>(p[u])));
            return static_cast<uint>(p[u]);
        }

        bool merge(uint u, uint v) {
            if ((u = find(u)) == (v = find(v)))
                return false;
            if (p[u] > p[v])
                swap(u, v);
            p[u] += p[v]; // Prefer wide trees to deep trees.
            p[v] = static_cast<int>(u);
            return true;
        }
    };

    struct Node {
        Node *ch[2]{};
        uint src, dst, weight;
        uint delta = 0;

        Node() {}

        explicit Node(const uint &src, const uint &dst, const uint &w) : src(src), dst(dst), weight(w) {}

        ~Node() {}
    };


    Node *new_node(const uint &src, const uint &dst, const uint &w) {
        // Works well because all allocations are done in the beginning (no fragmentation)
        nodes_alloc[curr_alloc] = Node(src, dst, w);
        return &nodes_alloc[curr_alloc++];
    }

    Node *new_nodes(const uint &n) {
        uint old_alloc = curr_alloc;
        curr_alloc += n;
        return &nodes_alloc[old_alloc];
    }

    struct SkewHeap {
        Node *root;

        SkewHeap() : root(nullptr) {}

        ~SkewHeap() {}

        void clear() {
            root = nullptr;
        }

        SkewHeap &operator=(SkewHeap &&obj) noexcept {
            root = obj.root;
            obj.root = nullptr;
            return *this;
        }

        SkewHeap(const SkewHeap &p1) = delete;

        void propagate(Node *a) {
            // Updates a's weight (lazily for the subtree)
            a->weight += a->delta;
            if (a->ch[0])
                a->ch[0]->delta += a->delta;
            if (a->ch[1])
                a->ch[1]->delta += a->delta;
            a->delta = 0;
        }

        template<bool prop = true>
        Node *merge(Node *a, Node *b) {
            if (!a || !b)
                return a ? a : b;
            if (prop) {
                propagate(a);
                propagate(b);
            }
            if (a->weight > b->weight)
                swap(a, b);
            a->ch[1] = merge<prop>(b, a->ch[1]);
            swap(a->ch[0], a->ch[1]);
            return a;
        }

        void push(Node *node) {
            // assumes that all nodes have delta=0, propagation is disabled.
            root = merge<false>(root, node);
        }

        void pop() {
            propagate(root);
            root = merge(root->ch[0], root->ch[1]);
        }

        const Node &top() const {
            //propagate(root);
            // Does not propagate. weight of the node is equal to root->weight+root->delta
            return *root;
        }

        bool empty() const {
            return !root;
        }

        void add(const uint &delta) {
            root->delta += delta;
        }

        void merge(SkewHeap &x) {
            root = merge(root, x.root);
        }
    };

    uint n_;
    uint curr_alloc;
    array<SkewHeap, MAX_SIZE> heaps{};
    array<Node, MAX_SIZE * MAX_SIZE> nodes_alloc;
public:

    explicit MSA(const uint &n) {
        init(n);
    }

    void init(const uint &n) {
        n_ = n;
        curr_alloc = 0;
        for (auto &h: heaps)
            h.clear();
    }

    static MSA<MAX_SIZE> *instances_[MAX_THREAD_NUMBER];

    static MSA<MAX_SIZE> *get_instance(const uint &i = 0) {
        if (MSA<MAX_SIZE>::instances_[i] == nullptr)
            MSA<MAX_SIZE>::instances_[i] = new MSA<MAX_SIZE>(0);
        return MSA<MAX_SIZE>::instances_[i];
    }

    void add_edge(const uint &src, const uint &dst, const uint &weight) {
        heaps[dst].push(new_node(src, dst, weight));
    }

    // Row contains [(src, weight), ...] - in min-heap order relative to weight.
    void add_edges(const uint &dst, const static_array<pair<uint, uint>, MAX_SIZE> &row) {
        const uint &sz = row.size();
        Node *nodes = new_nodes(sz);

        {
            Node *cur_node = nodes;
            for (uint i = 0; i < sz; i++, cur_node++) {
                cur_node->delta = 0;
                cur_node->ch[0] = nullptr;
                cur_node->ch[1] = nullptr;
            }
        }

        Node *cur_node = nodes;
        for (uint i = 0; i < sz; i++, cur_node++) {
            const auto &[src, w] = row[i];
            const uint left = (i << 1) | 1;

            cur_node->src = src;
            cur_node->dst = dst;
            cur_node->weight = w;
            if (left < sz) {
                cur_node->ch[0] = nodes + left;
                if (left + 1 < sz)
                    cur_node->ch[1] = nodes + left + 1;
            }
        }

        heaps[dst].push(nodes);
    }

    uint solve(const uint &root, const uint &ub, const Bitset &nodes) {
        UnionFind uf(n_);
        uint score = 0;
        static_array<int, MAX_SIZE> seen(n_, -1);
        seen[root] = static_cast<int>(root); // seen : if <0, has no incoming arc & needs it

        for (const auto &s: nodes) { // For each node (arbitrary order)
            static_array<uint, MAX_SIZE> path;
            for (uint u = s;
                 seen[u] < 0;) { // for each node that has no selected incoming arc & needs one
                path.push_back(u);
                seen[u] = static_cast<int>(s);// u was reached from s
                if (heaps[u].empty())
                    return ub; // can't reach anything else - no solution

                // Use the cheapest edge coming to u & increment score
                const Node &min_e = heaps[u].top();
                const uint src(min_e.src);
                const uint weight = min_e.weight + min_e.delta; // avoid one propagate operation
                score += weight;

                // decrease all keys by min_e's cost
                heaps[u].add(-weight);
                heaps[u].pop(); // pop this edge (and propagate delta)

                const uint v = uf.find(src);
                if (seen[v] == static_cast<int>(s)) { // cycle detected
                    SkewHeap new_heap;
                    while (true) { // merge current path to a supernode
                        const uint w = path.back();
                        path.pop_back();
                        new_heap.merge(heaps[w]); // clears heaps[w]
                        if (!uf.merge(v, w)) // until they both belong to the same tree
                            break;
                    }
                    heaps[uf.find(v)] = std::move(new_heap); // incoming edges to the supernode (with updated costs)
                    seen[uf.find(v)] = -1; // the supernode has to be reached again
                }
                u = uf.find(v); // next = v if it has no parent, otherwise its parent
            }
        }
        return score;
    }
};

template<uint MAX_SIZE>
MSA<MAX_SIZE> *MSA<MAX_SIZE>::instances_[] = {};


#endif //TDTSPTW_MSA_HPP
