#ifndef TDTSPTW_STATIC_ARRAY_HPP
#define TDTSPTW_STATIC_ARRAY_HPP

#include <array>
#include <vector>
#include <iostream>
#include "constants.h"

template<typename T, size_t N>
class static_array {
    std::array<T, N> a_;
    uint size_;
public:
    explicit static_array(const uint &n = 0, const T &v = T()) {
        fill(n, v);
    }

    static_array(const std::vector<T> &other) {
        size_ = other.size();
        for (uint i = 0; i < size_; i++)
            a_[i] = other[i];
    }

    void fill(const uint &n, const T &v = T()) {
        size_ = n;
        std::fill(a_.begin(), a_.begin() + n, v);
    }

    void resize(const uint &n) {
        // Note: it's the caller's responsibility to initialize the cells.
        size_ = n;
    }

    void push_back(const T &elem) {
        a_[size_++] = elem;
    }

    template<typename... Args>
    void emplace_back(Args &&... args) {
        a_[size_++] = T(std::forward<Args>(args)...);
    }

    const T &back() const {
        return a_[size_ - 1];
    }

    T &back() {
        return a_[size_ - 1];
    }

    void pop_back() {
        --size_;
    }

    void clear() {
        size_ = 0;
    }

    const uint &size() const {
        return size_;
    }

    const bool empty() const {
        return size_ == 0;
    }

    typedef typename std::array<T, N>::const_iterator const_iterator;

    const_iterator begin() const { return a_.cbegin(); }

    const_iterator end() const { return a_.cbegin() + size_; }

    typedef typename std::array<T, N>::iterator iterator;

    iterator begin() { return a_.begin(); }

    iterator end() { return a_.begin() + size_; }

    T &operator[](std::size_t idx) { return a_[idx]; }

    const T &operator[](std::size_t idx) const { return a_[idx]; }
};

using Path = static_array<uint, MAX_SIZE>;

#endif //TDTSPTW_STATIC_ARRAY_HPP
