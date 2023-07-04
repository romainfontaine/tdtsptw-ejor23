#ifndef TDTSPTW_BITSETS_HPP
#define TDTSPTW_BITSETS_HPP

#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <bitset>
#include <unordered_map>
#include <climits>
#include "constants.h"

using namespace std;

template<typename T>
struct BSIterator // input iterator
{
    explicit BSIterator(const T &bs) : _bs(bs), _c(bs.ctz()) {}

    // Prefix increment
    BSIterator &operator++() {
        _bs.flip(_c);
        _c = _bs.ctz();
        return *this;
    }

    const uint &operator*() const {
        return _c;
    }

private:
    T _bs;
    uint _c;

    friend bool operator!=(const BSIterator &a, const BSIterator &b) {
        return a._bs != b._bs;
    };

    friend bool operator==(const BSIterator &a, const BSIterator &b) {
        return a._bs == b._bs;
    };
};

template<typename T>
struct std::iterator_traits<BSIterator<T>> {
    typedef std::input_iterator_tag iterator_category;
};

class StaticBitsetIterator {
    const uint *ptr;
public:
    explicit StaticBitsetIterator(const uint *p) : ptr(p) {}

    StaticBitsetIterator &operator++() {
        ptr++;
        return *this;
    }

    const uint &operator*() {
        return *ptr;
    }

    friend bool operator!=(const StaticBitsetIterator &a, const StaticBitsetIterator &b) {
        return a.ptr != b.ptr;
    };
};

template<uint MAX_SIZE>
class StaticBitsetEnumerator {
    std::array<uint, MAX_SIZE> array;
    uint size;
public:
    template<typename BSType>
    explicit StaticBitsetEnumerator(const BSType &bs) :size(0) {
        for (const uint &v: bs)
            array[size++] = v;
    }

    StaticBitsetEnumerator(const StaticBitsetEnumerator &) = delete;

    StaticBitsetEnumerator(StaticBitsetEnumerator &&) = delete;

    StaticBitsetIterator begin() const { return StaticBitsetIterator(array.data()); }

    StaticBitsetIterator end() const { return StaticBitsetIterator(array.data() + size); }
};

class Bitset64 {
    uint64_t bs;

    Bitset64(const uint64_t &bitset) : bs(bitset) {}

public:

    Bitset64() : Bitset64(0ull) {}

    const uint64_t &get_bs() const {
        return bs;
    }

    static Bitset64 full_set(const uint &i) {
        return Bitset64((1ull << i) - 1);
    }

    static Bitset64 unit_set(const uint &i) {
        return Bitset64(1ull << i);
    }

    static Bitset64 from_uint(const uint &i) {
        return Bitset64(i);
    }

    Bitset64 &add(const uint &i) {
        bs |= Bitset64::unit_set(i).bs;
        return *this;
    }

    Bitset64 &remove(const uint &i) {
        bs &= ~Bitset64::unit_set(i).bs;
        return *this;
    }

    void flip(const uint &i) {
        bs ^= Bitset64::unit_set(i).bs;
    }

    Bitset64 operator|(const Bitset64 &o) const {
        return Bitset64(bs | o.bs);
    }

    Bitset64 &operator|=(const Bitset64 &rhs) {
        this->bs |= rhs.bs;
        return *this;
    }

    Bitset64 &operator&=(const Bitset64 &rhs) {
        this->bs &= rhs.bs;
        return *this;
    }

    Bitset64 operator&(const Bitset64 &o) const {
        return Bitset64(bs & o.bs);
    }

    Bitset64 operator~() const {
        return Bitset64(~bs);
    }

    bool is_empty() const {
        return !bs;
    }

    bool any() const {
        return !is_empty();
    }

    void clear() {
        bs = 0;
    }

    uint count() const {
        return static_cast<uint>(__builtin_popcountll(bs));
    }

    bool operator!=(const Bitset64 &o) const {
        return bs != o.bs;
    }

    bool operator==(const Bitset64 &o) const {
        return bs == o.bs;
    }

    uint ctz() const {
        // Undefined if _bs == 0 - be careful
        return static_cast<uint>(__builtin_ctzll(bs));
    }

    inline bool contains(const uint &i) const {
        return bs & Bitset64::unit_set(i).bs;
    }

    inline bool operator[](const uint &i) const {
        return contains(i);
    }

    BSIterator<Bitset64> begin() const { return BSIterator<Bitset64>(*this); }

    BSIterator<Bitset64> end() const { return BSIterator<Bitset64>(Bitset64()); }

    friend ostream &operator<<(ostream &, const Bitset64 &);
};

inline ostream &operator<<(ostream &o, const Bitset64 &b) {
    o << "Bitset(" << bitset<64>(b.bs) << ')';
    return o;
}

template<uint MAX_SIZE>
class LargeBitset;

template<uint MAX_SIZE>
ostream &operator<<(ostream &, const LargeBitset<MAX_SIZE> &);

template<uint MAX_SIZE>
class LargeBitset {
public:
    static constexpr uint BITS_PER_ELEM = sizeof(uint64_t) * CHAR_BIT;
    static constexpr uint ARRAY_SIZE = MAX_SIZE / BITS_PER_ELEM;
private:
    uint64_t bs[ARRAY_SIZE] = {};
public:
    LargeBitset() {}

    const uint64_t &get_bs(const uint &i) const {
        return this->bs[i];
    }

    static LargeBitset unit_set(const uint &i) {
        LargeBitset bs;
        return bs.add(i);
    }

    static LargeBitset full_set(const uint &N) {
        LargeBitset bs;
        const uint k = N / BITS_PER_ELEM;
        for (uint i = 0; i < k; i++) {
            bs.bs[i] = ~0ul;
        }
        bs.bs[k] = (1ull << (N % BITS_PER_ELEM)) - 1;
        return bs;
    }

    LargeBitset operator|(const LargeBitset &o) const {
        LargeBitset result;
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            result.bs[i] = this->bs[i] | o.bs[i];
        }
        return result;
    }

    LargeBitset operator&(const LargeBitset &o) const {
        LargeBitset result;
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            result.bs[i] = this->bs[i] & o.bs[i];
        }
        return result;
    }

    LargeBitset &operator|=(const LargeBitset &rhs) {
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            bs[i] |= rhs.bs[i];
        }
        return *this;
    }

    LargeBitset &operator&=(const LargeBitset &rhs) {
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            bs[i] &= rhs.bs[i];
        }
        return *this;
    }

    LargeBitset operator~() const {
        LargeBitset result;
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            result.bs[i] = ~this->bs[i];
        }
        return result;
    }

    bool is_empty() const {
        bool empty = true;
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            empty = empty && (this->bs[i] == 0ull);
        }
        return empty;
    }

    bool any() const {
        return !is_empty();
    }

    void clear() {
        for (uint i = 0; i < ARRAY_SIZE; i++)
            this->bs[i] = 0;
    }

    uint count() const {
        uint sum(0);
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            sum += static_cast<uint>(__builtin_popcountll(bs[i]));
        }
        return sum;
    }

    bool operator==(const LargeBitset &o) const {
        bool same = true;
        for (uint i = 0; i < ARRAY_SIZE; i++) {
            same = same && bs[i] == o.bs[i];
        }
        return same;
    }

    bool operator!=(const LargeBitset &o) const {
        return !(*this == o);
    }

    LargeBitset<MAX_SIZE> &add(const uint &i) {
        bs[i / BITS_PER_ELEM] |= (1ull << (i % BITS_PER_ELEM));
        return *this;
    }

    LargeBitset<MAX_SIZE> &remove(const uint &i) {
        bs[i / BITS_PER_ELEM] &= ~(1ull << (i % BITS_PER_ELEM));
        return *this;
    }

    bool contains(const uint &i) const {
        return bs[i / BITS_PER_ELEM] & (1ull << (i % BITS_PER_ELEM));
    }

    inline bool operator[](const uint &i) const {
        return contains(i);
    }

    void flip(const uint &i) {
        bs[i / BITS_PER_ELEM] ^= (1ull << (i % BITS_PER_ELEM));
    }

    uint ctz() const {
        uint i = 0;
        while (bs[i] == 0ull && i < ARRAY_SIZE - 1) { // get the first bs != 0 (or the last one, if all are == 0)
            ++i;
        }
        // __builtin_ctzll may be called when bs[i] == 0, but the value will not be read
        return i * BITS_PER_ELEM + static_cast<uint>(__builtin_ctzll(bs[i]));
    }

    BSIterator<LargeBitset<MAX_SIZE>> begin() const { return BSIterator<LargeBitset<MAX_SIZE>>(*this); }

    BSIterator<LargeBitset<MAX_SIZE>> end() const { return BSIterator<LargeBitset<MAX_SIZE>>(LargeBitset<MAX_SIZE>()); }

    friend ostream &operator
    <<<>(ostream &, const LargeBitset &);
};

template<uint MAX_SIZE>
ostream &operator<<(ostream &o, const LargeBitset<MAX_SIZE> &b) {
    o << "LargeBitset<" << MAX_SIZE << ">(";
    for (int i = LargeBitset<MAX_SIZE>::ARRAY_SIZE - 1; i > 0; i--) {
        o << bitset<LargeBitset<MAX_SIZE>::BITS_PER_ELEM>(b.bs[i]) << " ";
    }
    o << bitset<LargeBitset<MAX_SIZE>::BITS_PER_ELEM>(b.bs[0]) << ")";
    return o;
}

namespace std {
    template<>
    struct hash<Bitset64> {
        inline std::size_t operator()(const Bitset64 &b) const {
            return b.get_bs();
        }
    };

    template<uint MAX_SIZE>
    struct hash<LargeBitset<MAX_SIZE>> {
        std::size_t operator()(const LargeBitset<MAX_SIZE> &b) const {
            size_t tmp(0);
            for (uint i = 0; i < LargeBitset<MAX_SIZE>::ARRAY_SIZE; i++) {
                tmp ^= b.get_bs(i);
            }
            return tmp;
        }
    };
}

using Bitset = conditional<MAX_SIZE <= 64, Bitset64, LargeBitset<MAX_SIZE>>::type;

#endif //TDTSPTW_BITSETS_HPP
