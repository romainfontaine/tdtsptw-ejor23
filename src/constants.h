#ifndef TDTSPTW_CONSTANTS_H
#define TDTSPTW_CONSTANTS_H

#include <iostream>
#include <chrono>


#ifndef MAX_INSTANCE_SIZE
const uint MAX_SIZE = 64;
#else
const uint MAX_SIZE = MAX_INSTANCE_SIZE;
#endif
const uint MAX_THREAD_NUMBER = 64;

using Action = uint;

using real_time_point = std::chrono::time_point<std::chrono::steady_clock>;


inline void ensure_max_size(const uint &N) {
    if (N + 1 > MAX_SIZE) {
        std::cerr << "Error: N+1>MAX_SIZE (" << N + 1 << ">" << MAX_SIZE << ")" << std::endl;
        exit(1);
    }
}

enum class LBCostType {
    naive = 0, medium = 1
};

enum class LDTFilteringMode {
    NONE = 0, CONDITION = 1, PRECOMPUTED = 2
};

struct TimeWindow {
    uint early;
    uint late;

    explicit TimeWindow(const uint &e = 0, const uint &l = 0) : early(e), late(l) {};
};


#endif //TDTSPTW_CONSTANTS_H
