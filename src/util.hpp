#ifndef TDTSPTW_UTIL_HPP
#define TDTSPTW_UTIL_HPP

#include <queue>
#include <fstream>
#include <unistd.h> // hostname
#include <sys/time.h> // CPUTime
#include <sys/resource.h> // CPUTime
#include <vector>
#include <queue>
#include <unordered_map>
#include <chrono>
#include "constants.h"
#include "unordered_dense.h"

using namespace std::string_literals;

inline std::string
peak_memory() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("VmPeak:", 0) == 0)
            return line;
    return "VmPeak not found in /proc/self/status";
}

inline std::string hostname() {
    char buf[512];
    if (gethostname(buf, 512) == 0) { // success = 0, failure = -1
        return {buf};
    }
    return {};
}

inline void limit_memory_use(unsigned long lim_gb) {
    struct rlimit memlimit;
    unsigned long bytes = lim_gb * (1024 * 1024 * 1024);
    memlimit.rlim_cur = bytes;
    memlimit.rlim_max = bytes;
    setrlimit(RLIMIT_AS, &memlimit);
}

inline std::pair<double, double> CPUTime() {
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return std::make_pair(static_cast<double>(ru.ru_utime.tv_sec) + static_cast<uint>(ru.ru_utime.tv_usec) * 1e-6,
                          static_cast<double>(ru.ru_stime.tv_sec) + static_cast<uint>(ru.ru_stime.tv_usec) * 1e-6);
}

template<class T, class S, class C>
void shrinkToFit(std::priority_queue<T, S, C> &q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
        static S &Container(std::priority_queue<T, S, C> &q) {
            return q.*&HackedQueue::c;
        }
    };
    HackedQueue::Container(q).shrink_to_fit();
}

template<class T, class S, class C>
void clear(std::priority_queue<T, S, C> &q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
        static S &Container(std::priority_queue<T, S, C> &q) {
            return q.*&HackedQueue::c;
        }
    };
    HackedQueue::Container(q).clear();
}

class stopwatch {
    using tp = std::chrono::time_point<std::chrono::steady_clock>;
    ankerl::unordered_dense::map<std::string, std::tuple<tp, double, double>> m_; // key -> (time point, lap, global)

    template<bool total = false>
    void display(const std::string &prefix = ""s) {
        for (const auto &[s, x]: m_)
            std::cout << prefix << s << "\t" << std::get<total ? 2 : 1>(x) << std::endl;
    }

public:
    void tic(const std::string &s) {
        std::get<0>(m_[s]) = std::chrono::steady_clock::now();
    }

    void toc(const std::string &s) {
        const auto it = m_.find(s);
        std::get<1>(it->second) += std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::steady_clock::now() - std::get<0>(it->second)).count();
    }

    void display_total(const std::string &prefix = ""s) {
        display<true>(prefix);
    }

    void display_lap(const std::string &prefix = ""s) {
        display<false>(prefix);
    }

    void clear() {
        m_.clear();
    }

    void lap() {
        for (auto &[s, x]: m_) {
            std::get<2>(x) += std::get<1>(x);
            std::get<1>(x) = 0;
        }
    }
};

#endif //TDTSPTW_UTIL_HPP
