#ifndef TDTSPTW_INSTANCES_H
#define TDTSPTW_INSTANCES_H

#include <fstream>
#include <vector>
#include <sstream>
#include "constants.h"

using namespace std;

struct ConstantInstance {
    uint N;
    vector<vector<uint>> mat;
    vector<TimeWindow> tws;

    explicit ConstantInstance(const string &filename) {
        // Precondition : file exists.
        ifstream f(filename);
        uint n;
        f >> n;
        N = n;
        ensure_max_size(N);
        mat = vector<vector<uint >>(n, vector<uint>(n));
        tws = vector<TimeWindow>(n);

        for (uint i = 0; i < n; i++)
            for (uint j = 0; j < n; j++)
                f >> mat[i][j];

        for (uint i = 0; i < n; i++)
            f >> tws[i].early >> tws[i].late;

        if (!f.good()) {
            cerr << "Parse error for file " << filename << endl;
            exit(1);
        }
    }
};


struct IGPInstance {
    uint N;
    vector<vector<double>> distances;
    vector<TimeWindow> tws;
    vector<vector<uint>> C;
    vector<pair<uint, uint>> times;
    uint n_zones;
    uint n_timesteps;
    vector<vector<double>> speeds;
    const bool use_correct_rounding;
    const bool memorize_costs;

    IGPInstance(const string &main_file, const string &new_jams_file,
                const bool &use_correct_rounding, const bool &use_new_jams,
                const bool &memorize_costs) :
            use_correct_rounding(use_correct_rounding), memorize_costs(memorize_costs) {
        // Precondition : main_file exists, and new_jams_file exists (if use_new_jams)
        ifstream f(main_file);

        f >> N;
        ensure_max_size(N);

        distances = vector<vector<double >>(N, vector<double>(N));
        tws = vector<TimeWindow>(N);
        C = vector<vector<uint >>(N + 1, vector<uint>(N + 1));

        for (uint i = 0; i < N; i++)
            for (uint j = 0; j < N; j++)
                f >> distances[i][j];

        for (uint i = 0; i < N; i++)
            f >> tws[i].early >> tws[i].late;

        {
            string line;
            getline(f, line); // ignore last \r\n
            getline(f, line); // ignore C:
        }

        n_zones = 0;
        for (uint i = 0; i < N + 1; i++)
            for (uint j = 0; j < N + 1; j++) {
                f >> C[i][j];
                if (C[i][j] > 0)
                    C[i][j]--;
                n_zones = max(n_zones, C[i][j]);
            }
        n_zones++;

        {
            string line;
            getline(f, line); // ignore last \r\n
            getline(f, line); // ignore Time:
        }

        string line;
        while (getline(f, line) && line.find("Speed:") != 0) {
            uint i, j;
            istringstream(line) >> i >> j;
            times.emplace_back(i, j);
        }

        n_timesteps = times.size();
        speeds = vector<vector<double >>(n_zones, vector<double>(n_timesteps));
        for (uint z = 0; z < n_zones; z++)
            for (uint t = 0; t < n_timesteps; t++)
                f >> speeds[z][t];

        if (!f.good()) {
            cerr << "Parse error for file " << main_file << endl;
            exit(1);
        }

        if (use_new_jams) {
            f.close();
            f.open(new_jams_file);
            auto new_jams = vector<vector<double >>(n_zones, vector<double>(n_timesteps));
            for (uint z = 0; z < n_zones; z++)
                for (uint t = 0; t < n_timesteps; t++)
                    f >> new_jams[z][t];

            if (!f.good()) {
                cerr << "Parse error for file " << new_jams_file << endl;
                exit(1);
            }

            // Compute new speeds
            for (uint z = 0; z < n_zones; z++)
                for (uint t = 0; t < n_timesteps; t++)
                    speeds[z][t] *= new_jams[z][t];
        }
    }
};

struct PWConstantInstance {
    uint N, w_ts, n_ts;
    bool force_fifo;
    vector<vector<vector<uint>>>
            mat;
    vector<TimeWindow> tws;
    vector<uint> s;

    PWConstantInstance(const uint &n, const uint &w_ts, const uint &n_ts, const bool &force_fifo,
                       const string &cost_file, const string &tw_file, const string &st_file) :
            N(n), w_ts(w_ts), n_ts(n_ts), force_fifo(force_fifo),
            mat(N, vector<vector<uint >>(N, vector<uint>(n_ts))),
            tws(N),
            s(N) {
        ensure_max_size(n);
        // Precondition : all files exist
        ifstream f_costs(cost_file);

        for (uint i = 0; i < N; i++)
            for (uint j = 0; j < N; j++) {
                for (uint t = 0; t < n_ts; t++)
                    f_costs >> mat[i][j][t];
                mat[i][j].push_back(mat[i][j].back()); // guard to make later computation easier
            }

        if (!f_costs.good()) {
            cerr << "Parse error for file " << cost_file << endl;
            exit(1);
        }

        ifstream f_tws(tw_file);
        vector<TimeWindow> time_windows(N);
        for (uint i = 0; i < N; i++) {
            uint e, l;
            f_tws >> e >> l;
            tws[i] = TimeWindow(e, l);
        }

        if (!f_tws.good()) {
            cerr << "Parse error for file " << tw_file << endl;
            exit(1);
        }

        ifstream f_service_times(st_file);
        for (uint i = 0; i < N; i++) {
            f_service_times >> s[i];
        }

        if (!f_service_times.good()) {
            cerr << "Parse error for file " << st_file << endl;
            exit(1);
        }
    }
};

#endif //TDTSPTW_INSTANCES_H
