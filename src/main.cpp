#include "cost_models.h"
#include "model.h"
#include <chrono>
#include <queue>
#include "omp.h"
#include "arg_parser.hpp"
#include "unordered_dense.h"
#include "solving_algorithms.h"

void interactive(TDTSPTW &t, const bool &preprocess_tws) {
    if (t.process_tws_and_compute_lbs(preprocess_tws, true)) // DO compute cost LBs
        exit(0); // instance is infeasible
    State s = t.initial_state();
    while (!t.is_terminal(s)) {
        cout << s << " - Obj=" << s.t << " - Actions : " << endl;
        const Bitset as = t.A(s);
        if (as.is_empty()) {
            cout << "No feasible action available" << endl;
            exit(0);
        }
        for (const auto &a: as) {
            const State s_ = t.tau(s, a);
            cout << a << " (" << t.tws[a].early << ";" << t.tws[a].late << ") " <<
                 "cost=" << s_.t - s.t
                 << ", h_msa=" << t.h_msa<LBCostType::naive>(s_)
                 << ", h_oa=" << t.h_outgoing_arcs<LBCostType::naive>(s_) << endl;
        }
        uint a;
        cin >> a;
        if (!as.contains(a)) {
            cout << "infeasible action" << endl;
            exit(0);
        }
        const State s_ = t.tau(s, a);
        cout << "cost of action " << a << " = " << (s_.t - s.t) << endl;
        s = s_;
    }
    cout << "Done, obj=" << s.t << endl;
}

#ifdef _OPENMP
const bool OPENMP(true);
#else
const bool OPENMP(false);
#endif
#ifdef NDEBUG
const bool N_DEBUG(true);
#else
const bool N_DEBUG(false);
#endif

uint set_nthreads(const uint &n_threads) {
#ifdef _OPENMP
    if (n_threads == 0)
        omp_set_num_threads(omp_get_max_threads() / 2);
    else
        omp_set_num_threads(static_cast<int>(n_threads));
    const uint effective_threads(static_cast<uint>(omp_get_max_threads()));
    if (effective_threads > MAX_THREAD_NUMBER) {
        cerr << "Cannot use " << effective_threads << " threads as MAX_THREAD_NUMBER=" << MAX_THREAD_NUMBER << endl;
        exit(1);
    }
    return effective_threads;
#else
    cerr << "Cannot use parallel computation as program was not compiled with OpenMP" << endl;
    exit(1);
#endif
}

const unordered_map<string, Heuristic> hs = {
        {"0",            &TDTSPTW::h_seq<&TDTSPTW::h_0_out_feasibility_checks<false, LBCostType::medium>>},
        {"oia",          &TDTSPTW::h_seq<&TDTSPTW::h_outgoing_incoming_arcs<LBCostType::medium>>},
        {"oia-ldt",      &TDTSPTW::h_seq<&TDTSPTW::h_outgoing_incoming_arcs<LBCostType::medium, LDTFilteringMode::NONE>>},
        {"msa",          &TDTSPTW::h_seq<&TDTSPTW::h_msa<LBCostType::medium>>},
        {"msa_par",      &TDTSPTW::h_par<&TDTSPTW::h_msa<LBCostType::medium>>},
};

TDTSPTW *parseIGP(ArgParser &a, const uint &start) {
    string main_file, new_jams_file;
    const bool use_correct_rounding(a.getArgument("-correct-rounding", true));
    const bool use_new_jams(a.getArgument("-use-new-jams", true));
    const bool memorize_costs(a.getArgument("-memorize-costs", false));

    const uint args = a.remainingPositionals();
    if (args == 5) {
        uint n = a.getPositional<uint>("n");
        uint delta = a.getPositional<uint>("delta");
        char pattern = a.getPositional<char>("pattern", {'A', 'B'});
        uint width = a.getPositional<uint>("width");
        uint instance = a.getPositional<uint>("instance");

        string benchmark_dir = a.getArgument<string>("-benchmark-dir", "./benchmarks/Ari18Vu20/");
        const string base(benchmark_dir + to_string(n - 1) + "/Delta" +
                          to_string(delta) + "/Pattern" + pattern +
                          "/Width_" + to_string(width) + "/");

        char letter = 'A' + (instance - 1) / 10;
        main_file = base + to_string(n - 1) + letter + "nodi_" + to_string((instance - 1) % 10 + 1) + "_tw_new.txt";
        a.ensureFileExists(main_file, "file");
        if (use_new_jams) {
            new_jams_file = base + "new_Jams_" + to_string(delta) + "_" + (pattern == 'A' ? "C1" : "C3");
            a.ensureFileExists(new_jams_file, "jam_file");
        }
    } else if (args >= 1) {
        main_file = a.getPositionalFile("file");
        if (use_new_jams)
            new_jams_file = a.getPositionalFile("jam_file");
    } else {
        a.parsingError("unable to parse IGP benchmark arguments");
    }
    cout << "# use_new_jams=" << use_new_jams << ", use_correct_rounding=" << use_correct_rounding <<
         ", memorize_costs=" << memorize_costs << endl;
    return new TDTSPTW(IGPInstance(main_file, new_jams_file, use_correct_rounding,
                                   use_new_jams, memorize_costs), start);
}

TDTSPTW *parseConstant(ArgParser &a, const uint &start) {
    const string &file = a.getPositionalFile("file");
    const bool &triangle_ineq = a.getArgument<bool>("-triangle-inequality", false);
    return new TDTSPTW(ConstantInstance(file), start, triangle_ineq);
}

TDTSPTW *parsePWConstant(ArgParser &a, const uint &start) {
    const auto &strictly_positive = [](const uint &x) { return x > 0; };

    const uint &N = a.getPositional<uint>("n", strictly_positive,
                                          "n must be strictly positive");
    const string &cost_file = a.getPositionalFile("cost_file");
    const string &tw_file = a.getPositionalFile("tw_file");
    const string &st_file = a.getPositionalFile("service_times_file");
    const bool &make_fifo = a.getArgument<bool>("-force-fifo", true);
    const auto &ts_duration = a.getArgument<uint>("-timestep-duration", 360, false, strictly_positive,
                                                  "-timestep-duration must be strictly positive");
    const auto &ts_number = a.getArgument<uint>("-timestep-number", 120, false, strictly_positive,
                                                "-timestep-number must be strictly positive");
    cout << "# n=" << N << ", make_fifo=" << make_fifo << ", timestep_duration="
         << ts_duration << ", timestep_number=" << ts_number << endl;
    return new TDTSPTW(PWConstantInstance(N, ts_duration, ts_number, make_fifo, cost_file, tw_file, st_file),
                       start);
}

const string HELP = "See file README.org for help.";

int main(int argc, char *argv[]) {
    real_time_point starting_time = chrono::steady_clock::now();
    ArgParser a("tdtsptw", "1.0.0", HELP, argv[0]);
    cout << "# ";
    for (int i = 0; i < argc; i++) {
        cout << argv[i] << " ";
    }
    cout << endl;

    a.parse(argc - 1, argv + 1);

    if (a.getArgument("-help", false)) {
        a.help(true);
        exit(0);
    }
    {
        ostringstream info;
        time_t now = time(nullptr);
        char *dt = ctime(&now);
        info << strtok(dt, "\n"); // remove trailing '\n'

        info << ", MAX_SIZE=" << MAX_SIZE;
        info << ", NDEBUG=" << N_DEBUG;
        info << ", omp_max_threads=";
        if (OPENMP)
            info << omp_get_max_threads();
        else
            info << "None";
        info << ", hostname=" + hostname();
        cout << "# " << info.str() << endl;
    }
    const bool preprocess_tws = a.getArgument("-preprocess-tws", true);
    const uint mem_limit_gb(a.getArgument("-mem-limit-gb", 64u));

    limit_memory_use(mem_limit_gb);
    ostringstream info;
    info << "preprocess-tws=" << preprocess_tws;
    info << ", mem-limit-gb=" << mem_limit_gb;

    const auto &subcommand(a.getPositional<string>("subcommand", {"solve", "preprocess"}));

    info << ", subcommand=" << subcommand;
    string algorithm;
    if (subcommand == "solve") {
        algorithm = a.getPositional<string>("algorithm",
                                            {"ACS", "BU", "interactive"});
        info << ", algorithm=" << algorithm;
    }

    TDTSPTW *t = nullptr;
    const uint start(0);

    const auto instance_type(a.getPositional<string>("instance_type", {"IGP", "constant", "pw_constant"}));
    info << ", instance_type=" << instance_type;
    cout << "# " << info.str() << endl;
    if (instance_type == "IGP") {
        t = parseIGP(a, start);
    } else if (instance_type == "constant") {
        t = parseConstant(a, start);
    } else if (instance_type == "pw_constant") {
        t = parsePWConstant(a, start);
    } else {
        exit(1);
    }

    if (subcommand == "solve" && algorithm == "ACS") {
        info = ostringstream();
        string h(a.getArgument("-h", "oia"s, false));
        const uint tl(a.getArgument("-tl", 60u));
        const bool ls(a.getArgument("-ls", true));
        const bool process_tws(a.getArgument("-process-tws", true));
        const bool greedy_ub(a.getArgument("-greedy-ub", true));
        const bool gub_ls_twp(a.getArgument("-gub-ls-twp", false)); // force ls and twp, only on the greedy upper bound.

        if (!greedy_ub && gub_ls_twp)
            a.parsingError("-greedy-ub must be set in order to use -gub-ls-twp");

        info << "h=" << h << ", tl=" << tl << ", ls=" << ls << ", process_tws=" << process_tws;

        if (h == "fea") h = "0"; // Alias for fea

        if (h.find("_par") != string::npos) {
            const uint threads(set_nthreads(a.getArgument<uint>("-threads", 0u, false)));
            info << ", threads=" << threads;
        }
        if (hs.count(h) == 0)
            a.parsingError("unknown heuristic");

        if (algorithm == "ACS") {
            const uint w(a.getArgument<uint>("-w", 1u, false,
                                             [](const uint &x) { return x > 0; }, "w must be strictly positive"));
            info << ", w=" << w;
            const bool &fewer_tl_checks = a.getArgument<bool>("-fewer-tl-checks", false);
            info << ", fewer-tl-checks=" << fewer_tl_checks;
            cout << "# " << info.str() << endl;
            a.done();
	    ACS(hs.at(h), *t, w, greedy_ub, gub_ls_twp, ls, preprocess_tws, process_tws, tl,
		fewer_tl_checks, false, starting_time);

        }
    } else if (subcommand == "solve" && algorithm == "BU") {
        const bool use_heuristic(a.getArgument("-use-heuristic", false));
        const bool compute_path(a.getArgument("-compute-path", true));
        const int ub(a.getArgument("-ub", -1));
        a.done();
        if (use_heuristic)
            bottom_up<true>(*t, preprocess_tws, ub, starting_time, compute_path);
        else
            bottom_up<false>(*t, preprocess_tws, ub, starting_time, compute_path);
    } else if (subcommand == "solve" && algorithm == "interactive") {
        a.done();
        interactive(*t, preprocess_tws);
    } else if (subcommand == "preprocess") {
        a.done();
        t->process_tws_and_compute_lbs(preprocess_tws, false); // DO NOT compute cost LBs
        t->dump();
    }
    delete t;
    return 0;
}
