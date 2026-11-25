// main2d.cpp -- plasma-2d driver: read config, step the 2D fluid model, write CSVs.
//   usage: plasma2d <config.cfg> [history.csv] [fields.csv]
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include "config.hpp"
#include "fluid2d.hpp"

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: plasma2d <config.cfg> [history.csv] [fields.csv]\n"; return 2; }
    const std::string histPath = argc > 2 ? argv[2] : "history2d.csv";
    const std::string fieldPath = argc > 3 ? argv[3] : "fields2d.csv";
    try {
        Config cfg = load_config(argv[1]);
        Fluid2D sim(cfg);
        double dt = 1e-6;
        auto sv = cfg.extra.find("solver2d");
        if (sv != cfg.extra.end() && sv->second.count("dt_s")) dt = sv->second["dt_s"];

        std::ofstream hist(histPath);
        if (!hist) throw std::runtime_error("cannot write " + histPath);
        hist.precision(9);
        // Volume averages over the slab; units in the header so nothing is guessed.
        hist << "time_s,n_e_mean_m3";
        for (int idx : sim.speciesIdx) hist << ',' << cfg.species[idx].name << "_mean_m3";
        hist << ",charge_C_per_m\n";
        auto row = [&]() {
            hist << sim.time << ',' << sim.mean(0);
            for (size_t s = 1; s < sim.n.size(); ++s) hist << ',' << sim.mean((int)s);
            hist << ',' << sim.totalCharge() << '\n';
        };

        const int nOut = std::max(cfg.n_outputs, 2);
        const long nSteps = (long)std::llround(cfg.t_end / dt);
        const long every = std::max(1L, nSteps / nOut);
        row();
        for (long k = 1; k <= nSteps; ++k) {
            sim.step(dt);
            if (k % every == 0 || k == nSteps) row();
        }

        std::ofstream fl(fieldPath);
        if (!fl) throw std::runtime_error("cannot write " + fieldPath);
        fl.precision(9);
        fl << "x_m,y_m,n_e_m3";
        for (int idx : sim.speciesIdx) fl << ',' << cfg.species[idx].name << "_m3";
        fl << ",phi_V\n";
        for (int j = 0; j < sim.g.ny; ++j)
            for (int i = 0; i < sim.g.nx; ++i) {
                const int k = sim.g.idx(i, j);
                fl << sim.g.x(i) << ',' << sim.g.y(j);
                for (size_t s = 0; s < sim.n.size(); ++s) fl << ',' << sim.n[s][k];
                fl << ',' << sim.phi[k] << '\n';
            }
        std::cout << "plasma2d: " << nSteps << " steps of " << dt << " s on " << sim.g.nx << "x" << sim.g.ny
                  << " cells; wrote " << histPath << " and " << fieldPath << "\n";
        std::cout << "final means: n_e = " << sim.mean(0) << " m^-3";
        for (size_t s = 1; s < sim.n.size(); ++s) std::cout << ", " << cfg.species[sim.speciesIdx[s - 1]].name << " = " << sim.mean((int)s);
        std::cout << "\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
