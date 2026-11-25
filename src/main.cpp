// main.cpp -- plasma-0d driver: read config, integrate, write CSV.
//   usage: plasma0d <config.cfg> [output.csv]
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include "config.hpp"
#include "model.hpp"
#include "solver.hpp"

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: plasma0d <config.cfg> [output.csv]\n"; return 2; }
    const std::string out_path = argc > 2 ? argv[2] : "output.csv";
    try {
        Config cfg = load_config(argv[1]);
        Model model(cfg);

        Vec y = model.initial();
        const int n = model.n_state();

        Rosenbrock ros;
        ros.f = [&](const Vec& v, Vec& d) { model.rhs(v, d); };
        ros.atol.assign(n, cfg.atol_n);
        ros.atol[n - 1] = cfg.atol_Te;           // Te is in eV, not m^-3
        ros.rtol = cfg.rtol; ros.hmax = cfg.hmax; ros.h = cfg.h0;

        std::ofstream out(out_path);
        if (!out) throw std::runtime_error("cannot write " + out_path);
        out.precision(9);
        // Column units are in the header so a plot script or a person never guesses.
        out << "time_s";
        for (int i : model.evolved()) out << ',' << cfg.species[i].name << "_m3";
        out << ",n_e_m3,Te_eV\n";

        auto row = [&](double t) {
            out << t;
            for (int k = 0; k < n - 1; ++k) out << ',' << y[k];
            out << ',' << model.n_e(y) << ',' << y[n - 1] << '\n';
        };

        double t = 0;
        row(t);
        // Log-spaced outputs: the physics spans nanoseconds to milliseconds.
        const int N = std::max(cfg.n_outputs, 2);
        const double ratio = std::pow(cfg.t_end / cfg.t_first, 1.0 / (N - 1));
        double tt = cfg.t_first;
        for (int k = 0; k < N; ++k, tt *= ratio) {
            ros.advance(y, t, std::min(tt, cfg.t_end));
            row(t);
        }
        std::cout << "plasma0d: " << ros.steps << " steps, " << ros.rejected << " rejected, wrote "
                  << out_path << "\n";
        std::cout << "final: n_e = " << model.n_e(y) << " m^-3, Te = " << y[n - 1] << " eV\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
