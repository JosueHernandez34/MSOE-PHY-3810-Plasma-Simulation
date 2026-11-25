// Check 1 -- uniform conditions reduce to the 0D result.
//
// A closed box (mirror on every side), no beam, uniform start: nothing can make the plasma
// non-uniform, so every cell must follow the same ODE system the 0D model solves. The
// reference is that ODE system integrated with the v0.1.0 Rosenbrock solver (solver.hpp,
// reused unchanged), written independently here from the same config reactions.
//
// What this does NOT test: the 0D model's own Te equation and wall terms, because v0.2.0
// holds Te constant and has no Bohm wall. Te is set to the v0.1.0 steady-state value, so the
// reaction rates are the ones the 0D stand-in case actually sits at.
//
// The 2D step is first order in time, so a single dt cannot tell "close" from "right". The
// check runs three time steps, each 10x smaller, and requires the difference to the ODE to
// shrink about 10x each time (first-order convergence) and to end below 0.1 %.
//
// Tolerances (PROPOSED, Josue confirms): finest-dt difference below 0.1 %; each 10x cut in dt
// reduces the difference by at least 8x; cell-to-cell spread below 1e-9 relative.
// (A first draft assumed dt = 1 us would already reach 1 %; measured, it gives about 5 %.
// That was a wrong estimate in the plan, not a code fault; the scan below is the evidence.)
#include <cmath>
#include <cstdio>
#include <vector>
#include "../src/config.hpp"
#include "../src/fluid2d.hpp"
#include "../src/solver.hpp"
#include "checkutil.hpp"

int main(int argc, char** argv) {
    Checker ck;
    Config cfg = load_config(argc > 1 ? argv[1] : "config/argon2d_standin.cfg");
    cfg.pressure = 1.333;   // v0.1.0 gas density, so rates match the 0D case
    cfg.Te0 = 1.86;         // v0.1.0 steady-state Te
    cfg.extra["grid"]["nx"] = 6; cfg.extra["grid"]["ny"] = 6;
    for (const char* s : {"bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax"}) cfg.extra["grid"][s] = 1;
    const double tEnd = 5e-3;

    // Reference ODE: state = [n_e, evolved species in config order].
    std::vector<int> ev;
    for (size_t i = 0; i < cfg.species.size(); ++i) if (!cfg.species[i].fixed) ev.push_back((int)i);
    const double nGas = cfg.pressure / (phys2::kB * cfg.Tgas);
    int nfixed = 0; for (auto& s : cfg.species) if (s.fixed) ++nfixed;
    auto slot = [&](int sp) { if (sp == -1) return 0; for (size_t k = 0; k < ev.size(); ++k) if (ev[k] == sp) return (int)k + 1; return -1; };
    Rosenbrock ros;
    ros.f = [&](const Vec& y, Vec& d) {
        d.assign(y.size(), 0.0);
        for (const Reaction& r : cfg.reactions) {
            double rate = r.A * std::pow(cfg.Te0, r.B) * std::exp(-r.E / cfg.Te0);
            for (const Term& t : r.reactants) {
                const int s = slot(t.species);
                rate *= std::pow(s >= 0 ? std::max(y[s], 0.0) : nGas / nfixed, t.count);
            }
            for (const Term& t : r.reactants) { const int s = slot(t.species); if (s >= 0) d[s] -= t.count * rate; }
            for (const Term& t : r.products)  { const int s = slot(t.species); if (s >= 0) d[s] += t.count * rate; }
        }
    };
    Vec y(ev.size() + 1);
    double ne0 = 0;
    for (size_t k = 0; k < ev.size(); ++k) { y[k + 1] = cfg.species[ev[k]].n0; ne0 += cfg.species[ev[k]].charge * cfg.species[ev[k]].n0; }
    y[0] = ne0;
    ros.atol.assign(y.size(), 1e5); ros.rtol = 1e-9; ros.hmax = 1e-4; ros.h = 1e-9;
    double t = 0;
    ros.advance(y, t, tEnd);

    std::printf("check 1: uniform 2D vs 0D reaction ODE at t = %g s (Te = %g eV held)\n", tEnd, cfg.Te0);
    const double dts[3] = {1e-6, 1e-7, 1e-8};
    double worst[3] = {0, 0, 0}, spreadMax = 0;
    for (int d = 0; d < 3; ++d) {
        Fluid2D sim(cfg);
        sim.beamOn = false;
        const long steps = (long)std::llround(tEnd / dts[d]);
        for (long k = 0; k < steps; ++k) sim.step(dts[d]);
        for (size_t s = 0; s < sim.n.size(); ++s) {
            const double got = sim.mean((int)s);
            worst[d] = std::max(worst[d], std::fabs(got - y[s]) / std::fabs(y[s]));
            for (double v : sim.n[s]) spreadMax = std::max(spreadMax, std::fabs(v - got) / got);
        }
        std::printf("  dt = %.0e s: worst species difference to the ODE %.3e\n", dts[d], worst[d]);
    }
    ck.expect(worst[2] < 1e-3, "difference below 0.1 % at the smallest dt");
    ck.expect(worst[0] / worst[1] >= 8 && worst[1] / worst[2] >= 8, "difference shrinks at least 8x per 10x smaller dt (first order)");
    ck.expect(spreadMax < 1e-9, "all cells identical (uniform stays uniform)");
    return ck.finish("check1_reduction");
}
