// Check 5 -- Poisson against a known field, and charge is conserved.
//
// (a) Manufactured field: pick phi = sin(kx x) sin(ky y) and a variable coefficient
//     kappa = 2 + cos(kx x), compute the matching right-hand side  f = -div(kappa grad phi)
//     analytically, solve, compare. Variable kappa matters: the semi-implicit scheme
//     solves exactly this kind of equation. Tolerance (PROPOSED): observed order >= 1.8.
// (b) Gauss's law in the discrete solve: the sum of f over all cells must equal the flux
//     through the walls (what the solved field says leaves the box). Tolerance: 1e-9 relative.
//     That is the statement "the field accounts for every charge".
// (c) Charge bookkeeping in the transport: in a closed box (mirror on every side, reactions
//     and beam off) the net charge and each species' total must not change. Tolerance: 1e-12
//     relative. Walls are the only place charge may leave.
#include <cmath>
#include <cstdio>
#include <random>
#include "../src/config.hpp"
#include "../src/fluid2d.hpp"
#include "checkutil.hpp"

static const double kPi = 3.14159265358979323846;

static double poissonError(int n, double& gaussGap) {
    Grid g; g.init(n, n, 1.0, 1.0);
    const double kx = kPi, ky = kPi;
    Field kappa(g.size()), f(g.size()), exact(g.size()), phi(g.size(), 0.0);
    for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        const double x = g.x(i), y = g.y(j), k = g.idx(i, j);
        kappa[k] = 2.0 + std::cos(kx * x);
        exact[k] = std::sin(kx * x) * std::sin(ky * y);
        f[k] = kappa[k] * (kx * kx + ky * ky) * exact[k] + kx * kx * std::sin(kx * x) * std::cos(kx * x) * std::sin(ky * y);
    }
    Stencil A; assemblePoisson(g, kappa, f, A);
    solveStencil(g, A, phi, 1e-13);
    double err = 0; for (int k = 0; k < g.size(); ++k) err = std::max(err, std::fabs(phi[k] - exact[k]));
    double sumf = 0, absf = 0; for (double v : f) { sumf += v * g.cellArea(); absf += std::fabs(v) * g.cellArea(); }
    gaussGap = std::fabs(sumf - wallFlux(g, kappa, phi)) / absf;
    return err;
}

int main(int argc, char** argv) {
    Checker ck;
    const int N[3] = {16, 32, 64}; double err[3], gap[3];
    std::printf("check 5a/b: Poisson with variable kappa vs manufactured field\n");
    for (int k = 0; k < 3; ++k) { err[k] = poissonError(N[k], gap[k]); std::printf("  %3dx%-3d  max err %.3e   Gauss gap %.2e\n", N[k], N[k], err[k], gap[k]); }
    const double o1 = observedOrder(err[0], err[1]), o2 = observedOrder(err[1], err[2]);
    std::printf("  observed order: %.2f, %.2f\n", o1, o2);
    ck.expect(o1 >= 1.8 && o2 >= 1.8, "observed order at least 1.8");
    ck.expect(gap[0] < 1e-9 && gap[1] < 1e-9 && gap[2] < 1e-9, "sum of charge equals wall flux (Gauss) to 1e-9");

    // (c) closed-box charge bookkeeping
    Config cfg = load_config(argc > 1 ? argv[1] : "config/argon2d_standin.cfg");
    cfg.extra["grid"]["nx"] = 20; cfg.extra["grid"]["ny"] = 20;
    for (const char* s : {"bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax"}) cfg.extra["grid"][s] = 1;
    Fluid2D sim(cfg);
    sim.reactionsOn = false; sim.beamOn = false;
    std::mt19937 rng(12345); std::uniform_real_distribution<double> u(0.2, 1.0);
    for (size_t s = 0; s < sim.n.size(); ++s) for (double& v : sim.n[s]) v = 1e15 * u(rng);   // lumpy, unequal e and ion
    std::vector<double> tot0(sim.n.size()); double q0 = sim.totalCharge();
    for (size_t s = 0; s < sim.n.size(); ++s) { double t = 0; for (double v : sim.n[s]) t += v; tot0[s] = t; }
    for (int k = 0; k < 200; ++k) sim.step(1e-6);
    double worst = 0;
    for (size_t s = 0; s < sim.n.size(); ++s) { double t = 0; for (double v : sim.n[s]) t += v; worst = std::max(worst, std::fabs(t - tot0[s]) / tot0[s]); }
    const double qrel = std::fabs(sim.totalCharge() - q0) / (phys2::e * 1e15 * sim.g.cellArea() * sim.g.size());
    std::printf("check 5c: closed box, 200 steps: worst relative change of a species total %.2e, of net charge %.2e\n", worst, qrel);
    ck.expect(worst < 1e-12, "every species' total is conserved to 1e-12");
    ck.expect(qrel < 1e-12, "net charge is conserved to 1e-12");
    return ck.finish("check5_poisson");
}
