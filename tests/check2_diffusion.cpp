// Check 2 -- pure diffusion against the closed-form answer.
//
// A neutral species (no drift) on a square with absorbing walls, started as the lowest
// eigenmode  sin(pi x/Lx) sin(pi y/Ly).  The exact solution just decays:
//     n(t) = n(0) exp(-D (pi^2/Lx^2 + pi^2/Ly^2) t)
// We compare the numerical amplitude after one e-folding time on three grids.
//
// Tolerances (PROPOSED, Josue confirms): finest-grid amplitude error below 0.1 %, and the
// observed order of accuracy at least 1.8 (the scheme is second order in space; the time
// step is small enough that its first-order error stays below the spatial error).
#include <cmath>
#include <cstdio>
#include "../src/fluid2d.hpp"
#include "checkutil.hpp"

static double run(int n) {
    Grid g; g.init(n, n, 1.0, 1.0);
    const double pi = 3.14159265358979323846, D = 1.0;
    Transport T{0.0, D, 1.0};
    Field phi(g.size(), 0.0), loss(g.size(), 0.0), src(g.size(), 0.0), c(g.size());
    double amp0 = 0;
    for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        c[g.idx(i, j)] = std::sin(pi * g.x(i)) * std::sin(pi * g.y(j));
    }
    // Amplitude measured as the projection onto the mode, so it does not depend on one point.
    auto amp = [&](const Field& f) { double a = 0, b = 0; for (int k = 0; k < g.size(); ++k) { int i = k % n, j = k / n;
        double m = std::sin(pi * g.x(i)) * std::sin(pi * g.y(j)); a += f[k] * m; b += m * m; } return a / b; };
    amp0 = amp(c);
    const double lambda = D * 2 * pi * pi, tEnd = 1.0 / lambda;
    const int steps = 40000; const double dt = tEnd / steps;   // small enough that time error stays below the spatial error
    for (int s = 0; s < steps; ++s) {
        Field old = c; Stencil A;
        assembleTransport(g, T, phi, dt, old, src, loss, A);
        solveStencil(g, A, c);
    }
    return amp(c) / amp0;   // exact: exp(-1)
}

int main() {
    Checker ck;
    const double exact = std::exp(-1.0);
    double err[3]; const int N[3] = {16, 32, 64};
    std::printf("check 2: pure diffusion, amplitude after one e-folding (exact %.8f)\n", exact);
    for (int k = 0; k < 3; ++k) { const double a = run(N[k]); err[k] = std::fabs(a - exact) / exact; std::printf("  %3dx%-3d  amp %.8f  rel err %.3e\n", N[k], N[k], a, err[k]); }
    const double o1 = observedOrder(err[0], err[1]), o2 = observedOrder(err[1], err[2]);
    std::printf("  observed order: %.2f, %.2f\n", o1, o2);
    ck.expect(err[2] < 1e-3, "finest-grid error below 0.1 %");
    ck.expect(o2 >= 1.8, "observed order at least 1.8 on the fine pair");
    return ck.finish("check2_diffusion");
}
