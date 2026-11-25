// Check 4 -- Scharfetter-Gummel reproduces the exact steady state.
//
// Steady drift-diffusion in a uniform field between two walls with fixed densities has an
// exact solution (constant flux):
//     n(x) = nL + (nR - nL) * (e^{v x/D} - 1) / (e^{v L/D} - 1)    (v/D = field / thermal voltage)
// The Scharfetter-Gummel flux is built from exactly this solution, so on a 1D problem it
// should reproduce it to the linear-solver tolerance, even when drift is so strong that the
// profile is a thin boundary layer a plain central-difference scheme cannot represent.
// Here v L / D = 20 and the layer is about two cells wide.
//
// Tolerance (PROPOSED, Josue confirms): relative error below 1e-8 at every cell.
#include <cmath>
#include <cstdio>
#include "../src/fluid2d.hpp"
#include "checkutil.hpp"

int main() {
    Checker ck;
    const int nx = 40, ny = 3;
    const double L = 1.0, D = 1.0, VT = 1.0, E0 = 20.0;   // Peclet number E0 L / VT = 20
    Grid g; g.init(nx, ny, L, 1.0);
    g.bc[YMIN] = g.bc[YMAX] = SYM;
    g.nWall[XMIN] = 1.0; g.nWall[XMAX] = 0.0;
    g.phiWall[XMIN] = 0.0; g.phiWall[XMAX] = -E0 * L;
    Transport T{1.0, D, VT};
    Field phi(g.size()), n(g.size(), 0.0), loss(g.size(), 0.0), src(g.size(), 0.0);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) phi[g.idx(i, j)] = -E0 * g.x(i);
    // dt huge and the old density zero: the 1/dt terms vanish, leaving the steady equation.
    Stencil A;
    assembleTransport(g, T, phi, 1e12, n, src, loss, A);
    solveStencil(g, A, n, 1e-13);

    const double Pe = E0 * L / VT, nL = 1.0, nR = 0.0;
    double worst = 0;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const double x = g.x(i);
        const double ex = nL + (nR - nL) * std::expm1(Pe * x / L) / std::expm1(Pe);
        worst = std::max(worst, std::fabs(n[g.idx(i, j)] - ex) / std::max(std::fabs(ex), 1e-30));
    }
    std::printf("check 4: strong drift (Pe = %g), steady state vs exact; worst relative error %.3e\n", Pe, worst);
    ck.expect(worst < 1e-8, "every cell within 1e-8 of the exact profile");
    return ck.finish("check4_sgsteady");
}
