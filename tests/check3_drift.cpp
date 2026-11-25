// Check 3 -- drift (with diffusion) against the closed-form answer.
//
// A Gaussian bump of ions in a uniform field. The exact solution of
//     dn/dt + d/dx( v n - D dn/dx ) = 0
// is a Gaussian that moves at v and widens:
//     n(x,t) = A (s0/s) exp( -(x - x0 - v t)^2 / (2 s^2) ),   s^2 = s0^2 + 2 D t
// "Pure drift" is the limit D -> 0; D is kept small but non-zero because diffusion is
// what keeps the discrete solution smooth. The Scharfetter-Gummel flux has to get BOTH the
// speed and the widening right, so one comparison checks both.
//
// The field is imposed (a straight-line potential), not solved: this isolates the transport.
//
// Tolerances (PROPOSED, Josue confirms): relative L2 error below 0.5 % on the finest grid,
// and the error must fall as the grid is refined.
#include <cmath>
#include <cstdio>
#include "../src/fluid2d.hpp"
#include "checkutil.hpp"

static double run(int nx) {
    const int ny = 4;
    const double Lx = 4.0, Ly = 1.0, D = 0.02, VT = 1.0, E0 = 8.0;   // v = (D/VT) E0 = 0.16
    const double v = D / VT * E0, x0 = 1.2, s0 = 0.15, tEnd = 2.0;
    Grid g; g.init(nx, ny, Lx, Ly);
    g.bc[YMIN] = g.bc[YMAX] = SYM;                      // a 1D problem in a 2D code
    g.phiWall[XMIN] = 0.0; g.phiWall[XMAX] = -E0 * Lx;  // consistent with phi(x) = -E0 x
    Transport T{1.0, D, VT};
    Field phi(g.size()), n(g.size()), loss(g.size(), 0.0), src(g.size(), 0.0);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        phi[g.idx(i, j)] = -E0 * g.x(i);
        const double X = g.x(i) - x0;
        n[g.idx(i, j)] = std::exp(-X * X / (2 * s0 * s0));
    }
    const int steps = 8000; const double dt = tEnd / steps;
    for (int s = 0; s < steps; ++s) {
        Field old = n; Stencil A;
        assembleTransport(g, T, phi, dt, old, src, loss, A);
        solveStencil(g, A, n);
    }
    const double s2 = s0 * s0 + 2 * D * tEnd, sg = std::sqrt(s2);
    double num = 0, den = 0;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const double X = g.x(i) - x0 - v * tEnd;
        const double ex = (s0 / sg) * std::exp(-X * X / (2 * s2));
        const double d = n[g.idx(i, j)] - ex;
        num += d * d; den += ex * ex;
    }
    return std::sqrt(num / den);
}

int main() {
    Checker ck;
    const int N[3] = {32, 64, 128}; double err[3];
    std::printf("check 3: drifting, widening Gaussian vs closed form (v = 0.16, D = 0.02, cell Peclet 1, 0.5, 0.25)\n");
    for (int k = 0; k < 3; ++k) { err[k] = run(N[k]); std::printf("  %4d cells  rel L2 err %.3e\n", N[k], err[k]); }
    std::printf("  observed order: %.2f, %.2f\n", observedOrder(err[0], err[1]), observedOrder(err[1], err[2]));
    ck.expect(err[2] < 5e-3, "finest-grid relative L2 error below 0.5 %");
    ck.expect(err[0] > err[1] && err[1] > err[2], "error falls with every refinement");
    return ck.finish("check3_drift");
}
