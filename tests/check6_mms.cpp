// Check 6 -- manufactured solution with grid refinement.
//
// Method of manufactured solutions: choose smooth n(x,y) and phi(x,y), compute by calculus
// the source S that makes them an exact steady solution of
//     div( s mu n E - D grad n ) = S,   E = -grad phi,   mu = D / VT
// feed S to the code, and see whether it hands back n. Because the answer is known
// everywhere, the error is measured everywhere, and halving the cell size should cut it by
// the scheme's order (2 for a second-order method). This exercises drift, diffusion, the
// wall treatment and the sparse solve together, in 2D.
//
//   n = sin(kx x) sin(ky y)         phi = Phi0 sin(kx x) sin(ky y)    (both vanish on the walls)
//   div Gamma = -s (D/VT) (grad n . grad phi + n lap phi) - D lap n
//
// Tolerance (PROPOSED, Josue confirms): observed order at least 1.8 on 16 -> 32 and 32 -> 64.
#include <cmath>
#include <cstdio>
#include "../src/fluid2d.hpp"
#include "checkutil.hpp"

static const double kPi = 3.14159265358979323846;

static double run(int N) {
    Grid g; g.init(N, N, 1.0, 1.0);
    const double kx = kPi, ky = kPi, D = 1.0, VT = 1.0, s = 1.0, Phi0 = 2.0;
    Transport T{s, D, VT};
    Field phi(g.size()), src(g.size()), n(g.size(), 0.0), zero(g.size(), 0.0), exact(g.size());
    for (int j = 0; j < N; ++j) for (int i = 0; i < N; ++i) {
        const double x = g.x(i), y = g.y(j); const int k = g.idx(i, j);
        const double sx = std::sin(kx * x), cx = std::cos(kx * x), sy = std::sin(ky * y), cy = std::cos(ky * y);
        const double lap2 = kx * kx + ky * ky;
        exact[k] = sx * sy;
        phi[k] = Phi0 * sx * sy;
        const double gradDot = Phi0 * (kx * kx * cx * cx * sy * sy + ky * ky * sx * sx * cy * cy);
        const double lapN = -lap2 * exact[k], lapPhi = -lap2 * phi[k];
        src[k] = -s * (D / VT) * (gradDot + exact[k] * lapPhi) - D * lapN;
    }
    Stencil A;
    assembleTransport(g, T, phi, 1e12, zero, src, zero, A);   // steady: 1/dt ~ 0
    solveStencil(g, A, n, 1e-13);
    double err = 0; for (int k = 0; k < g.size(); ++k) err = std::max(err, std::fabs(n[k] - exact[k]));
    return err;
}

int main() {
    Checker ck;
    const int N[3] = {16, 32, 64}; double err[3];
    std::printf("check 6: manufactured drift-diffusion solution, max-norm error\n");
    for (int k = 0; k < 3; ++k) { err[k] = run(N[k]); std::printf("  %3dx%-3d  max err %.3e\n", N[k], N[k], err[k]); }
    const double o1 = observedOrder(err[0], err[1]), o2 = observedOrder(err[1], err[2]);
    std::printf("  observed order: %.2f, %.2f\n", o1, o2);
    ck.expect(o1 >= 1.8 && o2 >= 1.8, "observed order at least 1.8");
    return ck.finish("check6_mms");
}
