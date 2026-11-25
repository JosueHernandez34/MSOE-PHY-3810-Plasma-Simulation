// poisson2d.hpp -- the sparse linear solve (no libraries) and the Poisson equation.
//
// Both the Poisson equation and the implicit transport step reduce to the same
// shape: one equation per cell that links the cell to its four neighbors (a
// "5-point stencil"). So one small solver, BiCGStab, serves both. BiCGStab
// handles the non-symmetric transport matrix; for the symmetric Poisson matrix
// it simply works too, which is worth more than saving a few iterations.
//
// Poisson here has a variable coefficient:
//     -div( kappa grad phi ) = f
// because the semi-implicit scheme (fluid2d.hpp) folds the plasma's own response
// into kappa. With kappa = eps0 it is the plain Poisson equation.
#pragma once
#include <cmath>
#include <stdexcept>
#include "grid2d.hpp"

// Row k reads:  aC x_k + aE x_{east} + aW x_{west} + aN x_{north} + aS x_{south} = b_k
// Coefficients that point off the grid are simply never read.
struct Stencil {
    Field aC, aE, aW, aN, aS, b;
    void reset(int n) {
        aC.assign(n, 0.0); aE.assign(n, 0.0); aW.assign(n, 0.0);
        aN.assign(n, 0.0); aS.assign(n, 0.0); b.assign(n, 0.0);
    }
};

inline void applyStencil(const Grid& g, const Stencil& A, const Field& x, Field& y) {
    y.resize(x.size());
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const int k = g.idx(i, j);
            double s = A.aC[k] * x[k];
            if (i + 1 < g.nx) s += A.aE[k] * x[k + 1];
            if (i > 0)        s += A.aW[k] * x[k - 1];
            if (j + 1 < g.ny) s += A.aN[k] * x[k + g.nx];
            if (j > 0)        s += A.aS[k] * x[k - g.nx];
            y[k] = s;
        }
}

namespace lin_detail {
inline double dot(const Field& a, const Field& b) { double s = 0; for (size_t k = 0; k < a.size(); ++k) s += a[k] * b[k]; return s; }
inline double norm(const Field& a) { return std::sqrt(dot(a, a)); }
}

// Preconditioned BiCGStab (Jacobi preconditioner = divide by the diagonal).
// x holds the starting guess on entry and the solution on exit.
// Returns the iteration count; throws if it does not converge, because a silently
// wrong field would look like physics.
inline int solveStencil(const Grid& g, const Stencil& A, Field& x, double relTol = 1e-14, int maxIter = 5000) {
    using namespace lin_detail;
    const size_t n = A.b.size();
    if (x.size() != n) x.assign(n, 0.0);
    Field r(n), rhat, p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n), Ax;
    applyStencil(g, A, x, Ax);
    for (size_t k = 0; k < n; ++k) r[k] = A.b[k] - Ax[k];
    const double bnorm = std::max(norm(A.b), 1e-300);
    if (norm(r) <= relTol * bnorm) return 0;
    rhat = r;
    double rho = 1, alpha = 1, omega = 1;
    for (int it = 1; it <= maxIter; ++it) {
        const double rhoNew = dot(rhat, r);
        if (std::fabs(rhoNew) < 1e-300) break;
        const double beta = (rhoNew / rho) * (alpha / omega);
        for (size_t k = 0; k < n; ++k) p[k] = r[k] + beta * (p[k] - omega * v[k]);
        for (size_t k = 0; k < n; ++k) ph[k] = p[k] / A.aC[k];
        applyStencil(g, A, ph, v);
        alpha = rhoNew / dot(rhat, v);
        for (size_t k = 0; k < n; ++k) s[k] = r[k] - alpha * v[k];
        if (norm(s) <= relTol * bnorm) {
            for (size_t k = 0; k < n; ++k) x[k] += alpha * ph[k];
            return it;
        }
        for (size_t k = 0; k < n; ++k) sh[k] = s[k] / A.aC[k];
        applyStencil(g, A, sh, t);
        omega = dot(t, s) / dot(t, t);
        for (size_t k = 0; k < n; ++k) { x[k] += alpha * ph[k] + omega * sh[k]; r[k] = s[k] - omega * t[k]; }
        if (norm(r) <= relTol * bnorm) return it;
        rho = rhoNew;
    }
    throw std::runtime_error("linear solve did not converge");
}

// Assemble  -div(kappa grad phi) = f  on the grid.
// Interior face: flux = kappa_face (phi_C - phi_N) / h, kappa_face = mean of the two cells.
// WALL face: phi is fixed on the face, which sits h/2 away, hence the factor 2.
// SYM face: no flux, nothing to add.
inline void assemblePoisson(const Grid& g, const Field& kappa, const Field& f, Stencil& A) {
    A.reset(g.size());
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const int k = g.idx(i, j);
            A.b[k] = f[k];
            for (int fc = 0; fc < 4; ++fc) {
                const Face F = faceOf(g, i, j, fc);
                if (F.nb >= 0) {
                    const double c = 0.5 * (kappa[k] + kappa[F.nb]) / (F.h * F.h);
                    A.aC[k] += c;
                    (fc == 0 ? A.aE : fc == 1 ? A.aW : fc == 2 ? A.aN : A.aS)[k] -= c;
                } else if (g.bc[F.side] == WALL) {
                    const double c = 2.0 * kappa[k] / (F.h * F.h);
                    A.aC[k] += c;
                    A.b[k]  += c * g.phiWall[F.side];
                }
            }
        }
}

// Net outward flux of -kappa*grad(phi) summed over all WALL faces, per unit depth.
// By the divergence theorem this must equal  sum of (f * cell area): the test for
// "charge is conserved" in the Poisson solve uses exactly that identity.
inline double wallFlux(const Grid& g, const Field& kappa, const Field& phi) {
    double tot = 0;
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i)
            for (int fc = 0; fc < 4; ++fc) {
                const Face F = faceOf(g, i, j, fc);
                if (F.nb < 0 && g.bc[F.side] == WALL) {
                    const double len = (F.side == XMIN || F.side == XMAX) ? g.dy : g.dx;
                    const int k = g.idx(i, j);
                    tot += kappa[k] * (phi[k] - g.phiWall[F.side]) / (0.5 * F.h) * len;
                }
            }
    return tot;
}
