// solver.hpp -- adaptive Rosenbrock integrator (ROS2, Verwer et al. 1999).
//
// Why Rosenbrock: plasma chemistry is stiff (electron heating and ionization
// happen in nanoseconds, wall losses in milliseconds). An explicit method would
// be forced to the fastest timescale for the whole run. Rosenbrock methods are
// implicit-quality but need only a linear solve per stage, no Newton iteration,
// and the step-size estimate comes for free from the embedded first-order result.
#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

using Vec = std::vector<double>;
using Rhs = std::function<void(const Vec& y, Vec& dydt)>;

namespace solver_detail {

// Dense LU with partial pivoting. N is tiny (a handful of species), so clarity
// beats speed here.
struct LU {
    int n = 0;
    std::vector<double> a;
    std::vector<int> piv;
    void factor(std::vector<double> m, int N) {
        n = N; a = std::move(m); piv.assign(n, 0);
        for (int k = 0; k < n; ++k) {
            int p = k;
            for (int i = k + 1; i < n; ++i) if (std::fabs(a[i * n + k]) > std::fabs(a[p * n + k])) p = i;
            if (a[p * n + k] == 0.0) throw std::runtime_error("singular matrix in Rosenbrock step");
            piv[k] = p;
            if (p != k) for (int j = 0; j < n; ++j) std::swap(a[k * n + j], a[p * n + j]);
            for (int i = k + 1; i < n; ++i) {
                a[i * n + k] /= a[k * n + k];
                for (int j = k + 1; j < n; ++j) a[i * n + j] -= a[i * n + k] * a[k * n + j];
            }
        }
    }
    void solve(Vec& b) const {
        for (int k = 0; k < n; ++k) {
            std::swap(b[k], b[piv[k]]);
            for (int i = k + 1; i < n; ++i) b[i] -= a[i * n + k] * b[k];
        }
        for (int i = n - 1; i >= 0; --i) {
            for (int j = i + 1; j < n; ++j) b[i] -= a[i * n + j] * b[j];
            b[i] /= a[i * n + i];
        }
    }
};

} // namespace solver_detail

struct Rosenbrock {
    Rhs f;
    Vec atol;           // per-component absolute tolerance (densities and Te differ by 15 decades)
    double rtol = 1e-5, hmax = 1e-3;
    double h = 1e-12;   // carried between calls so each output interval starts at a good step
    long steps = 0, rejected = 0;

    // Integrate y from t to t_target (exactly), updating y in place.
    void advance(Vec& y, double& t, double t_target) {
        const int n = (int)y.size();
        const double g = 1.0 + 1.0 / std::sqrt(2.0);   // gamma of ROS2: L-stable
        Vec f0(n), f1(n), k1(n), k2(n), yt(n), ynew(n), err(n);
        while (t < t_target) {
            const double h_try = std::min(h, hmax);
            double hh = std::min(h_try, t_target - t);
            const bool truncated = hh < h_try;
            f(y, f0);

            // Finite-difference Jacobian. The perturbation is relative to the
            // component's own scale so densities (1e16) and Te (a few eV) both
            // get a sensible increment.
            std::vector<double> J(n * n);
            for (int j = 0; j < n; ++j) {
                double sc = std::max(std::fabs(y[j]), atol[j]);
                double d = 1e-7 * sc;
                Vec yp = y; yp[j] += d;
                Vec fp(n); f(yp, fp);
                for (int i = 0; i < n; ++i) J[i * n + j] = (fp[i] - f0[i]) / d;
            }

            bool accepted = false;
            while (!accepted) {
                std::vector<double> M(n * n);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) M[i * n + j] = (i == j ? 1.0 : 0.0) - g * hh * J[i * n + j];
                solver_detail::LU lu; lu.factor(M, n);

                k1 = f0; lu.solve(k1);
                for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * k1[i];
                f(yt, f1);
                for (int i = 0; i < n; ++i) k2[i] = f1[i] - 2.0 * k1[i];
                lu.solve(k2);
                for (int i = 0; i < n; ++i) {
                    ynew[i] = y[i] + hh * (1.5 * k1[i] + 0.5 * k2[i]);
                    err[i] = 0.5 * hh * (k1[i] + k2[i]);   // 2nd-order minus embedded 1st-order
                }
                double en = 0;
                for (int i = 0; i < n; ++i) {
                    double sc = atol[i] + rtol * std::max(std::fabs(y[i]), std::fabs(ynew[i]));
                    en += (err[i] / sc) * (err[i] / sc);
                }
                en = std::sqrt(en / n);
                bool finite = true;
                for (double v : ynew) if (!std::isfinite(v)) finite = false;
                // Order-1 error estimate: optimal step scales as err^(-1/2).
                double fac = (finite && en > 0) ? 0.9 / std::sqrt(en) : (finite ? 5.0 : 0.1);
                fac = std::min(5.0, std::max(0.1, fac));
                if (finite && en <= 1.0) {
                    y = ynew; t += hh; ++steps; accepted = true;
                    h = std::max(hh * fac, 1e-30);
                    // A step shortened only to land on t_target says nothing about
                    // the next step's size, so do not let it shrink h.
                    if (truncated) h = std::max(h, h_try);
                } else {
                    ++rejected; hh *= fac;
                    if (hh < 1e-30) throw std::runtime_error("step size underflow: the model is not integrable from here");
                }
            }
        }
    }
};
