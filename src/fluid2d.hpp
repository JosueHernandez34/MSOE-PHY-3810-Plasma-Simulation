// fluid2d.hpp -- the 2D drift-diffusion + Poisson model (v0.2.0).
//
// Per species k:   d n_k/dt + div(Gamma_k) = S_k
//   Gamma_k = s_k mu_k n_k E - D_k grad n_k        (s_k = sign of the charge, 0 for neutrals)
//   -div(eps0 grad phi) = rho,   E = -grad phi,   D_k = mu_k * T_k  (Einstein relation, T_k in volts)
//
// Why each numerical choice (Learning Guide, Stage 1):
//  * Scharfetter-Gummel flux: handles drift and diffusion in one formula that falls
//    back to central differences when diffusion wins and to upwinding when drift wins.
//  * Implicit in the density (backward Euler): no diffusion time-step limit.
//  * Semi-implicit Poisson: the field is solved FIRST from a Poisson equation whose
//    coefficient already includes how the charges will respond during the step. That
//    removes the dielectric-relaxation limit (dt < eps0 / conductivity), which would
//    otherwise be ~1e-10 s here. Cost: first-order accuracy in time.
//  * Electron temperature is a CONSTANT in this version (config Te0_eV). The electron
//    energy equation is an open decision (session page); nothing here depends on it.
#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "config.hpp"
#include "grid2d.hpp"
#include "poisson2d.hpp"

namespace phys2 {
constexpr double kB   = 1.380649e-23;
constexpr double e    = 1.602176634e-19;
constexpr double eps0 = 8.8541878128e-12;
}

// Bernoulli function B(x) = x / (e^x - 1), with B(0) = 1. expm1 keeps it accurate near 0.
inline double bernoulli(double x) {
    if (std::fabs(x) < 1e-9) return 1.0 - 0.5 * x;
    return x / std::expm1(x);
}

// What a species does in space.
struct Transport {
    double s  = 0;   // charge number z (0 = neutral: diffusion only)
    double D  = 0;   // diffusion coefficient, m^2/s
    double VT = 1;   // thermal voltage T/q in volts; D = mu * VT (Einstein)
};

// Assemble the implicit transport equation for ONE species on a given potential:
//   (n_new - n_old)/dt + div(Gamma_SG(n_new)) + loss * n_new = src
// loss >= 0 (1/s) and src >= 0 (m^-3 s^-1) keep the matrix an M-matrix, so densities
// cannot go negative.
inline void assembleTransport(const Grid& g, const Transport& T, const Field& phi, double dt,
                              const Field& nOld, const Field& src, const Field& loss, Stencil& A) {
    A.reset(g.size());
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const int k = g.idx(i, j);
            A.aC[k] = 1.0 / dt + loss[k];
            A.b[k]  = nOld[k] / dt + src[k];
            for (int fc = 0; fc < 4; ++fc) {
                const Face F = faceOf(g, i, j, fc);
                if (F.nb >= 0) {
                    // a = s (phi_C - phi_N) / VT : voltage drop across the face in thermal units.
                    const double a = T.s * (phi[k] - phi[F.nb]) / T.VT;
                    const double c = T.D / (F.h * F.h);
                    A.aC[k] += c * bernoulli(-a);
                    (fc == 0 ? A.aE : fc == 1 ? A.aW : fc == 2 ? A.aN : A.aS)[k] -= c * bernoulli(a);
                } else if (g.bc[F.side] == WALL) {
                    // Wall sits h/2 away with a fixed density: same flux law over half a cell.
                    const double a = T.s * (phi[k] - g.phiWall[F.side]) / T.VT;
                    const double c = 2.0 * T.D / (F.h * F.h);
                    A.aC[k] += c * bernoulli(-a);
                    A.b[k]  += c * bernoulli(a) * g.nWall[F.side];
                }
            }
        }
}

// div( D grad n ) at each cell, same face rules as above (used by the Poisson right-hand side).
inline void diffusionDivergence(const Grid& g, double D, const Field& n, Field& out) {
    out.assign(g.size(), 0.0);
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const int k = g.idx(i, j);
            double s = 0;
            for (int fc = 0; fc < 4; ++fc) {
                const Face F = faceOf(g, i, j, fc);
                if (F.nb >= 0) s += D * (n[F.nb] - n[k]) / (F.h * F.h);
                else if (g.bc[F.side] == WALL) s += 2.0 * D * (g.nWall[F.side] - n[k]) / (F.h * F.h);
            }
            out[k] = s;
        }
}

class Fluid2D {
public:
    Grid g;
    Field phi;
    std::vector<Field> n;            // n[0] = electrons, n[1+k] = k-th evolved species of the config
    std::vector<int>   speciesIdx;   // config species index of n[1+k]
    std::vector<Transport> tr;       // same ordering as n
    double Te = 2.0;                 // eV, held constant
    double time = 0;
    bool reactionsOn = true, beamOn = true, poissonOn = true;
    Field beam;                      // prescribed pair-production profile, m^-3 s^-1
    double beamExciteRatio = 0;
    int ionIdx = -1, metaIdx = -1;   // positions in n[] that the beam feeds

    explicit Fluid2D(const Config& c) : c_(c) {
        auto ex = [&](const char* sec, const char* key, double d) {
            auto s = c.extra.find(sec);
            if (s == c.extra.end()) return d;
            auto v = s->second.find(key);
            return v == s->second.end() ? d : v->second;
        };
        g.init((int)ex("grid", "nx", 40), (int)ex("grid", "ny", 40), ex("grid", "Lx_m", 0.02), ex("grid", "Ly_m", 0.02));
        const char* sides[4] = {"bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax"};
        for (int s = 0; s < 4; ++s) g.bc[s] = ex("grid", sides[s], 0) != 0 ? SYM : WALL;

        Te = c.Te0;
        const double mu_e = ex("transport", "mu_e_m2Vs", 40), mu_i = ex("transport", "mu_i_m2Vs", 0.16);
        const double Tgas_V = phys2::kB * c.Tgas / phys2::e;
        ngas_ = c.pressure / (phys2::kB * c.Tgas);

        int nfixed = 0;
        for (const Species& s : c.species) if (s.fixed) ++nfixed;
        if (nfixed == 0) throw std::runtime_error("need at least one fixed (background gas) species");
        nfix_ = ngas_ / nfixed;

        // Electrons first.
        tr.push_back({-1.0, mu_e * Te, Te});
        double ne0 = 0;
        for (size_t i = 0; i < c.species.size(); ++i) {
            const Species& s = c.species[i];
            if (s.fixed) continue;
            Transport t;
            t.s = s.charge;
            if (s.charge > 0) { t.D = mu_i * Tgas_V; t.VT = Tgas_V; if (ionIdx < 0) ionIdx = (int)tr.size(); }
            else {
                t.D = 0;
                for (const Wall& w : c.walls) if (w.species == (int)i && !w.bohm) t.D = w.Dp / c.pressure;
                if (metaIdx < 0) metaIdx = (int)tr.size();
            }
            speciesIdx.push_back((int)i);
            tr.push_back(t);
            ne0 += s.charge * s.n0;   // electrons start neutralizing the initial ions
        }
        n.assign(tr.size(), Field(g.size(), 0.0));
        std::fill(n[0].begin(), n[0].end(), std::max(ne0, 1.0));
        for (size_t k = 0; k < speciesIdx.size(); ++k)
            std::fill(n[k + 1].begin(), n[k + 1].end(), c.species[speciesIdx[k]].n0);
        phi.assign(g.size(), 0.0);

        // Beam: a prescribed source -- Gaussian across the slab (x), decaying with depth (y).
        // EVERY number here is a placeholder; the real profile is an open decision.
        beam.assign(g.size(), 0.0);
        const double S0 = ex("beam", "S0_m3s", 0), xc = ex("beam", "x_center_m", 0.5 * g.Lx);
        const double sx = ex("beam", "sigma_x_m", 0.001), dep = ex("beam", "depth_y_m", 0.005);
        beamExciteRatio = ex("beam", "excite_ratio", 0);
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                const double X = g.x(i) - xc;
                beam[g.idx(i, j)] = S0 * std::exp(-X * X / (2 * sx * sx)) * std::exp(-g.y(j) / dep);
            }
    }

    int nSpecies() const { return (int)n.size(); }
    double density(int cfgSpecies, int cell) const {
        if (c_.species[cfgSpecies].fixed) return nfix_;
        for (size_t k = 0; k < speciesIdx.size(); ++k) if (speciesIdx[k] == cfgSpecies) return std::max(n[k + 1][cell], 0.0);
        return 0.0;
    }

    // Production and loss-frequency of every species from the config's reaction list, using the
    // current densities. Loss is returned as a frequency so the step can treat it implicitly.
    void reactions(std::vector<Field>& P, std::vector<Field>& L) const {
        const int N = g.size();
        P.assign(n.size(), Field(N, 0.0));
        L.assign(n.size(), Field(N, 0.0));
        auto slot = [&](int cfgSp) {   // position in n[] of a config species; -1 = electron -> 0; fixed -> -1
            if (cfgSp == -1) return 0;
            for (size_t k = 0; k < speciesIdx.size(); ++k) if (speciesIdx[k] == cfgSp) return (int)k + 1;
            return -1;
        };
        for (int cell = 0; cell < N; ++cell) {
            const double ne = std::max(n[0][cell], 0.0);
            for (const Reaction& r : c_.reactions) {
                double rate = r.A * std::pow(Te, r.B) * std::exp(-r.E / Te);
                for (const Term& t : r.reactants) rate *= std::pow(t.species < 0 ? ne : density(t.species, cell), t.count);
                if (rate <= 0) continue;
                for (const Term& t : r.reactants) {
                    const int s = slot(t.species);
                    if (s < 0) continue;
                    L[s][cell] += t.count * rate / std::max(n[s][cell], 1.0);
                }
                for (const Term& t : r.products) {
                    const int s = slot(t.species);
                    if (s >= 0) P[s][cell] += t.count * rate;
                }
            }
        }
    }

    // One time step of length dt (see the file header for the order of operations).
    void step(double dt) {
        const int N = g.size();
        const double q = phys2::e;

        // 1. Field from a Poisson equation that already knows how the charges will move.
        if (poissonOn && g.hasWall()) {
            Field kappa(N, phys2::eps0), f(N, 0.0), lap;
            for (size_t s = 0; s < n.size(); ++s) {
                if (tr[s].s == 0) continue;
                const double mu = tr[s].D / tr[s].VT;
                for (int c = 0; c < N; ++c) {
                    kappa[c] += dt * q * std::fabs(tr[s].s) * mu * std::max(n[s][c], 0.0);
                    f[c]     += q * tr[s].s * n[s][c];
                }
                diffusionDivergence(g, tr[s].D, n[s], lap);
                for (int c = 0; c < N; ++c) f[c] += dt * q * tr[s].s * lap[c];
            }
            Stencil A;
            assemblePoisson(g, kappa, f, A);
            solveStencil(g, A, phi);
        }

        // 2. Sources and sinks from the reactions and the beam.
        std::vector<Field> P(n.size(), Field(N, 0.0)), L(n.size(), Field(N, 0.0));
        if (reactionsOn) reactions(P, L);
        if (beamOn) for (int c = 0; c < N; ++c) {
            P[0][c] += beam[c];
            if (ionIdx >= 0)  P[ionIdx][c]  += beam[c];
            if (metaIdx >= 0) P[metaIdx][c] += beamExciteRatio * beam[c];
        }

        // 3. Implicit transport of every species on the new field.
        for (size_t s = 0; s < n.size(); ++s) {
            Stencil A;
            assembleTransport(g, tr[s], phi, dt, n[s], P[s], L[s], A);
            solveStencil(g, A, n[s]);
        }
        time += dt;
    }

    double mean(int s) const { double t = 0; for (double v : n[s]) t += v; return t / g.size(); }
    // Net charge in the slab, coulombs per metre of depth.
    double totalCharge() const {
        double t = 0;
        for (size_t s = 0; s < n.size(); ++s) { double sum = 0; for (double v : n[s]) sum += v; t += tr[s].s * sum; }
        return t * phys2::e * g.cellArea();
    }

private:
    const Config& c_;
    double ngas_ = 0, nfix_ = 0;
};
