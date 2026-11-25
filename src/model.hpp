// model.hpp -- the physics: turns a Config into the right-hand side dy/dt.
//
// State vector y = [ n_s for every NON-fixed species (m^-3) ..., Te (eV) ].
// Electron density is not in the state: it follows from charge neutrality,
//     n_e = sum_s charge_s * n_s
// which removes one stiff equation and guarantees the plasma stays neutral.
//
// Electron energy balance (the Te equation), all in eV m^-3 s^-1:
//     d/dt( 3/2 n_e Te ) = P/(V e)  -  sum_r dE_r * R_r  -  wall losses
// Expanding the left side gives dTe/dt once dn_e/dt is known (see rhs()).
#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "config.hpp"
#include "solver.hpp"   // for Vec

namespace phys {
constexpr double kB   = 1.380649e-23;     // J/K
constexpr double e    = 1.602176634e-19;  // C (also J per eV)
constexpr double amu  = 1.66053906660e-27;// kg
constexpr double me   = 9.1093837015e-31; // kg
constexpr double pi   = 3.14159265358979323846;
}

class Model {
public:
    explicit Model(const Config& c) : c_(c) {
        const int ns = (int)c.species.size();
        idx_.assign(ns, -1);
        int nfixed = 0;
        for (int i = 0; i < ns; ++i) {
            if (c.species[i].fixed) ++nfixed;
            else { idx_[i] = (int)evolved_.size(); evolved_.push_back(i); }
        }
        if (nfixed == 0) throw std::runtime_error("need at least one fixed (background gas) species");
        // Background gas density from the ideal gas law. Several fixed species
        // share it equally; refine here if a mixture is ever needed.
        ngas_ = c.pressure / (phys::kB * c.Tgas);
        nfix_ = ngas_ / nfixed;
        V_ = phys::pi * c.R * c.R * c.L;
        lambda_i_ = 1.0 / (ngas_ * c.sigma_i);
        // Edge-to-center density ratios (Lieberman & Lichtenberg, low pressure).
        hL_ = 0.86 / std::sqrt(3.0 + c.L / (2.0 * lambda_i_));
        hR_ = 0.80 / std::sqrt(4.0 + c.R / lambda_i_);
        Aeff_ = 2.0 * phys::pi * c.R * c.R * hL_ + 2.0 * phys::pi * c.R * c.L * hR_;
        // Fundamental diffusion length of a cylinder: 1/L^2 = (pi/L)^2 + (2.405/R)^2.
        const double a = phys::pi / c.L, b = 2.405 / c.R;
        Lambda2_ = 1.0 / (a * a + b * b);
    }

    int n_state() const { return (int)evolved_.size() + 1; }
    const std::vector<int>& evolved() const { return evolved_; }

    Vec initial() const {
        Vec y(n_state());
        for (size_t k = 0; k < evolved_.size(); ++k) y[k] = c_.species[evolved_[k]].n0;
        y.back() = c_.Te0;
        return y;
    }

    // Full density of species i (fixed ones come from pressure).
    double density(const Vec& y, int i) const {
        return c_.species[i].fixed ? nfix_ : std::max(y[idx_[i]], 0.0);
    }

    double n_e(const Vec& y) const {
        double ne = 0;
        for (size_t i = 0; i < c_.species.size(); ++i) ne += c_.species[i].charge * density(y, (int)i);
        // Floor keeps Te's equation finite if every ion briefly hits zero.
        return std::max(ne, 1.0);
    }

    void rhs(const Vec& y, Vec& dy) const {
        const int ns = (int)c_.species.size();
        const double Te = std::max(y.back(), 0.05);  // keeps exp(-E/Te) and Te^B sane
        const double ne = n_e(y);

        std::vector<double> dn(ns, 0.0);
        double Eloss = 0;   // eV m^-3 s^-1 leaving the electron gas

        for (const Reaction& r : c_.reactions) {
            double rate = r.A * std::pow(Te, r.B) * std::exp(-r.E / Te);   // k(Te)
            bool has_e = false;
            for (const Term& t : r.reactants) {
                if (t.species < 0) { has_e = true; rate *= std::pow(ne, t.count); }
                else rate *= std::pow(density(y, t.species), t.count);
            }
            (void)has_e;
            // Electrons are derived, so only heavy species enter the bookkeeping.
            for (const Term& t : r.reactants) if (t.species >= 0) dn[t.species] -= t.count * rate;
            for (const Term& t : r.products)  if (t.species >= 0) dn[t.species] += t.count * rate;
            Eloss += (r.dE + r.dEte * Te) * rate;
        }

        for (const Wall& w : c_.walls) {
            const Species& s = c_.species[w.species];
            const double n = density(y, w.species);
            if (w.bohm) {
                // Flux n u_B through the effective wall area, per unit volume.
                const double M  = s.mass_amu * phys::amu;
                const double uB = std::sqrt(phys::e * Te / M);
                const double loss = n * uB * Aeff_ / V_;
                dn[w.species] -= loss;
                // Each ion lost takes 2 Te (electron) plus the sheath drop
                // Vs = (Te/2) ln(M / 2 pi me) of kinetic energy to the wall.
                const double Vs = 0.5 * Te * std::log(M / (2.0 * phys::pi * phys::me));
                Eloss += (2.0 * Te + Vs) * loss;
            } else {
                const double D = w.Dp / c_.pressure;           // m^2/s
                dn[w.species] -= n * D / Lambda2_;             // n / tau, tau = Lambda^2 / D
            }
        }

        dy.assign(n_state(), 0.0);
        for (size_t k = 0; k < evolved_.size(); ++k) dy[k] = dn[evolved_[k]];

        double dne = 0;
        for (int i = 0; i < ns; ++i) if (!c_.species[i].fixed) dne += c_.species[i].charge * dn[i];

        const double heat = c_.power / (V_ * phys::e);         // eV m^-3 s^-1
        // 3/2 (ne dTe + Te dne) = heat - Eloss
        dy.back() = ((heat - Eloss) / 1.5 - Te * dne) / ne;
    }

private:
    const Config& c_;
    std::vector<int> evolved_, idx_;
    double ngas_ = 0, nfix_ = 0, V_ = 0, lambda_i_ = 0, hL_ = 0, hR_ = 0, Aeff_ = 0, Lambda2_ = 0;
};
