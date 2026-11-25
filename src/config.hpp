// config.hpp -- reads the plain-text config into plain structs.
// Hand-written parser on purpose: no JSON/YAML library means nothing to install.
#pragma once
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Species {
    std::string name;
    double charge = 0, mass_amu = 0;
    bool fixed = false;   // density set from pressure, not integrated
    double n0 = 0;        // initial density (evolved species)
};

struct Term { int species; int count; };   // species index, or -1 for the electron

struct Reaction {
    std::string name;
    std::vector<Term> reactants, products;
    double A = 0, B = 0, E = 0;     // k = A Te^B exp(-E/Te)
    double dE = 0, dEte = 0;        // eV lost per event; extra loss = dEte * Te
};

struct Wall { int species; bool bohm; double Dp; };

struct Config {
    // [global]
    double pressure = 1, Tgas = 300, power = 100, R = 0.1, L = 0.2, sigma_i = 1e-18;
    double Te0 = 2, t_end = 1e-2, t_first = 1e-9, rtol = 1e-5, atol_n = 1e6, atol_Te = 1e-7;
    double h0 = 1e-12, hmax = 1e-3;
    int n_outputs = 200;
    std::vector<Species> species;
    std::vector<Reaction> reactions;
    std::vector<Wall> walls;
    // v0.2.0: sections only the 2D code reads ([grid], [transport], [beam], [solver2d]).
    // Kept as a generic map so the 0D parser and every 0D config behave exactly as before.
    std::map<std::string, std::map<std::string, double>> extra;

    int find(const std::string& n) const {
        for (size_t i = 0; i < species.size(); ++i) if (species[i].name == n) return (int)i;
        return -2;
    }
};

namespace cfg_detail {

inline std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

// key=value tokens separated by spaces, e.g. "A=2.3e-14 B=0.5"
inline std::map<std::string, std::string> kv(const std::string& s) {
    std::map<std::string, std::string> m;
    std::istringstream in(s);
    std::string tok;
    while (in >> tok) {
        size_t p = tok.find('=');
        if (p != std::string::npos) m[tok.substr(0, p)] = tok.substr(p + 1);
    }
    return m;
}

inline double num(const std::map<std::string, std::string>& m, const std::string& k, double d) {
    auto it = m.find(k);
    return it == m.end() ? d : std::stod(it->second);
}

// "Ar+ + 2 e" -> terms. Species may contain '+' in their names (Ar+), so we split
// on " + " with surrounding spaces only.
inline std::vector<Term> side(const std::string& s, const Config& c, int line) {
    std::vector<Term> out;
    size_t pos = 0;
    std::string str = trim(s);
    while (pos <= str.size()) {
        size_t nx = str.find(" + ", pos);
        std::string t = trim(str.substr(pos, nx == std::string::npos ? std::string::npos : nx - pos));
        if (!t.empty()) {
            int count = 1;
            std::istringstream in(t);
            std::string a, b;
            in >> a >> b;
            std::string name = a;
            if (!b.empty()) { count = std::stoi(a); name = b; }
            int idx = (name == "e") ? -1 : c.find(name);
            if (idx == -2) throw std::runtime_error("line " + std::to_string(line) + ": unknown species '" + name + "'");
            out.push_back({idx, count});
        }
        if (nx == std::string::npos) break;
        pos = nx + 3;
    }
    return out;
}

} // namespace cfg_detail

inline Config load_config(const std::string& path) {
    using namespace cfg_detail;
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open config file: " + path);
    Config c;
    std::string line, sect;
    int ln = 0;
    while (std::getline(f, line)) {
        ++ln;
        size_t h = line.find('#');
        if (h != std::string::npos) line.erase(h);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[') { sect = line; continue; }

        if (sect == "[global]") {
            size_t eq = line.find('=');
            if (eq == std::string::npos) throw std::runtime_error("line " + std::to_string(ln) + ": expected key = value");
            std::string k = trim(line.substr(0, eq)); double v = std::stod(trim(line.substr(eq + 1)));
            if (k == "pressure_Pa") c.pressure = v; else if (k == "gas_temp_K") c.Tgas = v;
            else if (k == "power_W") c.power = v; else if (k == "radius_m") c.R = v;
            else if (k == "length_m") c.L = v; else if (k == "ion_cross_section_m2") c.sigma_i = v;
            else if (k == "Te0_eV") c.Te0 = v; else if (k == "t_end_s") c.t_end = v;
            else if (k == "t_first_s") c.t_first = v; else if (k == "n_outputs") c.n_outputs = (int)v;
            else if (k == "rtol") c.rtol = v; else if (k == "atol_density") c.atol_n = v;
            else if (k == "atol_Te") c.atol_Te = v; else if (k == "h_initial_s") c.h0 = v;
            else if (k == "h_max_s") c.hmax = v;
            else throw std::runtime_error("line " + std::to_string(ln) + ": unknown [global] key '" + k + "'");
        } else if (sect == "[species]") {
            size_t co = line.find(':');
            if (co == std::string::npos) throw std::runtime_error("line " + std::to_string(ln) + ": expected name : key=value ...");
            Species s; s.name = trim(line.substr(0, co));
            auto m = kv(line.substr(co + 1));
            s.charge = num(m, "charge", 0); s.mass_amu = num(m, "mass", 0);
            s.fixed = num(m, "fixed", 0) != 0; s.n0 = num(m, "n0", 0);
            if (s.charge < 0) throw std::runtime_error("negative ions are not supported (n_e comes from charge neutrality)");
            c.species.push_back(s);
        } else if (sect == "[reactions]") {
            size_t co = line.find(':'), ar = line.find("=>"), sc = line.find(';');
            if (co == std::string::npos || ar == std::string::npos || sc == std::string::npos)
                throw std::runtime_error("line " + std::to_string(ln) + ": expected name : lhs => rhs ; A=.. B=.. E=..");
            Reaction r; r.name = trim(line.substr(0, co));
            r.reactants = side(line.substr(co + 1, ar - co - 1), c, ln);
            r.products = side(line.substr(ar + 2, sc - ar - 2), c, ln);
            auto m = kv(line.substr(sc + 1));
            r.A = num(m, "A", 0); r.B = num(m, "B", 0); r.E = num(m, "E", 0);
            r.dE = num(m, "dE", 0); r.dEte = num(m, "dEte", 0);
            c.reactions.push_back(r);
        } else if (sect == "[walls]") {
            size_t co = line.find(':');
            if (co == std::string::npos) throw std::runtime_error("line " + std::to_string(ln) + ": expected name : bohm | diffusion Dp=..");
            Wall w; std::string nm = trim(line.substr(0, co));
            w.species = c.find(nm);
            if (w.species < 0) throw std::runtime_error("line " + std::to_string(ln) + ": unknown species '" + nm + "'");
            std::string rest = trim(line.substr(co + 1));
            w.bohm = rest.rfind("bohm", 0) == 0;
            if (!w.bohm && rest.rfind("diffusion", 0) != 0)
                throw std::runtime_error("line " + std::to_string(ln) + ": wall type must be bohm or diffusion");
            w.Dp = num(kv(rest), "Dp", 0);
            c.walls.push_back(w);
        } else if (sect == "[grid]" || sect == "[transport]" || sect == "[beam]" || sect == "[solver2d]") {
            size_t eq = line.find('=');
            if (eq == std::string::npos) throw std::runtime_error("line " + std::to_string(ln) + ": expected key = value");
            c.extra[sect.substr(1, sect.size() - 2)][trim(line.substr(0, eq))] = std::stod(trim(line.substr(eq + 1)));
        } else {
            throw std::runtime_error("line " + std::to_string(ln) + ": text outside a known [section]");
        }
    }
    return c;
}
