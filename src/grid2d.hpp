// grid2d.hpp -- the 2D Cartesian slab grid (x-y) used by the v0.2.0 fluid model.
//
// Cell-centered finite volumes: every unknown lives at a cell center, every flux
// lives on a cell face. That choice is what makes "conservation" checkable --
// whatever leaves one cell enters its neighbor, so only the walls can lose charge.
//
// Per unit depth: the slab is infinite in z (a line source), so areas are just
// face lengths and volumes are dx*dy.
#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

typedef std::vector<double> Field;   // one value per cell, index = j*nx + i

enum Side   { XMIN = 0, XMAX = 1, YMIN = 2, YMAX = 3 };
enum BCType { WALL = 0, SYM = 1 };   // WALL: absorbing and grounded. SYM: mirror, no flux.

struct Grid {
    int nx = 0, ny = 0;
    double Lx = 0, Ly = 0, dx = 0, dy = 0;
    BCType bc[4]      = {WALL, WALL, WALL, WALL};
    double phiWall[4] = {0, 0, 0, 0};   // potential on a WALL face, V
    double nWall[4]   = {0, 0, 0, 0};   // density on a WALL face, m^-3 (0 = perfectly absorbing)

    void init(int nx_, int ny_, double Lx_, double Ly_) {
        if (nx_ < 2 || ny_ < 2) throw std::runtime_error("grid needs at least 2 x 2 cells");
        nx = nx_; ny = ny_; Lx = Lx_; Ly = Ly_; dx = Lx / nx; dy = Ly / ny;
    }
    int size() const { return nx * ny; }
    int idx(int i, int j) const { return j * nx + i; }
    double x(int i) const { return (i + 0.5) * dx; }
    double y(int j) const { return (j + 0.5) * dy; }
    double cellArea() const { return dx * dy; }
    bool hasWall() const { for (int s = 0; s < 4; ++s) if (bc[s] == WALL) return true; return false; }
};

// One face of one cell. nb >= 0: interior face, nb is the neighbor cell.
// nb < 0: boundary face on `side`. h is the center-to-center distance across an
// interior face, or the full cell width normal to a boundary face (the wall sits
// h/2 from the cell center).
struct Face { int nb; int side; double h; };

// f: 0 = east (+x), 1 = west (-x), 2 = north (+y), 3 = south (-y)
inline Face faceOf(const Grid& g, int i, int j, int f) {
    switch (f) {
        case 0: return i + 1 < g.nx ? Face{g.idx(i + 1, j), -1, g.dx} : Face{-1, XMAX, g.dx};
        case 1: return i > 0        ? Face{g.idx(i - 1, j), -1, g.dx} : Face{-1, XMIN, g.dx};
        case 2: return j + 1 < g.ny ? Face{g.idx(i, j + 1), -1, g.dy} : Face{-1, YMAX, g.dy};
        default: return j > 0       ? Face{g.idx(i, j - 1), -1, g.dy} : Face{-1, YMIN, g.dy};
    }
}
