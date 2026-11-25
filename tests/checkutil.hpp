// checkutil.hpp -- the tiny shared bit of every verification program in tests/.
// Each check prints its numbers, then PASS or FAIL, and returns non-zero on FAIL
// so a script (or CMake's ctest) can tell without reading.
//
// Tolerances live in each check's source, next to the reason for the number.
#pragma once
#include <cmath>
#include <cstdio>

struct Checker {
    int fails = 0;
    void expect(bool ok, const char* what) {
        std::printf("  [%s] %s\n", ok ? "pass" : "FAIL", what);
        if (!ok) ++fails;
    }
    int finish(const char* name) {
        std::printf("%s: %s\n", name, fails ? "FAIL" : "PASS");
        return fails ? 1 : 0;
    }
};

// Observed order of accuracy between a coarse and a fine grid (each refinement halves h).
inline double observedOrder(double errCoarse, double errFine) { return std::log2(errCoarse / errFine); }
