# plasma-0d

A 0D (global) model of a low-pressure plasma: species densities and electron temperature as coupled ODEs.
C++17 solver writes a CSV; Python plots it and overlays measured data.

> **Every number in `config/argon_standin.cfg` is a placeholder** (generic argon, textbook order of magnitude).
> Results are not the lab's plasma until Dr. Masoud's values replace them.

## Build

CMake:

    cmake -S . -B build
    cmake --build build

Or directly, no CMake (standard library only):

    g++ -std=c++17 -O2 -Wall -Wextra src/main.cpp -o plasma0d

## Run

    plasma0d config/argon_standin.cfg output.csv
    python python/plot.py output.csv                       # needs matplotlib
    python python/plot.py output.csv data/my_data.csv      # overlay measured data

## Files

| file | role |
|---|---|
| `src/config.hpp` | config parser (hand-written, no libraries) |
| `src/model.hpp` | the physics: reactions, wall losses, energy balance |
| `src/solver.hpp` | adaptive Rosenbrock (ROS2) stiff integrator |
| `src/main.cpp` | driver: read config, integrate, write CSV |
| `python/plot.py` | plots and overlays |
| `data/measured_template.csv` | measured-data format |
| `src/grid2d.hpp`, `src/poisson2d.hpp`, `src/fluid2d.hpp`, `src/main2d.cpp` | v0.2.0 2D model (see the 2D section below) |
| `config/argon2d_standin.cfg` | 2D stand-in case (all placeholders) |
| `python/plot2d.py` | field maps and history for the 2D run |
| `tests/check1` to `check6` | verification checks, run with `ctest` |
| `docs/walkthrough_v0.2.0.md` | code walkthrough: science to files |

## What the model solves

State: density of every non-fixed species, plus Te. Background neutrals are fixed from pressure and gas temperature (ideal gas, assumed undepleted).

- Electron density from charge neutrality: `n_e = sum(charge * n)`.
- Reaction rate: `k(Te) = A Te^B exp(-E/Te)`; rate = `k * product of reactant densities` (electrons use n_e).
- Electron energy: `d/dt(3/2 n_e Te) = P/(V e) - sum(dE * rate) - wall energy loss`.
- Ion wall loss: Bohm flux over an effective area (Lieberman & Lichtenberg edge factors hL, hR), each ion carries off `2 Te + Vs`, with `Vs = (Te/2) ln(M / 2 pi m_e)`.
- Neutral (metastable) wall loss: diffusion, `rate = n D / Lambda^2`, `D = Dp / pressure`, cylinder `Lambda`.

## Config format

Sections `[global]`, `[species]`, `[reactions]`, `[walls]`; `#` starts a comment.

**Add a species:** add a line in `[species]` (`name : charge=.. mass=.. n0=..`), then reactions and a wall line for it. No C++ change. Negative ions are not supported yet (n_e comes from neutrality).

## Units

Densities m^-3, Te and energies eV, rate constants m^3/s, pressure Pa, power W (absorbed by the plasma), lengths m, time s.

## Output CSV

`time_s, <species>_m3 ..., n_e_m3, Te_eV` (first row is t = 0, then log-spaced).

## Measured data format

See `data/measured_template.csv`. Same units as the output; **densities in m^-3, not cm^-3**. Columns that share a name with an output column are overlaid. Fitting parameters to data is a later stage and not built yet.

## Assumptions and limits

- 0D, uniform plasma; Maxwellian electrons; no spatial effects.
- Power is power absorbed by the plasma.
- Background gas undepleted; several fixed species would share the gas density equally.
- Te is floored at 0.05 eV inside the rate evaluation only to keep the math finite.
- Stand-in cross sections, rates and Dp are not from the lab.

## 2D model (v0.2.0)

> **Every number in `config/argon2d_standin.cfg` is a placeholder.** The pressure, mobilities and beam were chosen so the run is numerically comfortable. They are not the lab's plasma and not a model of the CNT beam plume.

A 2D Cartesian slab (x-y, per unit depth). Per species `n_k`: `dn/dt + div(Gamma) = S`, with `Gamma = s mu n E - D grad n`, `E = -grad phi`, `-div(eps0 grad phi) = rho`, `D = mu * (temperature in volts)`. Scharfetter-Gummel flux, implicit transport, field solved first with the plasma's response folded in (semi-implicit, first order in time).

Build and run (the static option avoids a Windows DLL mix-up on laptops with two MinGW installs):

    cmake -S . -B build2 -DPLASMA_STATIC_RUNTIME=ON
    cmake --build build2
    build2/plasma2d config/argon2d_standin.cfg history2d.csv fields2d.csv
    python python/plot2d.py fields2d.csv history2d.csv --save fields2d.png
    ctest --test-dir build2 --output-on-failure        # the six verification checks

Or directly: `g++ -std=c++17 -O2 -Wall -Wextra -static src/main2d.cpp -o plasma2d`.

New config sections (numeric `key = value`): `[grid]` (`nx`, `ny`, `Lx_m`, `Ly_m`, side types `bc_xmin` etc.: 0 wall, 1 mirror), `[transport]` (`mu_e_m2Vs`, `mu_i_m2Vs`), `[beam]` (`S0_m3s`, `x_center_m`, `sigma_x_m`, `depth_y_m`, `excite_ratio`), `[solver2d]` (`dt_s`). Species, reactions, `Te0_eV`, `pressure_Pa`, `t_end_s` use the 0D syntax; neutral diffusion comes from the `[walls]` line.

Output: `history2d.csv` (slab averages and net charge per time) and `fields2d.csv` (`x_m, y_m`, one column per species, `phi_V`).

Limits: Te is constant; no electron energy equation; constant mobilities (no Boltzmann-solver tables); cells larger than the Debye length near walls, so wall values are not physical; the metastable density is still rising at the end of the stand-in run; no comparison with measured data yet. Walkthrough: `docs/walkthrough_v0.2.0.md`.
