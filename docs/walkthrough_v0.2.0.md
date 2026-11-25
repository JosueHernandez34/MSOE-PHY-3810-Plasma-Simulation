# v0.2.0 code walkthrough: from the science to the files

Companion to the Learning Guide, Stage 1 (2D fluid model). Read the Stage 1 papers first, then
open the file named in each section. All numbers in the stand-in config are placeholders.

## The picture in one paragraph

v0.1.0 treated the plasma as one well-mixed box. v0.2.0 puts it on a grid. Every cell holds a
density for electrons, Ar+ and Ar*. Charged particles drift in the electric field and diffuse;
metastables only diffuse. The electric field comes from the charge itself (Poisson). The beam is
a prescribed source of electron-ion pairs. The electron temperature is a constant for now.

## Idea to file map

| Idea | Where |
|---|---|
| The grid: cells, faces, neighbors, wall or mirror sides | `src/grid2d.hpp` (`Grid`, `faceOf`) |
| One equation per cell linking it to 4 neighbors, and the solver that handles all of them | `src/poisson2d.hpp` (`Stencil`, `applyStencil`, `solveStencil`) |
| Poisson: charge makes the field, `-div(kappa grad phi) = f` | `src/poisson2d.hpp` (`assemblePoisson`) |
| Gauss's law as a check: charge inside equals flux through the walls | `src/poisson2d.hpp` (`wallFlux`), tests/check5 |
| Drift and diffusion through one face, the Scharfetter-Gummel flux | `src/fluid2d.hpp` (`bernoulli`, `assembleTransport`) |
| Continuity: what enters a cell minus what leaves it, plus sources and losses | `src/fluid2d.hpp` (`assembleTransport`) |
| Einstein relation `D = mu * T` | `src/fluid2d.hpp`, `Fluid2D` constructor (`Transport` values) |
| Reactions as sources and sinks per cell | `src/fluid2d.hpp` (`Fluid2D::reactions`) |
| The beam as a prescribed source | `src/fluid2d.hpp`, `Fluid2D` constructor (`beam`) |
| The time loop: field first, then sources, then transport | `src/fluid2d.hpp` (`Fluid2D::step`) |
| Reading the new config sections | `src/config.hpp` (`extra`), `config/argon2d_standin.cfg` |
| Running it and writing CSVs | `src/main2d.cpp` |
| Plotting | `python/plot2d.py` |

## The pieces, one at a time

### Cells and faces (`grid2d.hpp`)

Unknowns sit at cell centers; fluxes live on faces. Whatever leaves one cell enters its neighbor,
so the only places charge can disappear are wall faces. `faceOf` returns, for one face of one
cell, either the neighbor across it or the wall on that side, plus the distance `h`. Every loop
in the code is "for each cell, for each of 4 faces", using this one helper, so wall handling is
written once.

A `WALL` side is absorbing and grounded (potential fixed, density fixed, normally 0). A `SYM`
side is a mirror: nothing crosses it. Tests use mirrors to build 1D problems and closed boxes.

### The Scharfetter-Gummel flux (`assembleTransport`)

Across a face the particle flux is drift plus diffusion. A plain difference scheme goes wrong when
drift dominates (negative densities, wiggles). Scharfetter-Gummel solves the 1D drift-diffusion
equation exactly across each face and uses that:

```
flux out of cell C into neighbor N = (D / h) * [ B(-a) n_C - B(a) n_N ]
a = s (phi_C - phi_N) / V_T          B(x) = x / (e^x - 1)
```

`a` is the voltage drop across the face measured in thermal volts, so it does not need the
mobility separately (the Einstein relation `D = mu V_T` cancels it). When `a` is near 0, `B` is near
1 and this is the ordinary diffusion difference. When `a` is large, one `B` goes to 0 and the
flux becomes pure upwind drift. That is why one formula covers both regimes. Check 4 shows it is
exact in 1D even at drift 20 times stronger than diffusion over the box.

Each cell's equation is backward Euler (implicit): `(n_new - n_old)/dt + outflow + loss*n_new = source`.
Implicit means no limit on dt from diffusion. Sources are explicit and losses implicit, so every
off-diagonal is negative and the diagonal dominates: densities cannot go negative.

### Poisson, and why it is solved first (`Fluid2D::step`, step 1)

The field is the slow-to-get-right part. Electrons respond in about 1e-10 s, far faster than any
time step we want, so a naive scheme (compute field from old charge, move charge, repeat) blows
up. The semi-implicit trick: solve a Poisson equation whose coefficient `kappa` already includes
how much the charges will move during this step,

```
kappa = eps0 + dt * sum_k |q_k| mu_k n_k
f     = rho + dt * div( sum_k q_k D_k grad n_k )
```

This removes the limit on dt from dielectric relaxation. The price is first order accuracy in
time, which check 1 measures directly.

### Reactions (`Fluid2D::reactions`)

Same reaction list and same `A Te^B exp(-E/Te)` rate form as v0.1.0, read by the same config
parser. For each cell and reaction the rate is computed once, then added to the production of
the products and to the loss frequency of the reactants. Loss goes in as a frequency so the step
can treat it implicitly. Te is the constant `Te0_eV` from the config.

### The beam (`Fluid2D` constructor)

A Gaussian across x times an exponential decay in y, scaled by `S0_m3s`. It creates electron-ion
pairs, plus metastables at `excite_ratio` per pair. This is a placeholder shape. A real source
would come from an energy-deposition calculation, which is still an open item.

### The time loop (`Fluid2D::step`)

1. Poisson: new potential from the old densities (skipped if the box has no wall at all, since
   then the potential is undetermined).
2. Reactions and beam: production and loss in every cell.
3. Transport: one implicit solve per species on the new potential.

### The linear solve (`solveStencil`)

BiCGStab with a diagonal preconditioner. It serves both Poisson (symmetric) and transport
(not symmetric). It throws if it fails to converge, because a wrong field would look like physics.
Convergence tolerance is 1e-14 relative, which also sets how exactly charge is conserved.

## How the checks map to the science

| Check | Idea it protects |
|---|---|
| 1 uniform reduces to 0D | the reaction terms are the same chemistry as v0.1.0 |
| 2 pure diffusion | the diffusion difference and the wall treatment |
| 3 drift with diffusion | speed and spreading of a moving cloud (Scharfetter-Gummel) |
| 4 strong-drift steady state | Scharfetter-Gummel exactness |
| 5 Poisson, Gauss, charge | the field solve and charge bookkeeping |
| 6 manufactured solution | everything together, 2D, with the right order of accuracy |

## What to look for when you run it

```
build\plasma2d.exe config\argon2d_standin.cfg history2d.csv fields2d.csv
python python\plot2d.py fields2d.csv history2d.csv --save fields2d.png
```

- Electrons and Ar+ rise together, then level off within about a millisecond.
- The plasma potential rises a few volts above the walls: electrons are fast and leave first, the
  slower ions are held back, and the field that results keeps electrons in. That is the
  ambipolar field, and it is the main new thing 2D shows over 0D.
- Ar* is still rising at the end of the run: metastables diffuse slowly (about 0.1 s to settle).

## Not in v0.2.0 (so do not look for it)

- Electron energy equation: Te is constant. (Open decision.)
- Tabulated rate and transport coefficients from a Boltzmann solver.
- Sheath resolution: cells are larger than the Debye length near the walls, which is normal for a
  fluid code and why wall values should not be read as physical.
- Any comparison with the lab's data.
