# Changelog

## 0.2.0 - 2026-11-25
2D fluid model on a Cartesian slab (x-y), verified against known answers. All physical numbers are still stand-ins.
- New `plasma2d` program: drift-diffusion for electrons, Ar+ and Ar*, Poisson equation for the field, Scharfetter-Gummel flux, implicit transport, semi-implicit field solve (first order in time). Own C++, no libraries.
- Prescribed beam source (Gaussian across x, decaying in y) as a labeled placeholder; same species and reactions as 0.1.0 plus electrons as their own species.
- Te is held constant (`Te0_eV`); the electron energy equation is not built yet.
- `config/argon2d_standin.cfg`; `python/plot2d.py`; `docs/walkthrough_v0.2.0.md` maps the science to the files.
- Six verification checks in `tests/` (`ctest`): uniform reduces to the 0D reaction ODE (first-order convergence), pure diffusion, drift with diffusion, Scharfetter-Gummel steady state, Poisson with Gauss and charge conservation, manufactured solution (second order). All pass.
- `src/config.hpp`: additive only (`[grid]`, `[transport]`, `[beam]`, `[solver2d]` sections). The 0D output for `config/argon_standin.cfg` is byte-identical to 0.1.0.
- CMake: new targets, plus optional `-DPLASMA_STATIC_RUNTIME=ON` (also applies to `plasma0d`) for Windows setups where the C++ runtime DLLs clash.
- Still placeholders: every number in `config/argon2d_standin.cfg` (pressure, mobilities, beam amplitude and shape, grid size). Not the CNT beam plume, not the lab's plasma.
- Not built: electron energy equation, tabulated transport and rate coefficients, any comparison with measured data. Open: argon levels, whether the target is a CNT beam plume at hundreds of mbar.

## 0.1.0 - 2026-11-25
First runnable build. All physical numbers are stand-ins.
- Config-driven 0D global model: species, reactions, walls, discharge conditions in `config/argon_standin.cfg`.
- State: evolved species densities plus Te; n_e from charge neutrality; background gas fixed from pressure.
- Adaptive Rosenbrock (ROS2) integrator; log-spaced CSV output.
- `python/plot.py` plots model output and overlays measured data.
- `data/measured_template.csv` documents the measured-data format.
- Stand-in argon case reaches steady state (n_e about 4e17 m^-3, Te about 1.9 eV) in 0.5 s of model time. Not the lab's plasma.
