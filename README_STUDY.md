# Quasi-1D nozzle grid/CFL study

The study uses `nozzle_fvm.c`, a first-order cell-centered finite-volume solver
for the quasi-one-dimensional compressible Euler equations. It uses a Rusanov
flux, a geometric momentum source term, and explicit CFL-controlled pseudo-time
marching.

## Build and run

```bash
cc -std=c11 -O2 -Wall -Wextra -pedantic nozzle_fvm.c -lm -o nozzle_fvm
./nozzle_fvm cases/baseline.cfg
```

Each run writes:

- `<output_prefix>.csv`: x, area, density, velocity, pressure, Mach number
- `<output_prefix>_residual.csv`: time step, residual, and time step size

## Generate the study data

```bash
for f in cases/grid_*.cfg; do ./nozzle_fvm "$f"; done
for f in cases/cfl_*.cfg; do ./nozzle_fvm "$f"; done
```

The grid study uses 100, 200, and 400 cells at CFL=0.5. The CFL study uses
200 cells at CFL=0.2, 0.5, 0.8, and 1.0. All cases use the same geometry and
thermodynamic conditions so that only the intended parameter changes.

## Input conditions

The working fluid is air with `gamma=1.4` and `R=287.05 J/(kg K)`. The inlet
stagnation state is `p0=101325 Pa`, `T0=300 K`; the outlet back pressure is
`p_back=30000 Pa`. The nozzle is 1 m long, with its throat at x=0.4 m and
areas `A_in=0.020 m^2`, `A_throat=0.010 m^2`, and `A_exit=0.016 m^2`.

## Quantities to extract

For each final solution, compare pressure and Mach profiles. Track the largest
Mach-number gradient as an approximate internal shock position when a shock is
present. From the residual files compare iteration count, residual decay, and
stability. For grid refinement, compare the 200-cell result with the 400-cell
result using an L1 or L2 difference in pressure and Mach number.
