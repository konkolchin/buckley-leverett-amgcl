# 1D Buckley–Leverett with AMGCL

Sequential implicit solution of a one-dimensional waterflood. Saturation is the Buckley–Leverett equation. Pressure is the incompressible total-velocity equation. Both linear systems are solved with [AMGCL](https://github.com/ddemidov/amgcl) 1.4.3 (header-only, MIT).

Dimensionless variables: \(x\in[0,1]\), \(t\) in pore volumes injected, total velocity fixed at 1 by the inlet rate.

\[
\frac{\partial S}{\partial t}+\frac{\partial f(S)}{\partial x}=0,
\qquad
\frac{d}{dx}\!\left(k(x)\,\lambda_t(S)\frac{dp}{dx}\right)=0,
\qquad
-k\lambda_t\frac{dp}{dx}=1,
\qquad
p(1)=0.
\]

Corey relative permeabilities, \(S_{wc}=S_{or}=0.2\), \(n_w=n_o=2\), \(\mu_o/\mu_w=2\). A tight streak \(k=0.1\) on \(x\in[0.45,0.55]\) steepens the pressure gradient. With a prescribed inlet rate the streak does not change saturation: in 1D incompressible flow the total velocity is fixed by the boundary flux.

## Notes

Two LaTeX notes in `docs/` cover the same ground in more detail:

- `docs/theory.tex` — Darcy fractional flow, Corey $f'(S)$, the Welge tangent, and the bidiagonal Jacobian
- `docs/solver.tex` — the sequential time step, the AMGCL types, and the checks

Compile with `pdflatex` from `docs/`.

## Discrete scheme

Saturation is a fully implicit, first-order upwind finite volume. The Newton matrix is lower bidiagonal:

\[
R_i=\frac{\Delta x}{\Delta t}(S_i-S_i^n)+f(S_i)-f(S_{i-1}),
\quad
J_{i,i}=\frac{\Delta x}{\Delta t}+f'(S_i),
\quad
J_{i,i-1}=-f'(S_{i-1}).
\]

The inlet flux is \(f(S_{\max})=1\). This matrix is nonsymmetric, so AMGCL solves it with **BiCGStab** preconditioned by **smoothed-aggregation AMG** and **ILU(0)** relaxation.

Pressure is a standard two-point flux cell-centered scheme with harmonic transmissibilities and \(p=0\) on the outlet face. The matrix is SPD, so AMGCL solves it with **CG** preconditioned by **smoothed-aggregation AMG** and **SPAI(0)** relaxation. The same AMG types are what a reservoir simulator uses on the pressure block; a fully coupled CPR preconditioner (`amgcl::preconditioner::cpr`) is the industrial wrapper around that pressure AMG, and is not required for this sequential split.

`coarse_enough` is set to 60 so the hierarchy printed for this 1D grid actually has several levels. AMGCL’s default would keep coarsening until a direct skyline LU is cheap, which on a few hundred unknowns can collapse to one level.

On this grid the transport solve reports one BiCGStab iteration. That is the expected result: a lower bidiagonal matrix is factored exactly by ILU(0), so the AMG V-cycle is an exact preconditioner. The pressure matrix is tridiagonal, ILU(0)/SPAI(0) is not an exact inverse, and CG takes on the order of ten iterations. Inlet pressure rises during the flood because \(\lambda_t(S)\) has a minimum between the endpoints (both phases are partially mobile), so a fixed rate needs a larger pressure drop. The tight streak adds a steep gradient on \(x\in[0.45,0.55]\).

## Checks

- \(f'(S)\) against a finite-difference derivative.
- Welge tangent condition \(f'(S_f)=(f(S_f)-f(S_{wc}))/(S_f-S_{wc})\).
- Saturation against that analytical shock-plus-rarefaction solution.
- Pressure against the exact solution of the discrete TPFA system (face flux identically 1).
- Conservation: water gained in the core equals injected minus produced. The implicit finite-volume scheme satisfies this by construction.

## Build

AMGCL 1.4.3 is the Git submodule in `third_party/amgcl`. Clone the project and that submodule together:

```bat
git clone --recurse-submodules https://github.com/konkolchin/buckley-leverett-amgcl.git
cd buckley-leverett-amgcl
```

If the project is already cloned and `third_party/amgcl` is empty, fetch the submodule with:

```bat
git submodule update --init
```

From a Visual Studio 2022 x64 developer prompt, in the project directory:

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
build\Release\buckley.exe
```

Optional arguments: `buckley [cells] [pore_volumes] [cfl]`. Defaults are 800 cells, a time that places the shock near \(x=0.55\), and CFL = 1 based on the shock speed. Boost is not required.

Outputs in `output/`:

- `report.txt` — AMG hierarchy, iteration counts, errors
- `profile.csv` — \(S\), Welge \(S\), \(p\), exact discrete \(p\)
- `fractional_flow.csv` — \(f(S)\) and the Welge convex hull
- `profile.svg` — saturation and pressure plots

## What this is

A study solver for the classical 1D waterflood, small enough to derive on a whiteboard and large enough that the AMGCL calls are the real ones: `amgcl::make_solver`, `amgcl::amg`, smoothed aggregation, SPAI(0), ILU(0), CG, and BiCGStab, with the matrix passed as a CRS tuple (`std::tie(n, ptr, col, val)`). It is not a multi-well black-oil simulator.
