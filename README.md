# IMPS: Iterative Methods for Pressure Solvers
### Comparing Iterative Methods on Pressure Solver Performanc e<br />for an Incompressible Lid-driven Cavity Flow Simulation
<p align="center">
<img src="https://github.com/mw6136/IMPS/assets/144184708/e160933c-7c15-4d0b-a7f8-ac9eb382a8c5" width="35%" height="35%">
</p>

## Lid-driven Cavity FLow
<p align="center">
<img src="https://github.com/mw6136/IMPS/assets/144184708/91fe936b-cf6c-4379-8975-3bc750c4e65d" width="70%" height="70%">
<img src="https://github.com/mw6136/IMPS/assets/144184708/346fae4f-bb4c-46bf-9fcb-cc5455e0a624" width="35%" height="35%">
</p>

## Solver information (and citation)

The solver implemented in `IMPS` is an explicit Euler solverme as described in [1] using a second-order central-difference approach. Iterative methods for solving the Poisson equation are described in [1] and the lid-driven cavity flow problem and associated boundary conditions are shown in [1]. This work was adapted from the work of S. T. Fush in [3].

[1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery. Numerical Recipes 3rd Edition: The Art of Scientific
Computing. Cambridge University Press, USA, 3 edition, 2007.

[2] C. Hirsch. Numerical computation of internal and external flows: The fundamentals of computational fluid dynamics.
Elsevier, 2007.

[3] S. T. Fush. VIPS: Varying Incompressible flow Pressure Solver. `https://github.com/trevorfush/vips`, 2023.
