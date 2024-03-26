# IMPS: Iterative Methods for Pressure Solvers
### Comparing Iterative Methods on Pressure Solver Performance<br />for an Incompressible Lid-driven Cavity Flow Simulation
<p align="left">
<img src="https://github.com/mw6136/IMPS/assets/144184708/b4cd463d-1b50-4618-9e75-38729b456c86" width="60%">
</p>

## Lid-driven Cavity FLow
<p align="left">
<img src="https://github.com/mw6136/IMPS/assets/144184708/91fe936b-cf6c-4379-8975-3bc750c4e65d" width="75%">
<img src="https://github.com/mw6136/IMPS/assets/144184708/346fae4f-bb4c-46bf-9fcb-cc5455e0a624" width="35%" height="35%">
</p>

## Running the code
### Interactive Jobs
To run a simulation with an interactive slurm job, `run.sh` is a file that provides the ability to run as a simple bash script. At the top of the file, the following variables can be found and set.

```
RUNNAME="example_runname"
PSOLVER=1                # Pressure solver: 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
TOL=0.00001              # Tolerance for the pressure solver
OMEGA=1.8                # Omega for the SOR method

NX1=128                   # Number of grid cells in the x1 direction
NX2=128                   # Number of grid cells in the x2 direction

TMAX=0.004                # Max simulation time
SAVEDT=0.0001           # How often to save the output h5 files
```

`RUNNAME` is the name of the directory where the executable will be compiled to. Once everything is set here, it can simply be run with `./run.sh`.

### Slurm Submission
There are also varying example slurm scripts that can be run by submitting to the queue. For example, one simulation run that was used in the final report was `sor_omega_1.8.sh`, which has a beginning that looks very similar to `run.sh`:

```
#!/bin/bash
#SBATCH --job-name=om1.8        
#SBATCH --nodes=1              
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=8G         
#SBATCH --time=03:00:00 

####################################################################################

RUNNAME="sor_omega_1.8"
PSOLVER=2                # Pressure solver: 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
TOL=0.0000001              # Tolerance for the pressure solver
OMEGA=1.8                # Omega for the SOR method

NX1=64                   # Number of grid cells in the x1 direction
NX2=64                   # Number of grid cells in the x2 direction

TMAX=0.004                # Max simulation time
SAVEDT=0.0001           # How often to save the output h5 files

####################################################################################
```

However now, the sbatch directives can be set for the desired simulation.

## Analysis and Plotting
Within the `analysis` directory, there are various scripts that are used for plotting and making figures.

This is done by gathering directories within the root directory with the naming convention `sor_omega_?.?` and parsing the files within these for the pressure solve time at each iteration. 

A script that is more useful for visulaizing the flow is the `generateplot.py` script. This contains the following CL arguments:

```
-h, --help            show this help message and exit
  -f FIELD, --field FIELD
  -o OUTPUT, --output OUTPUT
```

where field is the field you want to be visualized (has to be one of the ones saved within the IMPS output h5 files), as well as the output the data will be loaded from.

## Solver information (and citation)

The solver implemented in `IMPS` is an explicit Euler solverme as described in [1] using a second-order central-difference approach. Iterative methods for solving the Poisson equation are described in [1] and the lid-driven cavity flow problem and associated boundary conditions are shown in [2]. This work was adapted from the work of S. T. Fush in [3].

[1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery. *Numerical Recipes 3rd Edition: The Art of Scientific
Computing*. Cambridge University Press, USA, 3 edition, 2007.

[2] C. Hirsch. *Numerical computation of internal and external flows: The fundamentals of computational fluid dynamics*.
Elsevier, 2007.

[3] S. T. Fush. VIPS: Varying Incompressible flow Pressure Solver. `https://github.com/trevorfush/vips`, 2023.
