#!/bin/bash
#SBATCH --job-name=om1.5        
#SBATCH --nodes=1              
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=8G         
#SBATCH --time=03:00:00 

####################################################################################

RUNNAME="sor_omega_1.5"
PSOLVER=2                # Pressure solver: 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
TOL=0.0000001              # Tolerance for the pressure solver
OMEGA=1.5                # Omega for the SOR method

NX1=64                   # Number of grid cells in the x1 direction
NX2=64                   # Number of grid cells in the x2 direction

TMAX=0.004                # Max simulation time
SAVEDT=0.0001           # How often to save the output h5 files

####################################################################################

export HDF5_DISABLE_VERSION_CHECK=2

module load hdf5/gcc/1.10.6

if [ ! -d ${RUNNAME} ]
then
    mkdir ${RUNNAME}
fi 

cd ${RUNNAME}

echo "Compiling mesh.cpp and pressure.cpp..."
rm ../src/*.o

g++ -O3 -fopenmp -c -o ../src/mesh.o ../src/mesh.cpp
g++ -O3 -c -o ../src/pressure.o ../src/pressure.cpp

echo "Compiling main.cpp..."
g++ -O3 -fopenmp -o IMPS ../src/main.cpp ../src/mesh.o ../src/pressure.o -lhdf5 -lhdf5_cpp

echo "Running IMPS..."
time ./IMPS ${NX1} ${NX2} ${PSOLVER} ${TOL} ${OMEGA} ${TMAX} ${SAVEDT} > ${RUNNAME}.info

echo "Simulation complete!"