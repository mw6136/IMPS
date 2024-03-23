#include <vector>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "H5Cpp.h"

#ifndef MESH_HPP_
#define MESH_HPP_

struct GridInfo {
  int nx1;
  int nx2;
  int ng;
  double x1min;
  double x1max;
  double x2min;
  double x2max;
};

// Utility functions for saving data
std::string ZeroPadNumber(int num);

std::string getFname(int iter, std::string prefix);

class VArray {
  /*
      Class that contains a 2D vector represented by 1D vector.

      The operators allow the data to be accessed easily by only 
      needing one set of parentheses to specify index location.
  */
  private: int nx1,
  nx2;
  public: std::vector < double > m_data;

  VArray(int i_nx1, int i_nx2, double v0 = 0.0): m_data(i_nx1 * i_nx2, v0) {
    nx1 = i_nx1;
    nx2 = i_nx2;
  }

  double & operator()(int ix, int iy) {
    return m_data[iy + (ix * nx2)];
  }

  double operator()(int ix, int iy) const {
    return m_data[iy + (ix * nx2)];
  }
};

class Mesh {
  /*
      Class containing the simulation grids. These grids can be accessed 
      and updated through the simulation, and are saved using HDF5.

      Also contains empty problem generator function to be specified later
      depending on the problem setup.
  */
  private:

    GridInfo info;

  public:
    // Number of cells in x1, x2 direction, number of ghost zones (should be 1)
    int nx1_,
  nx2_,
  ng_;

  // Grid extremity and spacing values
  double x1max_,
  x1min_,
  x2max_,
  x2min_,
  dx1_,
  dx2_;

  // Velocities at current timestep and next timestep
  VArray un,
  vn,
  unp1,
  vnp1;

  // Pressure at old and new timestep, pressure gradients
  VArray p,
  pnew,
  dpdx,
  dpdy;

  // Arrays for storing part of RHS for pressure solver
  VArray Hx,
  Hy,
  RHS;

  // Wall velocity
  double Uw;

  // Constructor that initializes shapes of arrays
  Mesh(int nx1, int nx2, int ng, double x1max, double x1min, double x2max, double x2min, double iUw): un(nx1 + 2 * ng, nx2 + 2 * ng),
  vn(nx1 + 2 * ng, nx2 + 2 * ng),
  unp1(nx1 + 2 * ng, nx2 + 2 * ng),
  vnp1(nx1 + 2 * ng, nx2 + 2 * ng),
  Hx(nx1 + 2 * ng, nx2 + 2 * ng),
  Hy(nx1 + 2 * ng, nx2 + 2 * ng),
  p(nx1 + 2 * ng, nx2 + 2 * ng),
  dpdx(nx1 + 2 * ng, nx2 + 2 * ng),
  dpdy(nx1 + 2 * ng, nx2 + 2 * ng),
  pnew(nx1 + 2 * ng, nx2 + 2 * ng),
  RHS(nx1 + 2 * ng, nx2 + 2 * ng) {
    info.nx1 = nx1;
    info.nx2 = nx2;
    info.ng = ng;
    info.x1max = x1max;
    info.x1min = x1min;
    info.x2max = x2max;
    info.x2min = x2min;

    nx1_ = nx1;
    nx2_ = nx2;
    ng_ = ng;
    x1max_ = x1max;
    x1min_ = x1min;
    x2max_ = x2max;
    x2min_ = x2min;

    dx1_ = (x1max - x1min) / nx1;
    dx2_ = (x2max - x2min) / nx2;

    Uw = iUw;
  }

  GridInfo getGridInfo() {
    return info;
  }

  // Generates the initial conditions for the mesh
  void ProblemGenerator(double p0, double T0, double R, double gamma);

  // Enforce wall boundary conditions (moving top wall)
  void enforceBCs();

  // Fill the Hx and Hy values (AND RHS WHEN THE TIME COMES)
  void fillRHS(double nu, double rho);

  // Update the velocities with corrected pressure
  void update(double dt, double rho);

  double Jacobi(double rho, double tol);
  double GaussSeidel(double rho, double tol);
  double SOR(double rho, double tol, double omega);

  // Save velocity and pressure fields with grid info
  void save(int iter, double t, std::string fpref, double presstime, double itertime);

};

#endif
