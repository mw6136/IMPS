#include <iostream>
#include <iomanip>
#include <cmath>
#include "omp.h"
#include <chrono>
#include "mesh.hpp"


void Mesh::ProblemGenerator(double p0, double T0, double R, double gamma) {

  double u0 = 0.0;
  double v0 = 0.0;

  double rho0 = p0 / (R * T0);

  for (int i = 0; i < nx1_ + 2 * ng_; ++i) {
    for (int j = 0; j < nx2_ + 2 * ng_; ++j) {

      un(i, j) = u0;
      vn(i, j) = v0;
      p(i, j) = 0.0;

    }
  }

}

int main(int argc, char * argv[]) {

  int nx1, nx2;
  double tmax, savedt, tol, omega;
  std::string solver;

  if (argc == 8) {

    nx1 = std::atoi(argv[1]);
    nx2 = std::atoi(argv[2]);

    if (std::atoi(argv[3]) == 0) {
      solver = "jacobi";
    } else if (std::atoi(argv[3]) == 1) {
      solver = "gs";
    } else if (std::atoi(argv[3]) == 2) {
      solver = "sor";
    } else {
      std::cout << "Please use an implemented pressure solver." << std::endl;
    }

    tol = std::atof(argv[4]);
    omega = std::atof(argv[5]);

    tmax = std::atof(argv[6]);
    savedt = std::atof(argv[7]);

  } else {
    std::cout << "Please ensure to specify the correct variables in the run.sh script." << std::endl;
  }

  // ----------------------------------------------------------------- //
  // Grid dimensions
  // ----------------------------------------------------------------- //
  int ng = 1;
  // ----------------------------------------------------------------- //
  // Fluid properties
  // ----------------------------------------------------------------- //
  double gamma = 1.4;
  double R = 287;

  double Re = 1000; // Reynolds Number
  double T0 = 300; // K
  double Ma = 0.5; // Mach number, just for setting wall velocity
  double P0 = 101330;
  double rho = P0 / (R * T0);

  double L = 1.0 / (Ma / 0.025);
  double x1max = L;
  double x1min = 0.0;
  double x2max = L;
  double x2min = 0.0;

  double a = std::sqrt(gamma * R * T0);
  double Uw = a * Ma;
  double nu = Uw * L / Re;
  // ----------------------------------------------------------------- //
  // Time info
  // ----------------------------------------------------------------- //
  double CFL = 0.3;
  std::string savepref = "IMPS.output.";
  double dt = CFL * ((x1max - x1min) / nx1) / (Uw);
  int printfreq = 10;

  std::cout << tol << std::endl;

  GridInfo info;

  Mesh * pmesh;
  pmesh = new Mesh(nx1, nx2, ng, x1max, x1min, x2max, x2min, Uw);

  info = pmesh -> getGridInfo();

  pmesh -> ProblemGenerator(P0, T0, R, gamma);

  double t = 0;
  int saveiter = 0;
  int iter = 0;
  int savefreq = (int) std::ceil(savedt / dt);
  double presstime, avg_iter_time;
  double itertime = 0.0;
  double pitertime = 0.0;

  while (t <= tmax) {

    auto iterstart = std::chrono::high_resolution_clock::now();
    auto piterstart = std::chrono::high_resolution_clock::now();

    if (iter % printfreq == 0) {
      avg_iter_time = pitertime / printfreq;
      std::cout << "Iteration : " << iter << ", Time: " << t << ", dt = " << dt << ", Avg. Iteration Time (ms) : " << avg_iter_time << std::endl;
      pitertime = 0.0;
    }

    // std::cout << "BCs1" << std::endl;
    // Enforce BCs
    pmesh -> enforceBCs();

    // std::cout << "RHS" << std::endl;
    // Fill Hx, Hy
    pmesh -> fillRHS(nu, rho);

    // std::cout << "BCs2" << std::endl;
    pmesh -> enforceBCs();

    // std::cout << "Press" << std::endl;
    // Do the pressure correction   
    if (solver == "jacobi") {
      presstime = pmesh -> Jacobi(rho, tol);
    } else if (solver == "gs") {
      presstime = pmesh -> GaussSeidel(rho, tol);
    } else if (solver == "sor") {
      presstime = pmesh -> SOR(rho, tol, omega);
    }

    // std::cout << "upd" << std::endl;
    // Update velocities
    pmesh -> update(dt, rho);

    // std::cout << "BCs3" << std::endl;
    pmesh -> enforceBCs();

    auto iterstop = std::chrono::high_resolution_clock::now();
    auto iterduration = std::chrono::duration_cast < std::chrono::milliseconds > (iterstop - iterstart);
    itertime += (double) iterduration.count();

    auto piterstop = std::chrono::high_resolution_clock::now();
    auto piterduration = std::chrono::duration_cast < std::chrono::milliseconds > (iterstop - iterstart);
    pitertime += (double) iterduration.count();

    // Save data
    if (iter % savefreq == 0) {
      itertime = itertime / savefreq;
      std::cout << "Iteration : " << iter << ", Time: " << t << ", dt = " << dt << ", Data saved" << std::endl;
      pmesh -> save(saveiter, t, savepref, presstime, itertime);
      saveiter += 1;
      itertime = 0.0;
    }

    t += dt;
    iter += 1;

  }

  return 0;
}