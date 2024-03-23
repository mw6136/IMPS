#include <iostream>
#include "omp.h"
#include "mesh.hpp"

std::string ZeroPadNumber(int num) {
  std::ostringstream ss;
  ss << std::setw(4) << std::setfill('0') << num;
  std::string result = ss.str();
  if (result.length() > 7) {
    result.erase(0, result.length() - 7);
  }
  return result;
}

std::string getFname(int iter, std::string prefix) {
  std::string fname;
  std::stringstream ss;
  std::stringstream ssfname;
  std::string newnum = ZeroPadNumber(iter);

  ss << prefix << newnum << ".h5";

  ss >> fname;

  return fname;
}

void Mesh::enforceBCs() {

  int nx1t = nx1_ + 2 * ng_;
  int nx2t = nx2_ + 2 * ng_;

  // Left wall
  for (int j = 1; j <= nx2_; ++j) {

    // No slip wall condition
    un(0, j) = -un(1, j);
    vn(0, j) = -vn(1, j);

    // Constant pressure gradient
    p(0, j) = p(1, j);
    pnew(0, j) = pnew(1, j);

    Hx(0, j) = -Hx(1, j);

  }

  // Right wall 
  for (int j = 1; j <= nx2_; ++j) {

    un(nx1t - 1, j) = -un(nx1t - 2, j);
    vn(nx1t - 1, j) = -vn(nx1t - 2, j);

    p(nx1t - 1, j) = p(nx1t - 2, j);
    pnew(nx1t - 1, j) = pnew(nx1t - 2, j);

    Hx(nx1t - 1, j) = -Hx(nx1t - 2, j);

  }

  // Bottom wall
  for (int i = 0; i <= nx1_ + 1; ++i) {

    un(i, 0) = -un(i, 1);
    vn(i, 0) = -vn(i, 1);

    p(i, 0) = p(i, 1);
    pnew(i, 0) = pnew(i, 1);

    Hy(i, 0) = -Hy(i, 1);

  }

  // Top wall
  for (int i = 0; i <= nx1_ + 1; ++i) {

    un(i, nx2t - 1) = 2 * Uw - un(i, nx2t - 2);
    vn(i, nx2t - 1) = -vn(i, nx2t - 2);

    p(i, nx2t - 1) = p(i, nx2t - 2);
    pnew(i, nx2t - 1) = pnew(i, nx2t - 2);

    Hy(i, nx2t - 1) = -Hy(i, nx2t - 2);

  }
}

void Mesh::fillRHS(double nu, double rho) {

  // Fill Hx and Hy
  double uipoh, uimoh, ujpoh, ujmoh;
  double vipoh, vimoh, vjpoh, vjmoh;
  int i, j;

  // #pragma omp parallel for collapse(2)
  for (i = 1; i <= nx1_; ++i) {
    for (j = 1; j <= nx2_; ++j) {

      uipoh = (un(i + 1, j) + un(i, j)) / 2;
      uimoh = (un(i, j) + un(i - 1, j)) / 2;
      ujpoh = (un(i, j + 1) + un(i, j)) / 2;
      ujmoh = (un(i, j) + un(i, j - 1)) / 2;

      vipoh = (vn(i + 1, j) + vn(i, j)) / 2;
      vimoh = (vn(i, j) + vn(i - 1, j)) / 2;
      vjpoh = (vn(i, j + 1) + vn(i, j)) / 2;
      vjmoh = (vn(i, j) + vn(i, j - 1)) / 2;

      Hx(i, j) = -((uipoh * uipoh - uimoh * uimoh) / dx1_) - ((ujpoh * vjpoh - ujmoh * vjmoh) / dx2_) +
        nu * ((un(i + 1, j) - 2 * un(i, j) + un(i - 1, j)) / (dx1_ * dx1_) + (un(i, j + 1) - 2 * un(i, j) + un(i, j - 1)) / (dx2_ * dx2_));

      Hy(i, j) = -((vjpoh * vjpoh - vjmoh * vjmoh) / dx2_) - ((uipoh * vipoh - uimoh * vimoh) / dx1_) +
        nu * ((vn(i + 1, j) - 2 * vn(i, j) + vn(i - 1, j)) / (dx1_ * dx1_) + (vn(i, j + 1) - 2 * vn(i, j) + vn(i, j - 1)) / (dx2_ * dx2_));

    }
  }

}

void Mesh::update(double dt, double rho) {

  double pipoh, pimoh, pjpoh, pjmoh;

  // #pragma omp parallel for collapse(2)
  for (int i = 1; i <= nx1_; ++i) {
    for (int j = 1; j <= nx2_; ++j) {

      pipoh = (p(i + 1, j) + p(i, j)) / 2;
      pimoh = (p(i, j) + p(i - 1, j)) / 2;
      pjpoh = (p(i, j + 1) + p(i, j)) / 2;
      pjmoh = (p(i, j) + p(i, j - 1)) / 2;

      un(i, j) += dt * Hx(i, j) - dt / rho * ((pipoh - pimoh) / dx1_);
      vn(i, j) += dt * Hy(i, j) - dt / rho * ((pjpoh - pjmoh) / dx2_);

    }
  }

}

void Mesh::save(int iter, double t, std::string fpref, double presstime, double itertime) {

  std::string fname = getFname(iter, fpref);

  int nx1 = nx1_ + 2 * ng_;
  int nx2 = nx2_ + 2 * ng_;

  std::vector < int > griddim = {
    nx1,
    nx2,
    ng_
  };
  std::vector < double > gridinfo = {
    x1min_,
    x1max_,
    x2min_,
    x2max_
  };
  std::vector < double > pressinfo = {
    t,
    itertime,
    presstime
  };

  H5::H5File file(fname.c_str(), H5F_ACC_TRUNC);

  hsize_t dims[1];
  dims[0] = nx2 * nx1;

  H5::DataSpace dataspace(1, dims);

  H5::DataSet XVEL = file.createDataSet("VelX1", H5::PredType::NATIVE_DOUBLE, dataspace);
  XVEL.write(un.m_data.data(), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet YVEL = file.createDataSet("VelX2", H5::PredType::NATIVE_DOUBLE, dataspace);
  YVEL.write(vn.m_data.data(), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet PRES = file.createDataSet("Press", H5::PredType::NATIVE_DOUBLE, dataspace);
  PRES.write(p.m_data.data(), H5::PredType::NATIVE_DOUBLE);

  hsize_t dims2[1];
  dims2[0] = 3;
  H5::DataSpace griddim_space(1, dims2);
  H5::DataSet GDIM = file.createDataSet("Dim", H5::PredType::NATIVE_INT, griddim_space);
  GDIM.write(griddim.data(), H5::PredType::NATIVE_INT);

  hsize_t dims3[1];
  dims3[0] = 4;
  H5::DataSpace gridinfo_space(1, dims3);
  H5::DataSet GINF = file.createDataSet("Info", H5::PredType::NATIVE_DOUBLE, gridinfo_space);
  GINF.write(gridinfo.data(), H5::PredType::NATIVE_DOUBLE);

  hsize_t dims4[1];
  dims4[0] = 3;
  H5::DataSpace pressinfo_space(1, dims4);
  H5::DataSet PINFO = file.createDataSet("TimeInfo", H5::PredType::NATIVE_DOUBLE, pressinfo_space);
  PINFO.write(pressinfo.data(), H5::PredType::NATIVE_DOUBLE);

}