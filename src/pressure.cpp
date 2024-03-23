# include "mesh.hpp"
# include < chrono >
# include <iostream>

  double Mesh::Jacobi(double rho, double tol) {

    double diff = 1e3;
    double maxdiff, iterdiff, newp;
    double calctime;
    int numiters = 0;
    int tot_iters = 40000;
    int i, j;

    auto start = std::chrono::high_resolution_clock::now();

    while (diff > tol) {

      maxdiff = 0.0;

      // For everywhere not on a boundary
      for (i = 2; i <= nx1_ - 1; ++i) {
        for (j = 2; j <= nx2_ - 1; ++j) {

          pnew(i, j) = 0.25 * (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * ((Hx(i + 1, j) - Hx(i - 1, j)) + (Hy(i, j + 1) - Hy(i, j - 1))));

          iterdiff = std::abs((pnew(i, j) - p(i, j)));

          maxdiff = std::max(iterdiff, maxdiff);
        }
      }

      // Left wall 
      i = 1;
      for (j = 2; j <= nx2_ - 1; ++j) {
        pnew(i, j) = 0.3333 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) - Hy(i, j - 1)));

        iterdiff = std::abs((pnew(i, j) - p(i, j)));

        maxdiff = std::max(iterdiff, maxdiff);
      }

      // Right wall
      i = nx1_;
      for (j = 2; j <= nx2_ - 1; ++j) {
        pnew(i, j) = 0.3333 * (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (-(Hx(i, j) + Hx(i - 1, j)) + Hy(i, j + 1) - Hy(i, j - 1)));

        iterdiff = std::abs((pnew(i, j) - p(i, j)));

        maxdiff = std::max(iterdiff, maxdiff);
      }

      // Bottom wall
      j = 1;
      for (i = 2; i <= nx1_ - 1; ++i) {
        pnew(i, j) = 0.3333 * (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) - Hx(i - 1, j) + Hy(i, j + 1) + Hy(i, j)));

        iterdiff = std::abs((pnew(i, j) - p(i, j)));

        maxdiff = std::max(iterdiff, maxdiff);
      }

      // Top wall
      j = nx2_;
      for (i = 2; i <= nx1_ - 1; ++i) {
        pnew(i, j) = 0.3333 * (p(i + 2, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) - Hx(i - 1, j) - (Hy(i, j) + Hy(i, j - 1))));

        iterdiff = std::abs((pnew(i, j) - p(i, j)));

        maxdiff = std::max(iterdiff, maxdiff);
      }

      // Top left corner
      i = 1;
      j = nx2_;

      pnew(i, j) = 0.5 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) - (Hy(i, j) + Hy(i, j - 1))));
      iterdiff = std::abs((pnew(i, j) - p(i, j)));
      maxdiff = std::max(iterdiff, maxdiff);

      // Top right corner
      i = nx1_;
      j = nx2_;

      pnew(i, j) = 0.5 * (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) + 4. * rho * dx1_ * 0.5 * (Hx(i, j) + Hx(i - 1, j) + Hy(i, j) + Hy(i, j - 1)));
      iterdiff = std::abs((pnew(i, j) - p(i, j)));
      maxdiff = std::max(iterdiff, maxdiff);

      // Bottom left corner
      i = 1;
      j = 1;

      pnew(i, j) = 0.5 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) + Hy(i, j)));
      iterdiff = std::abs((pnew(i, j) - p(i, j)));
      maxdiff = std::max(iterdiff, maxdiff);

      // Bottom right corner
      i = nx1_;
      j = 1;

      pnew(i, j) = 0.5 * (p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - p(i + 1, j) + p(i - 1, j) + p(i - 2, j) - 4. * rho * dx1_ * 0.5 * (Hy(i, j + 1) + Hy(i, j) - (Hx(i, j) + Hx(i - 1, j))));
      iterdiff = std::abs((pnew(i, j) - p(i, j)));
      maxdiff = std::max(iterdiff, maxdiff);

      diff = maxdiff;

      enforceBCs();

      // if (numiters % 1000 == 0) {
      //     std::cout << "Pressure correction error : " << diff << std::endl;
      // } 

      numiters += 1;

      // Swap p and pnew values
      for (i = 1; i <= nx1_; ++i) {
        for (j = 1; j <= nx2_; ++j) {
          p(i, j) = pnew(i, j);
        }
      }

    }

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);

    calctime = (double) duration.count();

    return calctime;

  }

double Mesh::GaussSeidel(double rho, double tol) {

  double diff = 1e3;
  double maxdiff, iterdiff, newp;
  double calctime;
  int numiters = 0;
  int tot_iters = 40000;
  int i, j;

  auto start = std::chrono::high_resolution_clock::now();

  while (diff > tol) {

    maxdiff = 0.0;

    // For everywhere not on a boundary
    for (i = 2; i <= nx1_ - 1; ++i) {
      for (j = 2; j <= nx2_ - 1; ++j) {

        newp = 0.25 * (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * ((Hx(i + 1, j) - Hx(i - 1, j)) + (Hy(i, j + 1) - Hy(i, j - 1))));

        iterdiff = std::abs((newp - p(i, j)));

        p(i, j) = newp;

        maxdiff = std::max(iterdiff, maxdiff);
        // if (iterdiff > maxdiff) {
        //     maxdiff = iterdiff;
        // }
      }
    }

    // Left wall 
    i = 1;
    for (j = 2; j <= nx2_ - 1; ++j) {
      newp = 0.3333 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) - Hy(i, j - 1)));

      iterdiff = std::abs((newp - p(i, j)));

      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
      // if (iterdiff > maxdiff) {
      //     maxdiff = iterdiff;
      // }
    }

    // Right wall
    i = nx1_;
    for (j = 2; j <= nx2_ - 1; ++j) {
      newp = 0.3333 * (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (-(Hx(i, j) + Hx(i - 1, j)) + Hy(i, j + 1) - Hy(i, j - 1)));

      iterdiff = std::abs((newp - p(i, j)));

      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
      // if (iterdiff > maxdiff) {
      //     maxdiff = iterdiff;
      // }
    }

    // Bottom wall
    j = 1;
    for (i = 2; i <= nx1_ - 1; ++i) {
      newp = 0.3333 * (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) - Hx(i - 1, j) + Hy(i, j + 1) + Hy(i, j)));

      iterdiff = std::abs((newp - p(i, j)));

      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
      // if (iterdiff > maxdiff) {
      //     maxdiff = iterdiff;
      // }
    }

    // Top wall
    j = nx2_;
    for (i = 2; i <= nx1_ - 1; ++i) {
      newp = 0.3333 * (p(i + 2, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) - Hx(i - 1, j) - (Hy(i, j) + Hy(i, j - 1))));

      iterdiff = std::abs((newp - p(i, j)));

      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
      // if (iterdiff > maxdiff) {
      //     maxdiff = iterdiff;
      // }
    }

    // Top left corner
    i = 1;
    j = nx2_;

    newp = 0.5 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) - (Hy(i, j) + Hy(i, j - 1))));
    iterdiff = std::abs((newp - p(i, j)));
    maxdiff = std::max(iterdiff, maxdiff);
    // if (iterdiff > maxdiff) {
    //     maxdiff = iterdiff;
    // }
    p(i, j) = newp;

    // Top right corner
    i = nx1_;
    j = nx2_;

    newp = 0.5 * (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) + 4. * rho * dx1_ * 0.5 * (Hx(i, j) + Hx(i - 1, j) + Hy(i, j) + Hy(i, j - 1)));
    iterdiff = std::abs((newp - p(i, j)));
    maxdiff = std::max(iterdiff, maxdiff);
    // if (iterdiff > maxdiff) {
    //     maxdiff = iterdiff;
    // }
    p(i, j) = newp;

    // Bottom left corner
    i = 1;
    j = 1;

    newp = 0.5 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) + Hy(i, j)));
    iterdiff = std::abs((newp - p(i, j)));
    maxdiff = std::max(iterdiff, maxdiff);
    // if (iterdiff > maxdiff) {
    //     maxdiff = iterdiff;
    // }
    p(i, j) = newp;

    // Bottom right corner
    i = nx1_;
    j = 1;

    newp = 0.5 * (p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - p(i + 1, j) + p(i - 1, j) + p(i - 2, j) - 4. * rho * dx1_ * 0.5 * (Hy(i, j + 1) + Hy(i, j) - (Hx(i, j) + Hx(i - 1, j))));
    iterdiff = std::abs((newp - p(i, j)));
    maxdiff = std::max(iterdiff, maxdiff);
    // if (iterdiff > maxdiff) {
    //     maxdiff = iterdiff;
    // }
    p(i, j) = newp;

    diff = maxdiff;
    enforceBCs();

    numiters += 1;

  }

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);

  calctime = (double) duration.count();

  return calctime;

}

double Mesh::SOR(double rho, double tol, double omega) {

  double diff = 1e3;
  double maxdiff, iterdiff, tempp, newp, calctime;
  int numiters = 0;
  int i, j;

  auto start = std::chrono::high_resolution_clock::now();

  while (diff > tol) {

    maxdiff = 0.0;

    // For everywhere not on a boundary
    // #pragma omp parallel for collapse(2) private(i,j) reduction(max:maxdiff)
    for (i = 2; i <= nx1_ - 1; ++i) {
      for (j = 2; j <= nx2_ - 1; ++j) {

        tempp = 0.25 * (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * ((Hx(i + 1, j) - Hx(i - 1, j)) + (Hy(i, j + 1) - Hy(i, j - 1))));
        newp = p(i, j) + omega * (tempp - p(i, j));
        iterdiff = std::abs((newp - p(i, j)));
        p(i, j) = newp;
        maxdiff = std::max(iterdiff, maxdiff);
      }
    }

    // Left wall 
    i = 1;
    for (j = 2; j <= nx2_ - 1; ++j) {
      tempp = 0.3333 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) - Hy(i, j - 1)));
      newp = p(i, j) + omega * (tempp - p(i, j));
      iterdiff = std::abs((newp - p(i, j)));
      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
    }

    // Right wall
    i = nx1_;
    for (j = 2; j <= nx2_ - 1; ++j) {
      tempp = 0.3333 * (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (-(Hx(i, j) + Hx(i - 1, j)) + Hy(i, j + 1) - Hy(i, j - 1)));
      newp = p(i, j) + omega * (tempp - p(i, j));
      iterdiff = std::abs((newp - p(i, j)));
      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
    }

    // Bottom wall
    j = 1;
    for (i = 2; i <= nx1_ - 1; ++i) {
      tempp = 0.3333 * (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) - Hx(i - 1, j) + Hy(i, j + 1) + Hy(i, j)));
      newp = p(i, j) + omega * (tempp - p(i, j));
      iterdiff = std::abs((newp - p(i, j)));
      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
    }

    // Top wall
    j = nx2_;
    for (i = 2; i <= nx1_ - 1; ++i) {
      tempp = 0.3333 * (p(i + 2, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) - Hx(i - 1, j) - (Hy(i, j) + Hy(i, j - 1))));
      newp = p(i, j) + omega * (tempp - p(i, j));
      iterdiff = std::abs((newp - p(i, j)));
      p(i, j) = newp;

      maxdiff = std::max(iterdiff, maxdiff);
    }

    // Top left corner
    i = 1;
    j = nx2_;

    tempp = 0.5 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) - (Hy(i, j) + Hy(i, j - 1))));
    newp = p(i, j) + omega * (tempp - p(i, j));
    iterdiff = std::abs((newp - p(i, j)));
    p(i, j) = newp;
    maxdiff = std::max(iterdiff, maxdiff);

    // Top right corner
    i = nx1_;
    j = nx2_;

    tempp = 0.5 * (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) + 4. * rho * dx1_ * 0.5 * (Hx(i, j) + Hx(i - 1, j) + Hy(i, j) + Hy(i, j - 1)));
    newp = p(i, j) + omega * (tempp - p(i, j));
    iterdiff = std::abs((newp - p(i, j)));
    p(i, j) = newp;
    maxdiff = std::max(iterdiff, maxdiff);

    // Bottom left corner
    i = 1;
    j = 1;

    tempp = 0.5 * (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - 4. * rho * dx1_ * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) + Hy(i, j)));
    newp = p(i, j) + omega * (tempp - p(i, j));
    iterdiff = std::abs((newp - p(i, j)));
    p(i, j) = newp;
    maxdiff = std::max(iterdiff, maxdiff);

    // Bottom right corner
    i = nx1_;
    j = 1;

    tempp = 0.5 * (p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - p(i + 1, j) + p(i - 1, j) + p(i - 2, j) - 4. * rho * dx1_ * 0.5 * (Hy(i, j + 1) + Hy(i, j) - (Hx(i, j) + Hx(i - 1, j))));
    newp = p(i, j) + omega * (tempp - p(i, j));
    iterdiff = std::abs((newp - p(i, j)));
    p(i, j) = newp;
    maxdiff = std::max(iterdiff, maxdiff);

    diff = maxdiff;
    enforceBCs();

    // if (numiters % 1000 == 0) {
    //     std::cout << "Pressure correction error : " << diff << std::endl;
    // } 

    numiters += 1;

  }

  // std::cout << "Finished P Iteration" << std::endl;

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);

  calctime = (double) duration.count();

  return calctime;

}