#pragma once

#include <Eigen/Dense>

using namespace Eigen;

namespace DenseForms {
  MatrixXd makeA (const int M, 
                  const int N,
                  const double h,
                  const double viscosity);

  MatrixXd makeLaplacianXBlock (const int M, 
                                const int N, 
                                const double h, 
                                const double viscosity);
  MatrixXd makeLaplacianYBlock (const int M,
                                const int N,
                                const double h,
                                const double viscosity);

  MatrixXd makeDivXBlock (const int M,
                          const int N,
                          const double h);
  MatrixXd makeDivYBlock (const int M,
                          const int N,
                          const double h);

  MatrixXd makeGradXBlock (const int M,
                           const int N,
                           const double h);
  MatrixXd makeGradYBlock (const int M,
                           const int N,
                           const double h);

  MatrixXd makeForcingMatrix (const int M,
                              const int N);

  MatrixXd makeBoundaryMatrix (const int M,
                               const int N,
                               const double h,
                               const double viscosity);

  MatrixXd makeBCLaplacianXBlock (const int M,
                                  const int N,
                                  const double h,
                                  const double viscosity);
  MatrixXd makeBCLaplacianYBlock (const int M,
                                  const int N,
                                  const double h,
                                  const double viscosity);

  MatrixXd makeBCDivXBlock (const int M,
                            const int N,
                            const double h);
  MatrixXd makeBCDivYBlock (const int M,
                            const int N,
                            const double h);
}
