#include <iostream>

#include <Eigen/Dense>

#include "matrixForms/denseForms.h"

#define DEBUG true

using namespace Eigen;
using namespace std;

namespace DenseForms {
  void makeA (Ref<MatrixXd> A,
              const int M,
              const int N,
              const double h,
              const double viscosity) {
    if (DEBUG) cerr << "Creating A." << endl;
    A = MatrixXd::Zero (3 * M * N - M - N, 3 * M * N - M - N);

    makeLaplacianXBlock (A.block (0,                 0,                 M * (N - 1), M * (N - 1)), M, N, h, viscosity);
    makeLaplacianYBlock (A.block (M * (N - 1),       M * (N - 1),       (M - 1) * N, (M - 1) * N), M, N, h, viscosity);
    makeGradXBlock      (A.block (0,                 2 * M * N - M - N, M * (N - 1), M * N),       M, N, h);
    makeGradYBlock      (A.block (M * (N - 1),       2 * M * N - M - N, (M - 1) * N, M * N),       M, N, h);
    makeDivXBlock       (A.block (2 * M * N - M - N, 0,                  M * N,      M * (N - 1)), M, N, h);
    makeDivYBlock       (A.block (2 * M * N - M - N, M * (N - 1),        M * N,      (M - 1) * N), M, N, h);
  }

  void makeLaplacianXBlock (Ref<MatrixXd> laplacian,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity) {
    if (DEBUG) cerr << "Creating LaplacianXBlock." << endl;
    laplacian               = MatrixXd::Zero (M * (N - 1), M * (N - 1));
    MatrixXd laplacianBlock = MatrixXd::Zero (N - 1,       N - 1);

    laplacianBlock.diagonal ()   = VectorXd::Constant (N - 1, viscosity * 4 / (h * h));
    laplacianBlock.diagonal (-1) = laplacianBlock.diagonal (1) = VectorXd::Constant (N - 2, -viscosity / (h * h));

    for (int i = 0; i < M; ++i) {
      laplacianBlock.diagonal() = VectorXd::Constant (N - 1, viscosity * 4 / (h * h));
      if ((i == 0) || (i == M - 1))
        laplacianBlock.diagonal () = VectorXd::Constant (N - 1, viscosity * 5 / (h * h));

      laplacian.block (i * (N - 1), i * (N - 1), (N - 1), (N - 1)) = laplacianBlock;
    }

    laplacian.diagonal (N - 1)    = VectorXd::Constant ((N - 1) * (M - 1), -viscosity / (h * h));
    laplacian.diagonal (-(N - 1)) = VectorXd::Constant ((N - 1) * (M - 1), -viscosity / (h * h));

  }

  void makeLaplacianYBlock (Ref<MatrixXd> laplacian,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity) {
    if (DEBUG) cerr << "Creating LaplacianYBlock." << endl;
    laplacian               = MatrixXd::Zero ((M - 1) * N, (M - 1) * N);
    MatrixXd laplacianBlock = MatrixXd::Zero (N,           N);

    laplacianBlock.diagonal () = VectorXd::Constant (N, viscosity * 4 / (h * h));
    laplacianBlock.diagonal (-1) = laplacianBlock.diagonal (1) = VectorXd::Constant (N - 1, -viscosity / (h * h));
    laplacianBlock (0, 0) = laplacianBlock (N - 1, N - 1) = viscosity * 5 / (h * h);

    for (int i = 0; i < M - 1; ++i)
      laplacian.block (i * N, i * N, N, N) = laplacianBlock;

    laplacian.diagonal (N) = laplacian.diagonal (-N) = VectorXd::Constant ((M - 2) * N, -viscosity / (h * h));
  }

  void makeGradXBlock (Ref<MatrixXd> grad,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradXBLock." << endl;
    MatrixXd gradBlock = MatrixXd::Zero (N - 1, N);

    gradBlock.diagonal()  = VectorXd::Constant (N - 1, -1 / h);
    gradBlock.diagonal(1) = VectorXd::Constant (N - 1,  1 / h);

    for (int i = 0; i < M; ++i)
      grad.block (i * (N - 1), i * N, N - 1, N) = gradBlock;

  }

  void makeGradYBlock (Ref<MatrixXd> grad,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradYBlock." << endl;

    grad.diagonal()  = VectorXd::Constant ((M - 1) * N, -1 / h);
    grad.diagonal(N) = VectorXd::Constant ((M - 1) * N,  1 / h); 
  }

  void makeDivXBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h) {
    if (DEBUG) cerr << "Creating DivXBLock." << endl;
    div               = MatrixXd::Zero (M * N, M * (N - 1));
    MatrixXd divBlock = MatrixXd::Zero (N,     N - 1);

    divBlock.diagonal ()   = VectorXd::Constant (N - 1,  1 / h);
    divBlock.diagonal (-1) = VectorXd::Constant (N - 1, -1 / h);

    for (int i = 0; i < M; ++i)
      div.block (i * N, i * (N - 1), N, N - 1) = divBlock;
  }

  void makeDivYBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h) {
    if (DEBUG) cerr << "Creating DivYBlock." << endl;
    div = MatrixXd::Zero (M * N, (M - 1) * N);

    div.diagonal ()   = VectorXd::Constant ((M - 1) * N,  1 / h);
    div.diagonal (-N) = VectorXd::Constant ((M - 1) * N, -1 / h);
  }

  void makeForcingMatrix (Ref<MatrixXd> forcingMatrix,
                          const int M,
                          const int N) {
    if (DEBUG) cerr << "Creating ForcingMatrix." << endl;
    forcingMatrix = MatrixXd::Zero (3 * M * N - M - N, 2 * M * N - M - N);

    forcingMatrix.block (0, 0, 2 * M * N - M - N, 2 * M * N - M - N) = MatrixXd::Identity (2 * M * N - M - N, 2 * M * N - M - N);
  }

  void makeBoundaryMatrix (Ref<MatrixXd> boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double viscosity) {
    if (DEBUG) cerr << "Creating BoundaryMatrix." << endl;
    boundaryMatrix = MatrixXd::Zero (3 * M * N - M - N, 2 * M + 2 * N);

    makeBCLaplacianXBlock (boundaryMatrix.block (0,                 0,     M * (N - 1), 2 * M), M, N, h, viscosity);
    makeBCLaplacianYBlock (boundaryMatrix.block (M * (N - 1),       2 * M, (M - 1) * N, 2 * N), M, N, h, viscosity);
    makeBCDivXBlock       (boundaryMatrix.block (2 * M * N - M - N, 0,     M * N,       2 * M), M, N, h);
    makeBCDivYBlock       (boundaryMatrix.block (2 * M * N - M - N, 2 * M, M * N,       2 * N), M, N, h);
  }

  void makeBCLaplacianXBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,
                              const double viscosity) {
    if (DEBUG) cerr << "Creating BCLaplacianXBlock." << endl;
    laplacianBC               = MatrixXd::Zero (M * (N - 1), 2 * M);
    MatrixXd laplacianBCBlock = MatrixXd::Zero (N - 1,       2);

    laplacianBCBlock (0, 0) = laplacianBCBlock (N - 2, 1) = viscosity / (h * h);

    for (int i = 0; i < M; ++i)
      laplacianBC.block (i * (N - 1), i * 2, N - 1, 2) = laplacianBCBlock;
  }

  void makeBCLaplacianYBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,
                              const double viscosity) {
    if (DEBUG) cerr << "Creating BCLaplacianYBlock." << endl;
    laplacianBC               = MatrixXd::Zero ((M - 1) * N, 2 * N);
    MatrixXd laplacianBCBlock = MatrixXd::Zero (N,           N);

    laplacianBCBlock.diagonal() = VectorXd::Constant (N, viscosity / (h * h));

    laplacianBC.block (0,           0, N, N) = laplacianBCBlock;
    laplacianBC.block ((M - 2) * N, N, N, N) = laplacianBCBlock;
  }

  void makeBCDivXBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivXBlock." << endl;
    divBC               = MatrixXd::Zero (M * N, 2 * M);
    MatrixXd divBCBlock = MatrixXd::Zero (N,     2);

    divBCBlock (0,     0) =  1 / h; 
    divBCBlock (N - 1, 1) = -1 / h;

    for (int i = 0; i < M; ++i)
      divBC.block (i * N, i * 2, N, 2) = divBCBlock;
  }

  void makeBCDivYBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivYBlock." << endl;
    divBC = MatrixXd::Zero (M * N, 2 * N);

    divBC.block (0,           0, N, N) =  1 / h * MatrixXd::Identity (N, N);
    divBC.block ((M - 1) * N, N, N, N) = -1 / h * MatrixXd::Identity (N, N);
  }
}

  
