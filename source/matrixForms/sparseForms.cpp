#include <iostream>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"

#define DEBUG true

using namespace Eigen;

namespace SparseForms {
  void makeStokesMatrix (SparseMatrix<double>& stokesMatrix,
                         const int M,
                         const int N,
                         const double h,
                         const double viscosity) {
    if (DEBUG) cerr << "Creating stokesMatrix." << endl;

    vector<Triplet<double> > tripletList;

    makeLaplacianXBlock (tripletList, 0,                 0,                 M, N, h, viscosity);
    makeLaplacianYBlock (tripletList, M * (N - 1),       M * (N - 1),       M, N, h, viscosity);
    makeGradXBlock      (tripletList, 0,                 2 * M * N - M - N, M, N, h);
    makeGradYBlock      (tripletList, M * (N - 1),       2 * M * N - M - N, M, N, h);
    makeDivXBlock       (tripletList, 2 * M * N - M - N, 0,                 M, N, h);
    makeDivYBlock       (tripletList, 2 * M * N - M - N, M * (N - 1),       M, N, h);
    
    stokesMatrix.setFromTriplets (tripletList.begin(), tripletList.end());
  }

  void makeLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity) {

    if (DEBUG) cerr << "Creating LaplacianXBlock." << endl;

    for (int i = 0; i < M; ++i) {
      for (int x = 0; x < (N - 1); ++x) {
        if ((i == 0) || (i == (M - 1)))
          tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * (N - 1) + x, viscosity * 5 / (h * h)));
        else
          tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * (N - 1) + x, viscosity * 4 / (h * h)));
      }
      for (int x = 0; x < N - 2; ++x) {
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + 1 + x, N0 + i * (N - 1) + x,     -viscosity / (h * h)));
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x,     N0 + i * (N - 1) + 1 + x, -viscosity / (h * h)));
      }
    }

    for (int i = 0; i < (N - 1) * (M - 1); ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,         N0 + N - 1 + i, -viscosity / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + N - 1 + i, N0 + i,         -viscosity / (h * h)));
    }
  }

  void makeLaplacianYBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity) {
    if (DEBUG) cerr << "Creating LaplacianYBlock." << endl;
    for (int i = 0; i < (M - 1); ++i) {
      for (int x = 0; x < N; ++x) {
        if ((x == 0) || (x == (N - 1)))
          tripletList.push_back (Triplet<double> (M0 + i * N + x, N0 + i * N + x, viscosity * 5 / (h * h)));
        else
          tripletList.push_back (Triplet<double> (M0 + i * N + x, N0 + i * N + x, viscosity * 4 / (h * h)));
      }
      for (int x = 0; x < (N - 1); ++x) {
        tripletList.push_back (Triplet<double> (M0 + i * N + x,     N0 + i * N + 1 + x, -viscosity / (h * h)));
        tripletList.push_back (Triplet<double> (M0 + i * N + 1 + x, N0 + i * N + x,     -viscosity / (h * h)));
      }
    }

    for (int i = 0; i < N * (M - 2); ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,     N0 + N + i, -viscosity / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + N + i, N0 + i,     -viscosity / (h * h)));
    }
  }

  void makeGradXBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradXBlock." << endl;
    for (int i = 0; i < M; ++i) {
      for (int x = 0; x < N - 1; ++x) {
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * N + x,     -1 / h));
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * N + 1 + x,  1 / h));
      }
    }
  }

  void makeGradYBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradYBlock." << endl;
    for (int i = 0; i < (M - 1) * N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i, N0 + i,     -1 / h));
      tripletList.push_back (Triplet<double> (M0 + i, N0 + N + i,  1 / h));
    }
  }

  void makeDivXBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h) {
    if (DEBUG) cerr << "Creating DivXBlock." << endl;
    for (int i = 0; i < M; ++i) {
      for (int x = 0; x < (N - 1); ++x) {
        tripletList.push_back (Triplet<double> (M0 + i * N + x,     N0 + i * (N - 1) + x,  1 / h));
        tripletList.push_back (Triplet<double> (M0 + i * N + 1 + x, N0 + i * (N - 1) + x, -1 / h));
      }
    }
  }

  void makeDivYBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h) {
    if (DEBUG) cerr << "Creating DivYBlock." << endl;
    for (int i = 0; i < (M - 1) * N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,     N0 + i,  1 / h));
      tripletList.push_back (Triplet<double> (M0 + i + N, N0 + i, -1 / h));
    }
  }

  void makeForcingMatrix (SparseMatrix<double>& forcingMatrix,
                          const int M,
                          const int N) {
    if (DEBUG) cerr << "Creating ForcingMatrix." << endl;
    vector<Triplet<double> > tripletList;

    for (int i = 0; i < 2 * M * N - M - N; ++i)
      tripletList.push_back (Triplet<double> (i, i, 1));
  
    forcingMatrix.setFromTriplets (tripletList.begin (), tripletList.end ());
  }

  void makeBoundaryMatrix (SparseMatrix<double>& boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double viscosity) {
    if (DEBUG) cerr << "Creating BoundaryMatrix." << endl;
    vector<Triplet<double> > tripletList;

    makeBCLaplacianXBlock (tripletList, 0,                 0,     M, N, h, viscosity);
    makeBCLaplacianYBlock (tripletList, M * (N - 1),       2 * M, M, N, h, viscosity);
    makeBCDivXBlock       (tripletList, 2 * M * N - M - N, 0,     M, N, h);
    makeBCDivYBlock       (tripletList, 2 * M * N - M - N, 2 * M, M, N, h);

    boundaryMatrix.setFromTriplets (tripletList.begin(), tripletList.end());
  }

  void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                              const int M0,
                              const int N0,
                              const int M,
                              const int N,
                              const double h,
                              const double viscosity) {
    if (DEBUG) cerr << "Creating BCLaplacianXBlock." << endl;
    for (int i = 0; i < M; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i * (N - 1),           N0 + i * 2,     viscosity / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + (i + 1) * (N - 1) - 1, N0 + i * 2 + 1, viscosity / (h * h)));
    }
  }

  void makeBCLaplacianYBlock (vector<Triplet<double> >& tripletList,
                              const int M0,
                              const int N0,
                              const int M,
                              const int N,
                              const double h,
                              const double viscosity) {
    if (DEBUG) cerr << "Creating BCLaplacianYBlock." << endl;
    for (int i = 0; i < N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,               N0 + i,     viscosity / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + (M - 2) * N + i, N0 + N + i, viscosity / (h * h)));
    }
  }

  void makeBCDivXBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivXBlock." << endl;
    for (int i = 0; i < M; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i * N,           N0 + i * 2,      1 / h));
      tripletList.push_back (Triplet<double> (M0 + (i + 1) * N - 1, N0 + i * 2 + 1, -1 / h));
    }
  }

  void makeBCDivYBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivYBlock." << endl;
    for (int i = 0; i < N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,               N0 + i,      1 / h));
      tripletList.push_back (Triplet<double> (M0 + (M - 1) * N + i, N0 + N + i, -1 / h));
    }
  }
}
