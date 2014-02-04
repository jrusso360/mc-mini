#pragma once

#include <Eigen/Sparse>

namespace SparseForms {
  SparseMatrix <double> makeA (const int M, 
                               const int N, 
                               const double h, 
                               const double viscosity);

  void makeLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity);

  void makeLaplacianYBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity);

  void makeGradXBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h);

  void makeGradYBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h);

  void makeDivXBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h);

  void makeDivYBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h);

  SparseMatrix<double> makeForcingMatrix (const int M,
                                          const int N);

  SparseMatrix<double> makeBoundaryMatrix (const int N,
                                           const int M,
                                           const double h,
                                           const double viscosity);

  void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                              const int M0,
                              const int N0,
                              const int M,
                              const int N,
                              const double h,
                              const double viscosity);

  void makeBCLaplacianYBlock (vector<Triplet>double> >& tripletList,
                              const int M0,
                              const int N0,
                              const int M,
                              const int N,
                              const double h,
                              const double viscosity);

  void makeBCDivXBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h);

  void makeBCDivYBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h);
}
