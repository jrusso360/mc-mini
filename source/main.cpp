#include <iostream>
#include <iomanip>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "paramParse/parser.h"

#include "matrixForms/denseForms.h"
#include "matrixForms/sparseForms.h"

using namespace Eigen;
using namespace std;
/*
SparseMatrix<double> makeA (const int, const int, const double, const double);
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
void makeGradXBlock      (vector<Triplet<double> >& tripletList,
                          const int M0,
                          const int N0,
                          const int M,
                          const int N,
                          const double h);
void makeGradYBlock      (vector<Triplet<double> >& tripletList,
                          const int M0,
                          const int N0,
                          const int M,
                          const int N,
                          const double h);
void makeDivXBlock       (vector<Triplet<double> >& tripletList,
                          const int M0,
                          const int N0,
                          const int M,
                          const int N,
                          const double h);
void makeDivYBlock       (vector<Triplet<double> >& tripletList,
                          const int M0,
                          const int N0,
                          const int M,
                          const int N,
                          const double h);


SparseMatrix<double> makeForcingMatrix (const int, const int);

SparseMatrix<double> makeBoundaryMatrix (const int, const int, const double, const double);
void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity);
void makeBCLaplacianYBlock (vector<Triplet<double> >& tripletList,
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
*/
int main(int argc, char ** argv) {

  int M, N;
  double viscosity, y_extent, h;

  string paramFile = argv[1];

  paramParser pp(paramFile);

  pp.getParamInt ("M", M); // Number of rows
  pp.getParamInt ("N", N); // Number of columns

  pp.getParamDouble ("viscosity", viscosity);

  pp.getParamDouble ("y_extent", y_extent);
  h = y_extent / N;

  // Set up locations of data in memory
  double * solnData     = new double[3 * M * N - M - N];
  double * boundaryData = new double[2 * M + 2 * N];
  double * forcingData  = new double[2 * M * N - M - N];
  double * rhsData      = new double[3 * M * N - M - N];
  // U data pointers
  double * uSolnData     = solnData;
  double * uBoundaryData = boundaryData;
  double * uForcingData  = forcingData;
  double * uRhsData      = rhsData;
  // V data pointers
  double * vSolnData     = solnData     + (M * (N - 1));
  double * vBoundaryData = boundaryData + 2 * M;
  double * vForcingData  = forcingData  + (M * (N - 1));
  double * vRhsData      = rhsData      + (M * (N - 1));
  // P data pointers
  double * pSolnData = solnData + (2 * M * N - M - N);
  double * pRhsData  = rhsData  + (N * (2 * N - 2));
  
  // Set up vector maps to data
  Map<VectorXd> solnVect     (solnData,     3 * M * N - M - N);
  Map<VectorXd> boundaryVect (boundaryData, 2 * M + 2 * N);
  Map<VectorXd> forcingVect  (forcingData,  2 * M * N - M - N);
  Map<VectorXd> rhsVect      (rhsData,      3 * M * N - M - N);

  // U vector maps
  Map<VectorXd> uSolnVect     (uSolnData,     M * (N - 1));
  Map<VectorXd> uBoundaryVect (uBoundaryData, 2 * M);
  Map<VectorXd> uForcingVect  (uForcingData,  M * (N - 1));
  Map<VectorXd> uVect         (uRhsData,      M * (N - 1));
  // U matrix maps
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uSolnMatrix     (uSolnData,     M, N - 1);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uBoundaryMatrix (uBoundaryData, M, 2);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uForcingMatrix  (uForcingData,  M, N - 1);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uRhsMatrix      (uRhsData,      M, N - 1);

  //V vector maps
  Map<VectorXd> vSolnVect     (vSolnData,     N * (M - 1));
  Map<VectorXd> vBoundaryVect (vBoundaryData, 2 * N);
  Map<VectorXd> vForcingVect  (vForcingData,  N * (M - 1));
  Map<VectorXd> vRhsVect      (vRhsData,      N * (M - 1));
  // V matrix maps
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vSolnMatrix     (vSolnData,     M - 1, N);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vBoundaryMatrix (vBoundaryData, 2,     N);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vForcingMatrix  (vForcingData,  M - 1, N);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vRhsMatrix      (vRhsData,      M - 1, N);

  // P vector maps
  Map<VectorXd> pSolnVect (pSolnData, M * N);
  Map<VectorXd> pRhsVect  (pRhsData,  M * N);
  // P matrix maps
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > pSolnMatrix (pSolnData, M, N);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > pRhsMatrix  (pRhsData,  M, N);

  MatrixXd denseA (3 * M * N - M - N, 3 * M * N - M - N);  
  DenseForms::makeA (denseA, M, N, h, viscosity);
  MatrixXd denseForcingMatrix (3 * M * N - M - N, 2 * M * N - M - N);
  DenseForms::makeForcingMatrix  (denseForcingMatrix, M, N); 
  MatrixXd denseBoundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N); 
  DenseForms::makeBoundaryMatrix (denseBoundaryMatrix, M, N, h, viscosity); 

  SparseMatrix<double> sparseA (3 * M * N - M - N, 3 * M * N - M - N);
  SparseForms::makeA (sparseA, M, N, h, viscosity);
  SparseMatrix<double> sparseForcingMatrix (3 * M * N - M - N, 2 * M * N - M - N);
  SparseForms::makeForcingMatrix (sparseForcingMatrix, M, N);
  
  MatrixXd differenceA = denseA + (MatrixXd::Identity (3 * M * N - M - N, 3 * M * N - M - N) * -1 * sparseA);
  MatrixXd differenceForcing = denseForcingMatrix + (MatrixXd::Identity (3 * M * N - M - N, 3 * M * N - M - N) * - 1 * sparseForcingMatrix);
  MatrixXd differenceBoundary = denseBoundaryMatrix;
  cout << differenceA << endl << endl
       << differenceForcing << endl << endl
       << differenceBoundary << endl << endl;
 
  /*
  for (int i = 0; i < M; ++i) 
    for (int j = 0; j < N - 1; ++j) 
      uForcingMatrix (i, j) = 3 * cos ((j + 1) * h) * sin ((i + 0.5) * h);
  for (int i = 0; i < M - 1; ++i)
    for (int j = 0; j < N; ++j) 
      vForcingMatrix (i, j) = -1 * sin ((j + 0.5) * h) * cos ((i + 1) * h);
  
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < 2; ++j)
      uBoundaryMatrix (i, j) = cos (j * N * h) * sin ((i + 0.5) * h);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < N; ++j)
      vBoundaryMatrix (i, j) = -1 * sin ((j + 0.5) * h) * cos (i * M * h);
 
  rhsVect = forcingMatrix * forcingVect + boundaryMatrix * boundaryVect;
  
  solnVect = A.householderQr().solve(rhsVect);

  cout << uSolnMatrix.colwise().reverse() << endl << endl
       << vSolnMatrix.colwise().reverse() << endl << endl
       << pSolnMatrix.colwise().reverse() << endl << endl;

  MatrixXd analyticU;
  analyticU.resizeLike(uSolnMatrix);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N - 1; ++j)
      analyticU (i, j) = cos ((j + 1) * h) * sin ((i + 0.5) * h);

  MatrixXd analyticV;
  analyticV.resizeLike(vSolnMatrix);
  for (int i = 0; i < M - 1; ++i)
    for (int j = 0; j < N; ++j)
      analyticV (i, j) = - sin ((j + 0.5) * h) * cos ((i + 1) * h);

  cout << sqrt((uSolnMatrix - analyticU).squaredNorm() * h * h) << "\t" << sqrt((vSolnMatrix - analyticV).squaredNorm() * h * h) << endl;
*/
  return 0;
}
/*
SparseMatrix<double> makeA (const int M, 
                            const int N, 
                            const double h, 
                            const double viscosity) {
  SparseMatrix<double> A (3 * M * N - M - N, 3 * M * N - M -N);
  A.reserve (21 * M * N - 7 * M - 7 * N);
  
  vector<Triplet<double> > tripletList;
  makeLaplacianXBlock (tripletList, 0,                 0,                 M, N, h, viscosity);
  makeLaplacianYBlock (tripletList, M * (N - 1),       M * (N - 1),       M, N, h, viscosity);
  makeGradXBlock      (tripletList, 0,                 2 * M * N - M - N, M, N, h);
  makeGradYBlock      (tripletList, M * (N - 1),       2 * M * N - M - N, M, N, h);
  makeDivXBlock       (tripletList, 2 * M * N - M - N, 0,                 M, N, h);
  makeDivYBlock       (tripletList, 2 * M * N - M - N, M * (N - 1),       M, N, h);

  A.setFromTriplets (tripletList.begin(), tripletList.end());

  return A;
}

void makeLaplacianXBlock (vector<Triplet<double> >& tripletList, 
                          const int M0, 
                          const int N0,
                          const int M,
                          const int N,
                          const double h,
                          const double viscosity) {

  for (int i = 0; i < M; ++i) {
    for (int x = 0; x < N - 1; ++x) {
      if ((i == 0) || (i == M - 1))
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * (N - 1) + x, viscosity * 5 / (h * h)));
      else
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * (N - 1) + x, viscosity * 4 / (h * h)));
    }
    for (int x = 0; x < N - 2; ++x) {
      tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + 1 + x, N0 + i * (N - 1) + x, -viscosity / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * (N - 1) + 1 + x, -viscosity / (h * h)));
    }
  }

  for (int i = 0; i < (N - 1) * (M - 1); ++i) {
    tripletList.push_back (Triplet<double> (M0 + i, N0 + N - 1 + i, -viscosity / (h * h)));
    tripletList.push_back (Triplet<double> (M0 + N - 1 + i, N0 + i, -viscosity / (h * h)));
  }
}

void makeLaplacianYBlock (vector<Triplet<double> >& tripletList,
                          const int M0,
                          const int N0,
                          const int M,
                          const int N,
                          const double h,
                          const double viscosity) {

  for (int i = 0; i < M - 1; ++i) {
    for (int x = 0; x < N; ++x) {
      if ((x == 0) || (x == (N - 1)))
        tripletList.push_back (Triplet<double> (M0 + i * N + x, N0 + i * N + x, viscosity * 5 / (h * h)));
      else
        tripletList.push_back (Triplet<double> (M0 + i * N + x, N0 + i * N + x, viscosity * 4 / (h * h)));
    }
    for (int x = 0; x < (N - 1); ++x) {
      tripletList.push_back (Triplet<double> (M0 + i * N + x, N0 + i * N + x + 1, -viscosity / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + i * N + x + 1, N0 + i * N + x, -viscosity / (h * h)));
    }
  }

  for (int i = 0; i < N * (M - 2); ++i) {
    tripletList.push_back (Triplet<double> (M0 + i, N0 + N + i, -viscosity / (h * h)));
    tripletList.push_back (Triplet<double> (M0 + N + i, N0 + i, -viscosity / (h * h)));
  }
}

void makeGradXBlock (vector<Triplet<double> >& tripletList,
                     const int M0,
                     const int N0,
                     const int M,
                     const int N,
                     const double h) {

  for (int i = 0; i < N; ++i) {
    for (int x = 0; x < M - 1; ++x) {
      tripletList.push_back (Triplet<double> (M0 + i * (M - 1) + x, N0 + i * M + x, -1 / h));
      tripletList.push_back (Triplet<double> (M0 + i * (M - 1) + x, N0 + i * M + x + 1, 1 / h));
    }
  }
}

void makeGradYBlock (vector<Triplet<double> >& tripletList,
                     const int M0,
                     const int N0,
                     const int M,
                     const int N,
                     const double h) {

  for (int i = 0; i < M * (N - 1); ++ i) {
    tripletList.push_back (Triplet<double> (M0 + i, N0 + i, -1 / h));
    tripletList.push_back (Triplet<double> (M0 + i, N0 + M + i, 1 / h));
  }
}

void makeDivXBlock (vector<Triplet<double> >& tripletList,
                    const int M0,
                    const int N0,
                    const int M,
                    const int N,
                    const double h) {

  for (int i = 0; i < M; ++i) {
    for (int x = 0; x < N - 1; ++x) {
      tripletList.push_back (Triplet<double> (M0 + i * N + x, N0 + i * (N - 1) + x, 1 / h));
      tripletList.push_back (Triplet<double> (M0 + i * N + x + 1, N0 + i * (N - 1) + x, -1 / h));
    }
  }
}

void makeDivYBlock (vector<Triplet<double> >& tripletList,
                    const int M0,
                    const int N0,
                    const int M,
                    const int N,
                    const double h) {

  for (int i = 0; i < M * (N - 1); ++i) {
    tripletList.push_back (Triplet<double> (M0 + i, N0 + i, 1 / h));
    tripletList.push_back (Triplet<double> (M0 + i + N, N0 + i, -1 / h));
  }
}

SparseMatrix<double> makeForcingMatrix (const int M,
                                        const int N) {
  SparseMatrix<double> forcingMatrix (3 * M * N - M - N, 2 * M * N - M - N);
  forcingMatrix.reserve (2 * M * N - M - N);

  vector<Triplet<double> > tripletList;

  for (int i = 0; i < 2 * M * N - M - N; ++i)
    tripletList.push_back (Triplet<double> (i, i, 1));

  forcingMatrix.setFromTriplets (tripletList.begin(), tripletList.end());

  return forcingMatrix;
}

SparseMatrix<double> makeBoundaryMatrix (const int M,
                                         const int N,
                                         const double h,
                                         const double viscosity) {
  SparseMatrix<double> boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  boundaryMatrix.reserve (3 * M * N - M - N);

  vector<Triplet<double> > tripletList;
  
  makeBCLaplacianXBlock (tripletList, 0,                 0,     M, N, h, viscosity);
  makeBCLaplacianYBlock (tripletList, M * (N - 1),       2 * M, M, N, h, viscosity);
  makeBCDivXBlock       (tripletList, 2 * M * N - M - N, 0,     M, N, h);
  makeBCDivYBlock       (tripletList, 2 * M * N - M - N, 2 * M, M, N, h);

  boundaryMatrix.setFromTriplets (tripletList.begin(), tripletList.end());

  return boundaryMatrix;
}

void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity) {

  for (int i = 0; i < M; ++i) {
    tripletList.push_back (Triplet<double> (M0 + i * (N - 1),           N0 + i * 2,     viscosity / (h * h)));
    tripletList.push_back (Triplet<double> (M0 + (i + 1) * (N - 1) - 1,       N0 + i * 2 + 1, viscosity / (h * h)));
  }
}

void makeBCLaplacianYBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double viscosity) {
  for (int i = 0; i < N; ++i) {
    tripletList.push_back (Triplet<double> (M0 + i,               N0 + i, viscosity / (h * h)));
    tripletList.push_back (Triplet<double> (M0 + (M - 2) * N + i, N0 + N + i, viscosity / (h * h)));
  }
}

void makeBCDivXBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h) {
  for (int i = 0; i < N; ++i) {
    tripletList.push_back (Triplet<double> (M0 + i * M,           N0 + i * 2,      1 / h));
    tripletList.push_back (Triplet<double> (M0 + (i + 1) * M - 1, N0 + i * 2 + 1, -1 / h));
  }
}

void makeBCDivYBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h) {

  for (int i = 0; i < N; ++i) {
    tripletList.push_back (Triplet<double> (M0 + i, N0 + i, 1 / h));
    tripletList.push_back (Triplet<double> (M0 + (M - 1) * N + i, N0 + N + i, -1 / h));
  }
}
*/
