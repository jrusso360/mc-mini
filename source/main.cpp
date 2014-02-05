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

int main(int argc, char ** argv) {

  int M, N;
  double viscosity, y_extent, h;

  string paramFile = argv[1];

  paramParser pp(paramFile);

  pp.getParamInt ("M", M); // Number of rows
  pp.getParamInt ("N", N); // Number of columns

  pp.getParamDouble ("viscosity", viscosity);

  pp.getParamDouble ("y_extent", y_extent);
  h = y_extent / M;

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
  
  SparseMatrix<double> A (3 * M * N - M - N, 3 * M * N - M - N);
  SparseForms::makeA (A, M, N, h, viscosity);
  SparseMatrix<double> forcingMatrix (3 * M * N - M - N, 2 * M * N - M - N);
  SparseForms::makeForcingMatrix (forcingMatrix, M, N);
  SparseMatrix<double> boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  SparseForms::makeBoundaryMatrix (boundaryMatrix, M, N, h, viscosity);
  A.makeCompressed(); 
  
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
  
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  solver.analyzePattern (A);
  solver.compute (A);
  solnVect = solver.solve(rhsVect);
  
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
  return 0;
}
