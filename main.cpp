#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

// THIS WORKS CORRECTLY! DO NOT TOUCH WITHOUT A BACKUP!

using namespace Eigen;
using namespace std;

MatrixXd makeDivX (const int, const int);
MatrixXd makeDivY (const int, const int);
MatrixXd makeGradX (const int, const int);
MatrixXd makeGradY (const int, const int);
MatrixXd makeLaplacianU (const int, const int);
MatrixXd makeLaplacianV (const int, const int);

MatrixXd makeBCLaplacianU (const int, const int);
MatrixXd makeBCLaplacianV (const int, const int);
MatrixXd makeBCDivX (const int, const int);
MatrixXd makeBCDivY (const int, const int);

int main () {
  const int M = 20; // Number of rows
  const int N = 40; // Number of columns

  const double viscosity = 1.0;

  const double h = M_PI / M;

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

  MatrixXd A              = MatrixXd::Zero (3 * M * N - M - N, 3 * M * N - M - N);
  MatrixXd forcingMatrix  = MatrixXd::Zero (3 * M * N - M - N, 2 * M * N - M - N);
  MatrixXd boundaryMatrix = MatrixXd::Zero (3 * M * N - M - N, 2 * M + 2 * N);

  A.block (0,                 0,                 M * (N - 1), M * (N - 1)) = -viscosity * makeLaplacianU (M, N) / (h * h);
  A.block (M * (N - 1),       M * (N - 1),       (M - 1) * N, (M - 1) * N) = -viscosity * makeLaplacianV (M, N) / (h * h);
  A.block (0,                 2 * M * N - M - N, M * (N - 1), M * N)       = makeGradX (M, N) / h;
  A.block (M * (N - 1),       2 * M * N - M - N, (M - 1) * N, M * N)       = makeGradY (M, N) / h;
  A.block (2 * M * N - M - N, 0,                 M * N,       M * (N - 1)) = makeDivX (M, N) / h;
  A.block (2 * M * N - M - N, M * (N - 1),       M * N,       (M - 1) * N) = makeDivY (M, N) / h;
  
  forcingMatrix.block (0, 0, 2 * M * N - M - N, 2 * M * N - M - N)  = MatrixXd::Identity (2 * M * N - M - N, 2 * M * N - M - N);

  boundaryMatrix.block (0,                 0,     M * (N - 1), 2 * M) = viscosity * makeBCLaplacianU (M, N) / (h * h);
  boundaryMatrix.block (M * (N - 1),       2 * M, (M - 1) * N, 2 * N) = viscosity * makeBCLaplacianV (M, N) / (h * h);
  boundaryMatrix.block (2 * M * N - M - N, 0,     M * N,       2 * M) = makeBCDivX (M, N) / h;
  boundaryMatrix.block (2 * M * N - M - N, 2 * M, M * N,       2 * N) = makeBCDivY (M, N) / h;
  
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

  /*
  cout << uSolnMatrix.colwise().reverse() << endl << endl
       << vSolnMatrix.colwise().reverse() << endl << endl
       << pSolnMatrix.colwise().reverse() << endl << endl;
*/
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

  cout << (uSolnMatrix - analyticU).norm() << "\t" << (vSolnMatrix - analyticV).norm() << endl;

  return 0;
}

// Create laplacian operator matrix for U
MatrixXd makeLaplacianU (const int M, const int N) {
  // Laplacian maps U -> U_0, so requires an (M * (N - 1))X(M * (N - 1)) matrix.
  MatrixXd laplacian = MatrixXd::Zero (M * (N - 1), M * (N - 1));
  MatrixXd laplacianBlock = MatrixXd::Zero (N - 1, N - 1);

  laplacianBlock.diagonal () = VectorXd::Constant (N - 1,  -4);
  laplacianBlock.diagonal (-1) = laplacianBlock.diagonal (1) = VectorXd::Constant (N - 2, 1);

  for (int i = 0; i < M; ++i) {
    laplacianBlock.diagonal () = VectorXd::Constant (N - 1, -4);
    if ((i == 0) || (i == M - 1))
      laplacianBlock.diagonal () = VectorXd::Constant (N - 1, -5);

    laplacian.block (i * (N - 1), i * (N - 1), (N - 1), (N - 1)) = laplacianBlock;
  }

  laplacian.diagonal (N - 1)    = VectorXd::Constant ((N - 1) * (M - 1), 1);
  laplacian.diagonal (-(N - 1)) = VectorXd::Constant ((N - 1) * (M - 1), 1);

  return laplacian;
}

// Create laplacian operator matrix for V.
MatrixXd makeLaplacianV (const int M, const int N) {
  // Laplacian maps V -> V_0, so requires an ((M - 1) * N)X((M - 1) * N) matrix.
  MatrixXd laplacian = MatrixXd::Zero ((M - 1) * N, (M - 1) * N);
  MatrixXd laplacianBlock = MatrixXd::Zero (N, N);
  
  laplacianBlock.diagonal () = VectorXd::Constant (N, -4);
  laplacianBlock.diagonal (-1) = laplacianBlock.diagonal (1) = VectorXd::Constant (N - 1, 1);
  laplacianBlock (0, 0) = laplacianBlock (N - 1, N - 1) = -5;

  for (int i = 0; i < M - 1; ++i) 
    laplacian.block (i * N, i * N, N, N) = laplacianBlock;

  laplacian.diagonal (N) = laplacian.diagonal (-N) = VectorXd::Constant ((M - 2) * N, 1);

  return laplacian;
}

// Create X-dimension gradient operator matrix
MatrixXd makeGradX (const int M, const int N) {
  
  return -1 * makeDivX (M, N).transpose();
}

// Create Y-dimension gradient operator matrix
MatrixXd makeGradY (const int M, const int N) {
  
  return -1 * makeDivY (M, N).transpose();
}

// Create X-dimension divergence operator matrix
MatrixXd makeDivX (const int M, const int N) {
  MatrixXd divX     = MatrixXd::Zero (M * N, M * (N - 1));
  MatrixXd divBlock = MatrixXd::Zero (N,     N - 1);

  divBlock.diagonal ()   = VectorXd::Constant (N - 1,  1);
  divBlock.diagonal (-1) = VectorXd::Constant (N - 1, -1);

  for (int i = 0; i < M; ++i) 
    divX.block (i * N, i * (N - 1), N, N - 1) = divBlock;

  return divX;
}

// Create Y-dimension divergence operator matrix
MatrixXd makeDivY (const int M, const int N) {
  MatrixXd divY = MatrixXd::Zero (M * N, (M - 1) * N);

  divY.diagonal ()   = VectorXd::Constant ((M - 1) * N,  1);
  divY.diagonal (-N) = VectorXd::Constant ((M - 1) * N, -1);

  return divY;
}



MatrixXd makeBCLaplacianU (const int M, const int N) {
  MatrixXd laplacianBC      = MatrixXd::Zero (M * (N - 1), 2 * M);
  MatrixXd laplacianBCBlock = MatrixXd::Zero (N - 1,       2);

  laplacianBCBlock (0, 0) = laplacianBCBlock (N - 2, 1) = 1;

  for (int i = 0; i < M; ++i)
    laplacianBC.block (i * (N - 1), i * 2, N - 1, 2) = laplacianBCBlock;

  return laplacianBC;
}

MatrixXd makeBCLaplacianV (const int M, const int N) {
  MatrixXd laplacianBC      = MatrixXd::Zero ((M - 1) * N, 2 * N);
  MatrixXd laplacianBCBlock = MatrixXd::Zero (N, N);

  laplacianBCBlock.diagonal () = VectorXd::Constant (N, 1);

  laplacianBC.block (0, 0, N, N) = laplacianBCBlock;
  laplacianBC.block ((M - 2) * N, N, N, N) = laplacianBCBlock;

  return laplacianBC;
}

MatrixXd makeBCDivX (const int M, const int N) {
  MatrixXd divBC      = MatrixXd::Zero (M * N, 2 * M);
  MatrixXd divBCBlock = MatrixXd::Zero (N,     2);

  divBCBlock (0, 0) = 1; divBCBlock (N - 1, 1) = -1;

  for (int i = 0; i < M; ++i)
    divBC.block (i * N, i * 2, N, 2) = divBCBlock;

  return divBC;
}

MatrixXd makeBCDivY (const int M, const int N) {
  MatrixXd divBC      = MatrixXd::Zero     (M * N, 2 * N);

  divBC.block (0,           0, N, N) =  MatrixXd::Identity (N, N);
  divBC.block ((M - 1) * N, N, N, N) =  -1 * MatrixXd::Identity (N, N);

  return divBC;
}
