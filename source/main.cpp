#include <iostream>
#include <iomanip>
#include <cassert>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "matrixForms/denseForms.h"
#include "paramParser/parser.h"
#include "geometry/geometry.h"
#include "problem/problem.h"


using namespace Eigen;
using namespace std;

int main(int argc, char ** argv) {

  assert (argc > 1);

  ParamParser pp(string{argv[1]});
  GeometryStructure geometry (pp);
  ProblemStructure problem (pp, geometry);

  problem.initializeProblem();

  problem.updateForcingTerms();
  problem.solveStokes();
  problem.solveAdvectionDiffusion();
  problem.advanceTimestep();
/*
  problem.outputTemperature();
  problem.outputForcing();
  problem.outputViscosity();
  problem.outputBoundaryVelocity();
  problem.outputPressure();
  problem.outputVelocity();

  problem.outputH5();
*/
  int M = geometry.getM();
  int N = geometry.getN();
  double h = problem.getH();

  MatrixXd analyticU;
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uSolnMatrix (geometry.getUVelocityData(), M, (N - 1));
  analyticU.resizeLike(uSolnMatrix);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N - 1; ++j)
      analyticU (i, j) = cos ((j + 1) * h) * sin ((i + 0.5) * h);

  MatrixXd analyticV;
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vSolnMatrix (geometry.getVVelocityData(), (M - 1), N);
  analyticV.resizeLike(vSolnMatrix);
  for (int i = 0; i < M - 1; ++i)
    for (int j = 0; j < N; ++j)
      analyticV (i, j) = - sin ((j + 0.5) * h) * cos ((i + 1) * h);

  cout << sqrt((uSolnMatrix - analyticU).squaredNorm() * h * h) 
       << "\t" 
       << sqrt((vSolnMatrix - analyticV).squaredNorm() * h * h) << endl;
  
  return 0;
}
