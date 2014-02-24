#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/geometry.h"
#include "problem/problem.h"

using namespace Eigen;
using namespace std;

void ProblemStructure::forwardEuler (double delta_t) {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> uTemperatureBoundaryVector (geometry.getUTemperatureBoundaryData(), M * 2);
  Map<VectorXd> vTemperatureBoundaryVector (geometry.getVTemperatureBoundaryData(), 2 * N);
  
  static SparseMatrix<double> B;

  static int oldSize;

  if (oldSize != M * N) {
    oldSize = M * N;
  }
}

void ProblemStructure::crankNicolson (double delta_t) {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);


}
