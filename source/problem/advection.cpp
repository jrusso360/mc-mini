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

// upwind method. Stable but inefficient.
void ProblemStructure::upwindMethod() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  double * uVelocityData = geometry.getUVelocityData();
  double * vVelocityData = geometry.getVVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();

  double leftVelocity, rightVelocity;
  double leftFlux, rightFlux;

  VectorXd temporaryVector = temperatureVector;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      leftVelocity = (j == 0)        ? 
                        uVelocityBoundaryData[2 * i]     : 
                        uVelocityData[i * (N - 1) + j];
      rightVelocity = (j == (N - 1)) ? 
                        uVelocityBoundaryData[2 * i + 1] : 
                        uVelocityData[i * (N - 1) + (j + 1)];

      leftFlux = rightFlux = 0;

      if (leftVelocity < 0) {
        if (j > 0) {
          leftFlux = temporaryVector(i * N + j) * leftVelocity * deltaT / (2 * h);
        }
      } else {
        if (j > 0) {
          leftFlux = temporaryVector(i * N + (j - 1)) * leftVelocity * deltaT / (2 * h);
        }
      }

      if (rightVelocity > 0) {
        if (j < (N - 1)) {
          rightFlux = -temporaryVector(i * N + j) * rightVelocity * deltaT / (2 * h);
        }
      } else {
        if (j < (N - 1)) {
          rightFlux = -temporaryVector(i * N + (j + 1)) * rightVelocity * deltaT / (2 * h);
        }
      }

      temperatureVector(i * N + j) = temporaryVector(i * N + j) + leftFlux + rightFlux;
    }
  }

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      leftVelocity = (i == 0) ?
                       vVelocityBoundaryData[j] :
                       vVelocityData[(i - 1) * N + j];
      rightVelocity = (i == (M - 1)) ?
                        vVelocityBoundaryData[N + j] :
                        vVelocityData[i * N + j];
      
      leftFlux = rightFlux = 0;

      if (leftVelocity < 0) {
        if (i > 0) {
          leftFlux = temperatureVector(i * N + j) * leftVelocity * deltaT / h;
        }
      } else {
        if (i > 0) {
          leftFlux = temperatureVector((i - 1) * N + j) * leftVelocity * deltaT / h;
        }
      }

      if (rightVelocity > 0) {
        if (i < (M - 1)) {
          rightFlux = -temperatureVector(i * N + j) * rightVelocity * deltaT / h;
        }
      } else {
        if (i < (M - 1)) {
          rightFlux = -temperatureVector((i + 1) * N + j) * rightVelocity * deltaT / h;
        }
      }

      temporaryVector (i * N + j) = temperatureVector (i * N + j) + leftFlux + rightFlux;
    }
  }

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      leftVelocity = (j == 0)        ? 
                        uVelocityBoundaryData[2 * i]     : 
                        uVelocityData[i * (N - 1) + j];
      rightVelocity = (j == (N - 1)) ? 
                        uVelocityBoundaryData[2 * i + 1] : 
                        uVelocityData[i * (N - 1) + (j + 1)];

      leftFlux = rightFlux = 0;

      if (leftVelocity < 0) {
        if (j > 0) {
          leftFlux = temporaryVector(i * N + j) * leftVelocity * deltaT / (2 * h);
        }
      } else {
        if (j > 0) {
          leftFlux = temporaryVector(i * N + (j - 1)) * leftVelocity * deltaT / (2 * h);
        }
      }

      if (rightVelocity > 0) {
        if (j < (N - 1)) {
          rightFlux = -temporaryVector(i * N + j) * rightVelocity * deltaT / (2 * h);
        }
      } else {
        if (j < (N - 1)) {
          rightFlux = -temporaryVector(i * N + (j + 1)) * rightVelocity * deltaT / (2 * h);
        }
      }

      temperatureVector(i * N + j) = temporaryVector(i * N + j) + leftFlux + rightFlux;
    }
  }
}
