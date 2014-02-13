#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"
#include "paramParser/parser.h"
#include "geometry/geometry.h"
#include "problem/problem.h"

using namespace Eigen;
using namespace std;

/*
 *
 * Main loop routines
 *
 */ 

// Update the forcing terms
// T -> F
void ProblemStructure::setForcingTerms() {
  double * uForcingData = geometry.getUForcingData();
  double * vForcingData = geometry.getVForcingData();
  
  if (forcingModel == "tauBenchmark") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        uForcingData [i * (N - 1) + j] = 3 * cos ((j + 1) * dx) * sin ((i + 0.5) * dx);

    for (int i = 0; i < (M - 1); ++i)
      for (int j = 0; j < N; ++j)
        vForcingData [i * N + j] = -sin ((j + 0.5) * dx) * cos ((i + 1) * dx);
  } 
}

// Solve the stokes equation
// F -> U X P
void ProblemStructure::solveStokes() {
  Map<VectorXd> stokesSolnVector (geometry.getStokesData(), M * (N - 1) + (M - 1) * N + M * N);
  
  double * viscosityData = geometry.getViscosityData();

  static double prevViscosity;

  static SparseMatrix<double> stokesMatrix   (3 * M * N - M - N, 3 * M * N - M - N);
  static SparseMatrix<double> forcingMatrix  (3 * M * N - M - N, 2 * M * N - M - N);
  static SparseMatrix<double> boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  if (prevViscosity != viscosity) {

    SparseForms::makeStokesMatrix   (stokesMatrix,   M, N, dx, viscosity);
    stokesMatrix.makeCompressed();
    SparseForms::makeForcingMatrix  (forcingMatrix,  M, N);
    forcingMatrix.makeCompressed();
    SparseForms::makeBoundaryMatrix (boundaryMatrix, M, N, dx, viscosity);
    boundaryMatrix.makeCompressed();

    solver.analyzePattern (stokesMatrix);
    solver.compute (stokesMatrix);
  
    prevViscosity = viscosity;
  }

  stokesSolnVector = solver.solve (forcingMatrix * Map<VectorXd>(geometry.getForcingData(), 2 * M * N - M - N) + boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * M + 2 * N));
}

// Solve the advection/diffusion equation
// U X T -> T
void ProblemStructure::solveAdvectionDiffusion() {
  crankNicolson (dt);
}
