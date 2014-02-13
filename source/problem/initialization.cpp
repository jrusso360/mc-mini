#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"
#include "paramParser/parser.h"
#include "geometry/geometry.h"
#include "problem/problem.h"

/*
 *
 * Data initialization routines
 *
 */ 

void ProblemStructure::initializeProblem() {
  initializeTimestep();
  initializeTemperature();
  initializeBoundary();
}

void ProblemStructure::initializeTimestep() {
  dt = cfl * dx / diffusivity;
  int nTimestep = (endTime - time) / dt;
  if (abs (nTimestep * dt + time - endTime) > 1E-06) 
    dt = (endTime - time) / ++nTimestep;
}

void ProblemStructure::initializeTemperature() {
  double * temperatureData = geometry.getTemperatureData();

  double scale;

  if (parser.push ("problemParams")) {
    if (parser.tryPush ("initialTemperature")) {
      parser.queryParamDouble ("temperatureScale", scale, 100.0);

      parser.pop();
    } else {
      scale = 100.0;
    }
    parser.pop();
  }
  
  if (temperatureModel == "constant") {
    
    for (int i = 0; i < N; ++i) 
      for (int j = 0; j < M; ++j) 
        temperatureData[j * N + i] = scale;
   
  } else if (temperatureModel == "sineWave") {
    int xModes;
    int yModes;

    if (parser.push ("problemParams")) {
      if (parser.tryPush ("initialTemperature")) {
        parser.queryParamInt ("xModes", xModes, 2);
        parser.queryParamInt ("yModes", yModes, 2);

        parser.pop();
      } else {
        xModes = 2;
        yModes = 2;
      }

      parser.pop();
    }

    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        temperatureData[j * N + i] = sin ((i + 0.5) * dx * xModes * M_PI / xExtent) * 
                                     sin ((j + 0.5) * dx * yModes * M_PI / yExtent) * scale;
  } else if (temperatureModel == "squareWave") {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        temperatureData[j * N + i] = 
            ((xExtent * 0.25 < (i + 0.5) * dx && (i + 0.5) * dx < xExtent * 0.75) &&
             (yExtent * 0.25 < (j + 0.5) * dx && (j + 0.5) * dx < yExtent * 0.75)) ? scale : 0;
  }
}

void ProblemStructure::initializeBoundary() {
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();

  if (boundaryModel == "tauBenchmark") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryData [i * 2 + j] = cos (j * N * dx) * sin ((i + 0.5) * dx);
  
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryData [i * N + j] = -sin ((j + 0.5) * dx) * cos (i * M * dx);
  }
}
