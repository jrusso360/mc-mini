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

ProblemStructure::ProblemStructure 
    (ParamParser&       pp,
     GeometryStructure& gs) :
    parser (pp),
    geometry (gs) {
  M = geometry.getM();
  N = geometry.getN();

  if (parser.push ("problemParams")) {
    parser.getParamDouble   ("cfl",              cfl);

    parser.getParamDouble   ("startTime",        time);
    parser.getParamDouble   ("endTime",          endTime);
    timestepNumber = 0;

    parser.queryParamDouble ("xExtent",          xExtent, 0.0);
    parser.queryParamDouble ("yExtent",          yExtent, 0.0);
    
    // Exactly one of xExtent or yExtent must be set. The other will be set
    // depending upon the first.
    assert (((xExtent == 0) ^ (yExtent == 0)));
    
    if (xExtent == 0) {
      dx      = yExtent / double(M);
      xExtent = dx * double(N);
    } else {
      dx      = xExtent / double(N);
      yExtent = dx * double(M);
    }

    parser.getParamDouble   ("diffusivity",      diffusivity);
    parser.getParamDouble   ("buoyancy",         buoyancy);

    parser.queryParamString ("forcingModel",     forcingModel,     "tauBenchmark");
    parser.queryParamString ("temperatureModel", temperatureModel, "constant");
    parser.queryParamString ("viscosityModel",   viscosityModel,   "constant");
    parser.queryParamString ("boundaryModel",    boundaryModel,    "tauBenchmark");

    parser.queryParamString ("outputFile",       outputFile,       "output.h5");

    parser.pop();
  }
}

ProblemStructure::~ProblemStructure() {
}                                   


void ProblemStructure::advanceTimestep() {
  time += dt;
  timestepNumber++;
}

double ProblemStructure::getH() {
  return dx;
}

double ProblemStructure::getTime() {
  return time;
}

double ProblemStructure::getEndTime() {
  return endTime;
}
