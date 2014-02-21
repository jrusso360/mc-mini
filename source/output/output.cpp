#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "output/output.h"

OutputStructure::OutputStructure (ParamParser&       pp,
                                  GeometryStructure& gs,
                                  ProblemStructure&  ps) :
    parser   (pp),
    geometry (gs),
    problem  (ps) {
  M  = geometry.getM(); 
  N  = geometry.getN();
  dx = problem.getH();

  if (parser.push ("outputParams")) {
    parser.queryParamString ("outputFilename", outputFilename, string("output.h5"));

    parser.pop();
  }

}

OutputStructure::~OutputStructure () {}
