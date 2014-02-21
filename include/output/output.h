#pragma once

#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"

class OutputStructure {
  public:
    OutputStructure (ParamParser& pp, GeometryStructure& gs, ProblemStructure& ps);
    ~OutputStructure();

    void writeFile ();

  private:
    ParamParser&       parser;
    GeometryStructure& geometry;
    ProblemStructure&  problem;

    int M;
    int N;

    double dx;

    string outputFormat;
    string outputFilename;
};
