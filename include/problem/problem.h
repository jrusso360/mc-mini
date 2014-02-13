#pragma once

#include "paramParser/parser.h"

class ProblemStructure {
  public:
    ProblemStructure (ParamParser& pp, GeometryStructure& gs);
    ~ProblemStructure();

    // Initialize problem (set dt, temperature, boundary.)
    void initializeProblem();
    void initializeTimestep();
    void initializeTemperature();
    void initializeBoundary();

    // Main loop equations.
    void setForcingTerms();
    void solveStokes();
    void solveAdvectionDiffusion();
    void advanceTimestep();

    // Advection/diffusion methods
    void forwardEuler  (double dt);
    void crankNicolson (double dt);

    void outputPressure();
    void outputVelocity();
    void outputBoundaryVelocity();
    void outputForcing();
    void outputTemperature();

    double getTime();
    double getEndTime();

  private:
    ParamParser& parser;
    GeometryStructure& geometry;

    string forcingModel;
    string temperatureModel;
    string boundaryModel;

    int M;
    int N;

    double cfl;

    double time;
    double endTime;
    double dt;
    int timestepNumber;

    double xExtent;
    double yExtent;
    double dx;

    double viscosity;
    double diffusivity;
    double buoyancy;
};
