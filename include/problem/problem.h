#pragma once

#include "parser/parser.h"

class ProblemStructure {
  public:
    ProblemStructure (ParamParser& pp, GeometryStructure& gs);
    ~ProblemStructure();

    // Initialize problem (set dt, temperature, boundary.)
    void initializeProblem();
    void initializeTimestep();
    void initializeViscosity();
    void initializeTemperature();
    void initializeBoundary();

    // Main loop equations.
    void updateForcingTerms();
    void updateViscosity();
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
    void outputViscosity();
    void outputTemperature();
    void outputBoundaryTemperature();

    void outputH5();

    double getH();
    double getTime();
    double getEndTime();

  private:
    ParamParser& parser;
    GeometryStructure& geometry;

    string forcingModel;
    string temperatureModel;
    string viscosityModel;
    string boundaryModel;
    string outputFile;

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
