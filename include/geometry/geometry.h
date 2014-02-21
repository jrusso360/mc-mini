#pragma once

#include "parser/parser.h"

class GeometryStructure {
  public:
    GeometryStructure (ParamParser& pp);
    ~GeometryStructure();
  
    int getM();
    int getN();

    // Stokes Data
    double * getStokesData();
    // Subsections of Stokes Data
    double * getVelocityData();
    double * getUVelocityData();
    double * getVVelocityData();
    double * getPressureData();

    // Velocity Boundary Data
    double * getVelocityBoundaryData();
    // Subsections of Velocity Boundary Data
    double * getUVelocityBoundaryData();
    double * getVVelocityBoundaryData();

    // Forcing Data
    double * getForcingData();
    // Subsections of Forcing Data
    double * getUForcingData();
    double * getVForcingData();

    // Viscosity Data 
    double * getViscosityData();

    // Temperature Data
    double * getTemperatureData();

    // Temperature Boundary Data
    double * getTemperatureBoundaryData();
    // Subsections of Temperature Boundary Data
    double * getUTemperatureBoundaryData();
    double * getVTemperatureBoundaryData();
  
  private:
    // Geometry sizes
    int M; // Rows in domain (y size)
    int N; // Columns in domain (x size)

    // Data pointers
    // Stokes solution data contains velocity and pressure solutions.
    double * stokesData;
    // Velocity boundary data contains prescribed velocity boundary values.
    double * velocityBoundaryData;

    // Forcing data contains stokes equation forcing terms.
    double * forcingData;

    // Viscosity data
    double * viscosityData;

    // Temperature data
    double * temperatureData;
    // Temperature Boundary Data
    double * temperatureBoundaryData;
};
