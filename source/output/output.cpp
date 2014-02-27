#include <iostream>

#include "boost/lexical_cast.hpp"

#include "H5Cpp.h"

#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "output/output.h"

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif
using namespace std;

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
    parser.queryParamString ("outputFormat", outputFormat, string("hdf5"));
    parser.queryParamString ("outputFilename", outputFilename, string("output"));

    parser.pop();
  }
}

OutputStructure::~OutputStructure () {}

void OutputStructure::writeDefaultFile() {
  if (outputFormat == "hdf5") {
    this->writeHDF5File();
  } else {
    cerr << "Unknown output format " << outputFormat << " specified in parameters." << 
            "Shutting down now." << endl << endl;
    exit (-1);
  }
}

void OutputStructure::writeDefaultFile (const int timestep) {
  if (outputFormat == "hdf5") {
    this->writeHDF5File (timestep);
  } else {
    cerr << "Unknown output format " << outputFormat << " specified in parameters." << 
            "Shutting down now." << endl << endl;
    exit (-1);
  }
}

void OutputStructure::writeHDF5File() {
  H5File outputFile (H5std_string (outputFilename + ".h5"), H5F_ACC_TRUNC);
  
  cout << " ==>> Outputting current data to file \"" << outputFilename << ".h5\"" << endl;
  
  DataSet dataset;
  hsize_t dimsf[2] = {0, 0};

  // Write x positions
  double *x = new double [N + 1];
  for (int i = 0; i < (N + 1); ++i)
    x [i] = dx * i;
  dimsf[0] = N + 1;
  dataset = outputFile.createDataSet ("X", 
                                      PredType::NATIVE_DOUBLE, 
                                      DataSpace (1, dimsf));
  dataset.write (x, 
                 PredType::NATIVE_DOUBLE);

  // Write y positions
  double *y = new double [M + 1];
  for (int j = 0; j < (M + 1); ++j)
    y [j] = dx * j;
  dimsf[0] = M + 1;
  dataset = outputFile.createDataSet ("Y",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (1, dimsf));
  dataset.write (y, 
                 PredType::NATIVE_DOUBLE);


  dimsf[0] = M; dimsf[1] = N;
  
  // Write temperature
  double * temperatureData = geometry.getTemperatureData();
  dataset = outputFile.createDataSet ("Temperature",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (temperatureData, 
                 PredType::NATIVE_DOUBLE);

  // Write pressure
  double * pressureData = geometry.getPressureData();
  dataset = outputFile.createDataSet ("Pressure",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (pressureData, 
                 PredType::NATIVE_DOUBLE);

  // Write U Velocity
  double * interpolatedUVelocityData = new double [M * N];
  double * uVelocityData = geometry.getUVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (j == 0) {
        interpolatedUVelocityData [i * N + j] = (uVelocityBoundaryData [i * 2]          +
                                                 uVelocityData         [i * (N - 1) + j]) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityData [i * N + j] = (uVelocityData         [i * (N - 1) + (j - 1)] +
                                                 uVelocityBoundaryData [i * 2       + 1      ]) / 2;
      } else {
        interpolatedUVelocityData [i * N + j] = (uVelocityData [i * (N - 1) + (j - 1)] +
                                                 uVelocityData [i * (N - 1) + j]      ) / 2;
      }
  dataset = outputFile.createDataSet ("UVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedUVelocityData,
                 PredType::NATIVE_DOUBLE);
  
  // Write V Velocity
  double * interpolatedVVelocityData = new double [M * N];
  double * vVelocityData             = geometry.getVVelocityData();
  double * vVelocityBoundaryData     = geometry.getVVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (i == 0) {
        interpolatedVVelocityData [i * N + j] = (vVelocityBoundaryData [j]               +
                                                 vVelocityData         [i       * N + j]) / 2;
      } else if (i == (M - 1)) {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityBoundaryData [N + j]          ) / 2;
      } else {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityData         [i       * N + j]) / 2;
      }
  dataset = outputFile.createDataSet ("VVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedVVelocityData,
                 PredType::NATIVE_DOUBLE);
                                                        
  // Write Viscosity
  double * viscosityData = geometry.getViscosityData();
  dataset = outputFile.createDataSet ("Viscosity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (viscosityData,
                 PredType::NATIVE_DOUBLE);
 
  delete[] x;
  delete[] y;
  delete[] interpolatedUVelocityData;
  delete[] interpolatedVVelocityData;
}

void OutputStructure::writeHDF5File (const int timestep) {
  H5File outputFile (H5std_string (outputFilename + 
                                   boost::lexical_cast<std::string> (timestep)  + 
                                   ".h5"),
                     H5F_ACC_TRUNC);

  cout << "==>> Outputting current data to \"" << outputFilename << timestep << ".h5\"" << endl;

  DataSet dataset;
  hsize_t dimsf[2] = {0, 0};

  // Write x positions
  double *x = new double [N + 1];
  for (int i = 0; i < (N + 1); ++i)
    x [i] = dx * i;
  dimsf[0] = N + 1;
  dataset = outputFile.createDataSet ("X", 
                                      PredType::NATIVE_DOUBLE, 
                                      DataSpace (1, dimsf));
  dataset.write (x, 
                 PredType::NATIVE_DOUBLE);

  // Write y positions
  double *y = new double [M + 1];
  for (int j = 0; j < (M + 1); ++j)
    y [j] = dx * j;
  dimsf[0] = M + 1;
  dataset = outputFile.createDataSet ("Y",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (1, dimsf));
  dataset.write (y, 
                 PredType::NATIVE_DOUBLE);


  dimsf[0] = M; dimsf[1] = N;
  
  // Write temperature
  double * temperatureData = geometry.getTemperatureData();
  dataset = outputFile.createDataSet ("Temperature",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (temperatureData, 
                 PredType::NATIVE_DOUBLE);

  // Write pressure
  double * pressureData = geometry.getPressureData();
  dataset = outputFile.createDataSet ("Pressure",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (pressureData, 
                 PredType::NATIVE_DOUBLE);

  // Write U Velocity
  double * interpolatedUVelocityData = new double [M * N];
  double * uVelocityData = geometry.getUVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (j == 0) {
        interpolatedUVelocityData [i * N + j] = (uVelocityBoundaryData [i * 2]          +
                                                 uVelocityData         [i * (N - 1) + j]) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityData [i * N + j] = (uVelocityData         [i * (N - 1) + (j - 1)] +
                                                 uVelocityBoundaryData [i * 2       + 1      ]) / 2;
      } else {
        interpolatedUVelocityData [i * N + j] = (uVelocityData [i * (N - 1) + (j - 1)] +
                                                 uVelocityData [i * (N - 1) + j]      ) / 2;
      }
  dataset = outputFile.createDataSet ("UVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedUVelocityData,
                 PredType::NATIVE_DOUBLE);
  
  // Write V Velocity
  double * interpolatedVVelocityData = new double [M * N];
  double * vVelocityData             = geometry.getVVelocityData();
  double * vVelocityBoundaryData     = geometry.getVVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (i == 0) {
        interpolatedVVelocityData [i * N + j] = (vVelocityBoundaryData [j]               +
                                                 vVelocityData         [i       * N + j]) / 2;
      } else if (i == (M - 1)) {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityBoundaryData [N + j]          ) / 2;
      } else {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityData         [i       * N + j]) / 2;
      }
  dataset = outputFile.createDataSet ("VVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedVVelocityData,
                 PredType::NATIVE_DOUBLE);
                                                        
  // Write Viscosity
  double * viscosityData = geometry.getViscosityData();
  dataset = outputFile.createDataSet ("Viscosity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (viscosityData,
                 PredType::NATIVE_DOUBLE);
 
  delete[] x;
  delete[] y;
  delete[] interpolatedUVelocityData;
  delete[] interpolatedVVelocityData;
}
