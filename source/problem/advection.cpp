#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "debug.h"

using namespace Eigen;
using namespace std;

// Upwind method. Stable but inefficient.
void ProblemStructure::upwindMethod() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), N - 1, M);
  DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), N, M - 1);
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;
  double G_Flux[], F_Flux[];

  VectorXd temporaryVector = temperatureVector;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j <  N; ++j) {
      // Find all four edge velocities for the current cell. 
      leftVelocity   = (j == 0) ?
                        uVelocityBoundaryWindow (0, i) :
                        uVelocityWindow (j - 1, i);
      rightVelocity  = (j == (N - 1)) ?
                        uVelocityBoundaryWindow (1, i) :
                        uVelocityWindow (j, i);
      bottomVelocity = (i == 0) ?
                        vVelocityBoundaryWindow (j, 0) :
                        vVelocityWindow (j, i - 1);
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryWindow (j, 1) :
                        vVelocityWindow (j, i);

		
      
      // Solve the Riemann problem on the neighboring velocities and calculate
      // the fluxes accross each edge.
		
        if (leftVelocity < 0) {
			F_Flux[j-1] = temporaryVector (i * N + j) * leftVelocity;
        } else {
			F_Flux[j-1] = temporaryVector (i * N + (j - 1)) * leftVelocity;
        }
		

	  
        if (rightVelocity > 0) {
			F_Flux[j] = temporaryVector (i * N + j) * rightVelocity;
        } else {
			F_Flux[j]= temporaryVector (i * N + (j + 1)) * rightVelocity;
        }
	  

		
        if (bottomVelocity < 0) {
          G_Flux[i-1] = temporaryVector (i * N + j) * bottomVelocity;
        } else {
          G_Flux[i-1] = temporaryVector ((i - 1) * N + j) * bottomVelocity;
        }
		  

		
        if (topVelocity > 0) {
          G_Flux[i] = temporaryVector (i * N + j) * topVelocity;
        } else {
          G_Flux[i] = temporaryVector ((i + 1) * N + j) * topVelocity;
        }
		

      // Standard Conservative Finite Difference Update

      temperatureVector (i * N + j) = temporaryVector (i * N + j) + 
										(deltaT / h) * (F_Flux(j-1,i) - F_Flux(j,i)) +
                                        (deltaT / h) * (G_Flux(j,i-1) - G_Flux(j,i));
    }
  }

	//Conservation check by summing all values in the temperature vector
	double total_temperature = std::accumulate(temperatureVector.begin(),temperatureVector.end(),0);
	cout << "\The total temperature is: " << total_temperature << endl;

}

void ProblemStructure::laxWendroff() {
}

void ProblemStructure::frommMethod() {
  // Temperature data (MxN cell-centered grid)
  DataWindow<double> temperatureWindow (geometry.getTemperatureData(), N, M);
  // Temperature boundary data (2xN transverse boundary grid)
  DataWindow<double> temperatureBoundaryWindow (geometry.getTemperatureBoundaryData(), N, 2);
  
  // U Velocity Data (Mx(N-1) lateral offset grid)
  DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), N - 1, M);
  // V Velocity Data ((M-1)xN transverse offset grid)
  DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), N, M - 1);
  
  // U Velocity Boundary Data (Mx2 lateral boundary grid)
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);
  // V Velocity Boundary Data (2xN transverse boundary grid)
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);

  // Initilize half-time data
  // Half-time temperature data (MxN cell-centered grid)
  static DataWindow<double> halfTimeTemperatureWindow (new double[M * N], N, M);
  // Half-time U-Offset temperature data (Mx(N-1) lateral offset grid)
  static DataWindow<double> halfTimeUOffsetTemperatureWindow (new double[M * (N - 1)], N - 1, M);
  // Half-time V-offset temperature data ((M-1)xN transverse offset grid)
  static DataWindow<double> halfTimeVOffsetTemperatureWindow (new double[(M - 1) * N], N, M - 1);

  // Half-time Forcing Data (for use in the Stokes solve)
  static double * halfTimeForcingData = new double[2 * M * N - M - N];
  // Half-time U forcing data (Mx(N-1) lateral offset grid)
  static double * halfTimeUForcingData = halfTimeForcingData;
  static DataWindow<double> halfTimeUForcingWindow (halfTimeUForcingData, N - 1, M); 
  // Half-time V forcing data ((M-1)xN transverse offset grid)
  static double * halfTimeVForcingData = halfTimeForcingData + M * (N - 1);
  static DataWindow<double> halfTimeVForcingWindow (halfTimeVForcingData, N, M - 1);
  
  // Half-time Stokes solution data (for use in the Stokes solve)
  static double * halfTimeStokesSolnData = new double[3 * M * N - M - N];
  // Half-time U velocity data (Mx(N-1) lateral offset grid)
  static double * halfTimeUVelocityData = halfTimeStokesSolnData;
  static DataWindow<double> halfTimeUVelocityWindow (halfTimeUVelocityData, N - 1, M);
  // Half-time V velocity data ((M-1)xN transverse offset grid)
  static double * halfTimeVVelocityData = halfTimeStokesSolnData + M * (N - 1);
  static DataWindow<double> halfTimeVVelocityWindow (halfTimeVVelocityData, N, M - 1);
  // Half-time pressure data (MxN cell-centered grid) (not technically necessary)
  static double * halfTimePressureData = halfTimeStokesSolnData + 2 * M * N - M - N;
  static DataWindow<double> halfTimePressureWindow (halfTimePressureData, N, M);

  // Cell-centered velocities
  static double * cellCenteredVelocityData = new double[2 * M * N];
  // Cell-centered U velocity (MxN cell-centered grid)
  static double * cellCenteredUVelocityData = cellCenteredVelocityData;
  static DataWindow<double> cellCenteredUVelocityWindow (cellCenteredUVelocityData, N, M);
  // Cell-centered V velocity (MxN cell-centered grid)
  static double * cellCenteredVVelocityData = cellCenteredVelocityData + M * N;
  static DataWindow<double> cellCenteredVVelocityWindow (cellCenteredVVelocityData, N, M);

  Map<VectorXd> halfTimeStokesSolnVector (halfTimeStokesSolnData, 3 * M * N - M - N);
  VectorXd temporaryTemperature = Map<VectorXd> (geometry.getTemperatureData(), N * M);
  
  double leftNeighborT, rightNeighborT, bottomNeighborT, topNeighborT;
  
  // Calculate cell-centered velocities (MxN cell-centered grid)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double leftVelocity   = (j == 0) ?
                              uVelocityBoundaryWindow (0, i) :
                              uVelocityWindow (j - 1, i);
      double rightVelocity  = (j == (N - 1)) ?
                              uVelocityBoundaryWindow (1, i) :
                              uVelocityWindow (j, i);
      double bottomVelocity = (i == 0) ?
                              vVelocityBoundaryWindow (j, 0) :
                              vVelocityWindow (j, i - 1);
      double topVelocity    = (i == (M - 1)) ?
                              vVelocityBoundaryWindow (j, 1) :
                              vVelocityWindow (j, i);
      
      cellCenteredVVelocityWindow (j, i) = (topVelocity + bottomVelocity) / 2;
      cellCenteredUVelocityWindow (j, i) = (leftVelocity + rightVelocity) / 2;
      // Check if the cell-centered velocity is machine-epsilon, and make it true zero if so.
      if (abs (cellCenteredUVelocityWindow (j, i)) < 1E-10) 
        cellCenteredUVelocityWindow (j, i) = 0;
      if (abs (cellCenteredVVelocityWindow (j, i)) < 1E-10)
        cellCenteredVVelocityWindow (j, i) = 0;
    }
  }

  #ifdef DEBUG
    cout << "<Half-time Cell-Centered Velocities Calculated>" << std::endl;
    cout << "<Cell-Centered U Velocity>" << endl;
    cout << cellCenteredUVelocityWindow.displayMatrix() << endl;
    cout << "<Cell-Centered V Velocity>" << endl;
    cout << cellCenteredVVelocityWindow.displayMatrix() << endl;
  #endif
  
  // Calculate temperatures at half-time level
  // Calculate half-time U-offset temperatures (Mx(N-1) lateral offset grid)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < (N - 1); ++j) {
      halfTimeUOffsetTemperatureWindow (j, i) = 0;
      // Flag to show which Riemann Problem solution to use. We use a bitflag
      // for the riemmanFlag since the case where velocities do not share the 
      // same sign is merely the average of the other two cases.
      int riemannFlag = 0b00;
      if (cellCenteredUVelocityWindow (j,     i) >= 0 && 
          cellCenteredUVelocityWindow (j + 1, i) >= 0) {
        riemannFlag = 0b01;
      } else if (cellCenteredUVelocityWindow (j,     i) <= 0 &&
                 cellCenteredUVelocityWindow (j + 1, i) <= 0) {
        riemannFlag = 0b10;
      } else {
        riemannFlag = 0b11;
      }

      // One or both velocities are positive. Take the flux from the left
      // neighboring cell.
      if (riemannFlag & 0b01) {
        if (j == 0) {
          leftNeighborT = temperatureWindow (j, i);
        } else {
          leftNeighborT = temperatureWindow (j - 1, i);
        }
        rightNeighborT = temperatureWindow (j + 1, i);

        halfTimeUOffsetTemperatureWindow (j, i) += 
            temperatureWindow (j, i) +
              (h / 2 - deltaT / 2 * cellCenteredUVelocityWindow (j, i)) * 
                (rightNeighborT - leftNeighborT) / (2 * h);
      }
      
      // One or both velocities are negative. Take the flux from the right
      // neighboring cell.
      if (riemannFlag & 0b10) {
        leftNeighborT = temperatureWindow (j, i);
        if (j == (N - 2)) {
          rightNeighborT = temperatureWindow (j + 1, i);
        } else {
          rightNeighborT = temperatureWindow (j + 2, i);
        }

        halfTimeUOffsetTemperatureWindow (j, i) += 
            temperatureWindow (j + 1, i) -
              (h / 2 + deltaT / 2 * cellCenteredUVelocityWindow (j + 1, i)) * 
                (rightNeighborT - leftNeighborT) / (2 * h); 
      
      }

      // Velocities are in opposing directions. Take the average of the fluxes
      // from the left and right neighboring cells
      if (riemannFlag == 0b11) halfTimeUOffsetTemperatureWindow (j, i) /= 2;

    }
  }

  #ifdef DEBUG
    cout << "<Half-Time Offset Temperatures>" << endl;
    cout << "<Half-Time U-Offset Temperature>" << endl;
    cout << halfTimeUOffsetTemperatureWindow.displayMatrix() << endl;
  #endif

  // Calculate half-time V-offset temperatures ((M-1)xN transverse offset grid)
  for (int i = 0; i < (M - 1); ++i) {
    for (int j = 0; j < N; ++j) {
      halfTimeVOffsetTemperatureWindow (j, i) = 0;
      // Flag to show which Riemann Problem solution to use. We use a bitflag
      // for the riemannFlag since the case where velocities do not share the
      // same sign is merely the average of the other two cases.
      int riemannFlag = 0b00;

      if (cellCenteredVVelocityWindow (j, i    ) >= 0 && 
          cellCenteredVVelocityWindow (j, i + 1) >= 0) {
        // One or both velocities are positive
        riemannFlag = 0b01;
      } else if (cellCenteredVVelocityWindow (j, i    ) <= 0 && 
                 cellCenteredVVelocityWindow (j, i + 1) <= 0) {
        // One or both velocities are negative
        riemannFlag = 0b10;
      } else {
        // One velocity is positive, one is negative.
        riemannFlag = 0b11;
      }

      // One or both velocities are positive. Take the flux from the bottom
      // neighboring cell
      if (riemannFlag & 0b01) {
        bottomNeighborT = temperatureWindow (j, i + 1);
        if (i == 0) {
          topNeighborT = temperatureWindow (j, i);
        } else {
          topNeighborT = temperatureWindow (j, i - 1);
        }

        halfTimeVOffsetTemperatureWindow (j, i) +=
            temperatureWindow (j, i) -
              (h / 2 + deltaT / 2 * cellCenteredVVelocityWindow (j, i)) * 
               (topNeighborT - bottomNeighborT) / (2 * h);
      }
      
      // One or both velocities are negative. Take the flux from the top
      // neighboring cell
      if (riemannFlag & 0b10) {
        if (i == (M - 2)) {
          bottomNeighborT = temperatureWindow (j, i + 1);
        } else {
          bottomNeighborT = temperatureWindow (j, i + 2);
        }
        topNeighborT = temperatureWindow (j, i);

        halfTimeVOffsetTemperatureWindow (j, i) += 
            temperatureWindow (j, i + 1) +
              (h / 2 - deltaT / 2 * cellCenteredVVelocityWindow (j, i + 1)) * 
               (topNeighborT - bottomNeighborT) / (2 * h);
      }

      // Velocities are in opposing directions. Take the average of the fluxes
      // from the left and right neighboring cells.
      if (riemannFlag == 0b11) {
        halfTimeVOffsetTemperatureWindow (j, i) /= 2;
      }
    }
  }

  #ifdef DEBUG
    cout << "<Half-Time V-Offset Temperature>" << endl;
    cout << halfTimeVOffsetTemperatureWindow.displayMatrix() << endl << endl;
  #endif

  double leftHalfTimeNeighborT, rightHalfTimeNeighborT, 
         bottomHalfTimeNeighborT, topHalfTimeNeighborT; 

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double diffusionWeight = 4;
      if (j == 0) {
        leftHalfTimeNeighborT = temperatureWindow (j, i);
        leftNeighborT         = temperatureWindow (j, i);
        diffusionWeight = 3;
      } else {
        leftHalfTimeNeighborT = halfTimeUOffsetTemperatureWindow (j - 1, i); 
        leftNeighborT         = temperatureWindow (j - 1, i);
      }
      if (j == (N - 1)) {
        rightHalfTimeNeighborT = temperatureWindow (j, i); 
        rightNeighborT         = temperatureWindow (j, i);
        diffusionWeight = 3;
      } else {
        rightHalfTimeNeighborT = halfTimeUOffsetTemperatureWindow (j, i);
        rightNeighborT         = temperatureWindow (j + 1, i);
      }

      if (i == 0) {
        bottomHalfTimeNeighborT = temperatureBoundaryWindow (j, 0);
        bottomNeighborT         = temperatureBoundaryWindow (j, 0);
      } else {
        bottomHalfTimeNeighborT = halfTimeVOffsetTemperatureWindow (j, i - 1);
        bottomNeighborT         = temperatureWindow (j, i - 1);
      }

      if (i == (M - 1)) {
        topHalfTimeNeighborT = temperatureBoundaryWindow (j, 1);
        topNeighborT         = temperatureBoundaryWindow (j, 1);
      } else {
        topHalfTimeNeighborT = halfTimeVOffsetTemperatureWindow (j, i);
        topNeighborT         = temperatureWindow (j, i + 1);
      }

      halfTimeTemperatureWindow (j, i) = 
        (rightHalfTimeNeighborT + leftHalfTimeNeighborT + 
         topHalfTimeNeighborT + bottomHalfTimeNeighborT) / 4 - 
         (diffusionWeight * temperatureWindow (j, i) + 
           leftNeighborT + rightNeighborT + bottomNeighborT + topNeighborT) 
         * deltaT * diffusivity / (h * h);
    }
  }
  
  #ifdef DEBUG
    cout << "<Half-Time Temperature Data>" << endl;
    cout << halfTimeTemperatureWindow.displayMatrix() << endl << endl;
  #endif
  
  // Calculate half-time forcing
  if (forcingModel == "tauBenchmark") {
    // Benchmark taken from Tau (1991; JCP Vol. 99)
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        halfTimeUForcingWindow (j, i) = 3 * cos ((j + 1) * h) * sin ((i + 0.5) * h);

    for (int i = 0; i < M - 1; ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (j, i) = -sin ((j + 0.5) * h) * cos ((i + 1) * h);

  } else if (forcingModel == "solCXBenchmark" ||
             forcingModel == "solKZBenchmark") {
    // solCX Benchmark taken from Kronbichler et al. (2011)
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        halfTimeUForcingWindow (j, i) = 0;

    for (int i = 0; i < M - 1; ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (j, i) = - sin((i + 0.5) * M_PI * h) * cos ((j + 1) * M_PI * h);

  } else if (forcingModel == "vorticalFlow") {
    for (int i = 0; i < M; ++i) 
      for (int j = 0; j < (N - 1); ++j)
        halfTimeUForcingWindow (j, i) = 3 * cos ((j + 1) * 2 * h) * sin ((i + 0.5) * 2 * h);

    for (int i = 0; i < (M - 1); ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (j, i) = -sin ((j + 0.5) * h) * cos ((i + 1) * h);

  } else if (forcingModel == "buoyancy") {
    double referenceTemperature;
    double densityConstant;
    double thermalExpansion;

    if (parser.push ("problemParams")) {
      if (parser.push ("buoyancyModelParams")) {
        parser.queryParamDouble ("referenceTemperature", referenceTemperature, 273.15);
        parser.queryParamDouble ("densityConstant",      densityConstant,      100.0);
        parser.queryParamDouble ("thermalExpansion",     thermalExpansion,       1.0);
        parser.pop();
      }
      parser.pop();
    }

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < (N - 1); ++j)
        halfTimeUForcingWindow (j, i) = 0;

    for (int i = 0; i < (M - 1); ++i)
      for (int j = 0; j < N; ++j) {
        halfTimeVForcingWindow (j, i) =  -1 * densityConstant *
                                          (1 - thermalExpansion * 
                                           ((halfTimeTemperatureWindow (j, i) + 
                                             halfTimeTemperatureWindow (j, i + 1)) / 2 -
                                            referenceTemperature));
      }
  } else {
    cout << "<Unexpected forcing model: \"" << forcingModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Half-Time Forcing Data>" << endl;
    cout << "<Half-Time U Forcing Data>" << endl;
    cout << halfTimeUForcingWindow.displayMatrix() << endl;
    cout << "<Half-Time V Forcing Data>" << endl;
    cout << halfTimeVForcingWindow.displayMatrix() << endl << endl;
  #endif

  // Initialize half-time solver
  static SparseMatrix<double> stokesMatrix   (3 * M * N - M - N, 3 * M * N - M - N);
  static SparseMatrix<double> forcingMatrix  (3 * M * N - M - N, 2 * M * N - M - N);
  static SparseMatrix<double> boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  static bool initialized;

  if (!(initialized) || !(viscosityModel=="constant")) {
    
    double * viscosityData = geometry.getViscosityData();

    SparseForms::makeStokesMatrix (stokesMatrix, M, N, h, viscosityData);
    stokesMatrix.makeCompressed();
    SparseForms::makeForcingMatrix (forcingMatrix, M, N);
    forcingMatrix.makeCompressed();
    SparseForms::makeBoundaryMatrix (boundaryMatrix, M, N, h, viscosityData);
    boundaryMatrix.makeCompressed();

    solver.analyzePattern (stokesMatrix);
    solver.factorize      (stokesMatrix);

    initialized = true;

    #ifdef DEBUG
      cout << "<Viscosity Data>" << endl;
      cout << DataWindow<double> (geometry.getViscosityData(), N + 1, M + 1).displayMatrix() << endl << endl;
    #endif
  }

  #ifdef DEBUG
    cout << "<Velocity Boundary Data>" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    cout << uVelocityBoundaryWindow.displayMatrix() << endl;
    cout << "<V Velocity Boundary Data>" << endl;
    cout << vVelocityBoundaryWindow.displayMatrix() << endl << endl;
  #endif

  // Solve stokes at the half-time to find velocities
  halfTimeStokesSolnVector = solver.solve 
           (forcingMatrix  * Map<VectorXd>(halfTimeForcingData, 2 * M * N - M - N) + 
            boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * M + 2 * N));

  for (int i = 0; i < M; i++)
    for (int j = 0; j < (N - 1); j++)
      if (abs(halfTimeUVelocityWindow (j, i)) < 10E-10)
        halfTimeUVelocityWindow (j, i) = 0;

  for (int i = 0; i < (M - 1); i++) 
    for (int j = 0; j < N; ++j)
      if (abs (halfTimeVVelocityWindow (j, i)) < 10E-10)
        halfTimeVVelocityWindow (j, i) = 0;

  #ifdef DEBUG
    cout << "<Half-Time Velocities>" << endl;
    cout << "<Half-Time U Velocity>" << endl;
    cout << halfTimeUVelocityWindow.displayMatrix() << endl;
    cout << "<Half-Time V Velocity>" << std::endl;
    cout << halfTimeVVelocityWindow.displayMatrix() << endl << endl;
  #endif

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;

  // Solve for full-time temperature
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        leftNeighborT = 0; 
        leftVelocity = 0;
      } else {
        leftNeighborT = halfTimeUOffsetTemperatureWindow (j - 1, i);
        leftVelocity = halfTimeUVelocityWindow (j - 1, i);
      }
      if (j == (N - 1)) {
        rightNeighborT = 0;
        rightVelocity = 0;
      } else {
        rightNeighborT = halfTimeUOffsetTemperatureWindow (j, i);
        rightVelocity = halfTimeUVelocityWindow (j, i);
      }

      if (i == 0) {
        bottomNeighborT = 0; 
        bottomVelocity = 0;
      } else {
        bottomNeighborT = halfTimeVOffsetTemperatureWindow (j, i - 1);
        bottomVelocity = halfTimeVVelocityWindow (j, i - 1);
      }

      if (i == (M - 1)) {
        topNeighborT = 0;
        topVelocity = 0;
      } else {
        topNeighborT = halfTimeVOffsetTemperatureWindow (j, i);
        topVelocity = halfTimeVVelocityWindow (j, i);
      }
 
      // Initialize the flux limiter to point to the desired limiter function.
      static double (ProblemStructure::*limiter) (double,double,double) = NULL;
      if (limiter == NULL) {
        if (fluxLimiter == "minmod") {
          limiter = &ProblemStructure::minmod;
        } else if (fluxLimiter == "superbee") {
          limiter = &ProblemStructure::superbee;
        } else if (fluxLimiter == "vanLeer") {
          limiter = &ProblemStructure::vanLeer;
        } else if (fluxLimiter == "none") {
          limiter = &ProblemStructure::minmod;
        }
      }

      double leftFirstOrderT = 0, rightFirstOrderT = 0, 
             bottomFirstOrderT = 0, topFirstOrderT = 0,
             secondLeftFirstOrderT = 0, secondBottomFirstOrderT = 0;
      double leftFirstOrderVelocity, rightFirstOrderVelocity,
             bottomFirstOrderVelocity, topFirstOrderVelocity;
      double lateralFlux = 0, transverseFlux = 0;

      if (fluxLimiter != "none") {
        leftFirstOrderVelocity   = (j == 0) ?
                                   uVelocityBoundaryWindow (0, i) :
                                   uVelocityWindow (j - 1, i);
        rightFirstOrderVelocity  = (j == (N - 1)) ?
                                   uVelocityBoundaryWindow (1, i) :
                                   uVelocityWindow (j, i);
        bottomFirstOrderVelocity = (i == 0) ?
                                   vVelocityBoundaryWindow (j, 0) :
                                   vVelocityWindow (j, i - 1);
        topFirstOrderVelocity    = (i == (M - 1)) ?
                                   vVelocityBoundaryWindow (j, 1) :
                                   vVelocityWindow (j, i);

        double leftFlux = 0, rightFlux = 0;
        double topFlux = 0, bottomFlux = 0;

        if (j > 0) {
          if (leftFirstOrderVelocity < 0) {
            leftFlux = temporaryTemperature (i * N + j) * leftFirstOrderVelocity * deltaT / h;
            leftFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            leftFlux = temporaryTemperature (i * N + (j - 1)) * leftFirstOrderVelocity * deltaT / h;
            leftFirstOrderT = temporaryTemperature (i * N + (j - 1));
          }
        }

        if (j > 1) {
          secondLeftFirstOrderT = temporaryTemperature (i * N + (j - 2));
        } else {
          // Handle the weird symmetric boundary condition using abs magic
          secondLeftFirstOrderT = temporaryTemperature (i * N + abs (1 - j));
        }

        if (j < (N - 1)) {
          if (rightFirstOrderVelocity > 0) {
            rightFlux = temporaryTemperature (i * N + j) * rightFirstOrderVelocity * deltaT / h;
            rightFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            rightFlux = temporaryTemperature (i * N + (j + 1)) * rightFirstOrderVelocity * deltaT / h;
            rightFirstOrderT = temporaryTemperature (i * N + (j + 1));
          }
        }

        if (i > 0) {
          if (bottomFirstOrderVelocity < 0) {
            bottomFlux = temporaryTemperature (i * N + j) * bottomFirstOrderVelocity * deltaT / h;
            bottomFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            bottomFlux = temporaryTemperature ((i - 1) * N + j) * bottomFirstOrderVelocity * deltaT / h;
            bottomFirstOrderT = temporaryTemperature ((i - 1) * N + j);
          }
        }

        if (i > 1) {
          secondBottomFirstOrderT = temporaryTemperature ((i - 2) * N + j);
        } else {
          secondBottomFirstOrderT = temperatureBoundaryWindow (0, j);
        }

        if (i < (M - 1)) {
          if (topFirstOrderVelocity > 0) {
            topFlux = temporaryTemperature (i * N + j) * topFirstOrderVelocity * deltaT / h;
            topFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            topFlux = temporaryTemperature ((i + 1) * N + j) * topFirstOrderVelocity * deltaT / h;
            topFirstOrderT = temporaryTemperature ((i + 1) * N + j);
          }
        }

        double leftPhi = (this->*limiter) (secondLeftFirstOrderT,
                                           leftFirstOrderT,
                                           temporaryTemperature (i * N + j));
        double rightPhi = (this->*limiter) (leftFirstOrderT, 
                                            temporaryTemperature (i * N + j), 
                                            rightFirstOrderT);
        double bottomPhi = (this->*limiter) (secondBottomFirstOrderT,
                                             bottomFirstOrderT,
                                             temporaryTemperature (i * N + j));
        double topPhi = (this->*limiter) (bottomFirstOrderT, 
                                          temporaryTemperature (i * N + j), 
                                          topFirstOrderT);
 
        lateralFlux = 
            ((1 - leftPhi) * leftFlux + leftPhi * (leftVelocity * leftNeighborT * deltaT / h)) - 
            ((1 - rightPhi) * rightFlux + rightPhi * rightVelocity * rightNeighborT * deltaT / h);
        transverseFlux = ((1 - bottomPhi) * bottomFlux + bottomPhi * bottomVelocity * bottomNeighborT * deltaT / h) -
                         ((1 - topPhi) * topFlux + topPhi * topVelocity * topNeighborT * deltaT / h);

        temperatureWindow (j, i) = temporaryTemperature (i * N + j) + transverseFlux + lateralFlux;
      } else 
        temperatureWindow (j, i) = temporaryTemperature (i * N + j) + deltaT / h * (leftVelocity * leftNeighborT - rightVelocity * rightNeighborT) + deltaT / h * (bottomVelocity * bottomNeighborT - topVelocity * topNeighborT);

      if (std::isnan((double)temperatureWindow (j, i))) {
        #ifdef DEBUG
          cout << "<NaN at " << i << "," << j << ">" << endl;
          cout << "<Neighbor values were: " << endl;
          cout << "\tleftNeighborT   = " << leftNeighborT << endl;
          cout << "\trightNeighborT  = " << rightNeighborT << endl;
          cout << "\tbottomNeighborT = " << bottomNeighborT << endl;
          cout << "\ttopNeighborT    = " << topNeighborT << ">" << endl;
          cout << "<Velocity values were: " << endl;
          cout << "\tleftVelocity    = " << leftVelocity << endl;
          cout << "\trightVelocity   = " << rightVelocity << endl;
          cout << "\tbottomVelocity  = " << bottomVelocity << endl;
          cout << "\ttopVelocity     = " << topVelocity << ">" << endl;
          cout << "<Shutting down now!>" << endl;
        #endif
        throw "found NaN";
      }
    }
  }

  #ifdef DEBUG
    cout << "<Full-Time Temperature Data>" << endl;
    cout << temperatureWindow.displayMatrix() << endl << endl;
  #endif
}
