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

void ProblemStructure::slopes() {
	
	
	// Compute the U-direction cell centered slopes for Temperature.
	
	for (int i = 0; i < M; ++i) {
	    for (int j = 0; j <  N; ++j) {
			// If both velocities move to the right (positive), use right neighbor temperature minus left neighbor temperature.
			if (cellCenteredUVelocityWindow (i, j) >= 0 &&
				cellCenteredUVelocityWindow (i + 1, j) >= 0){
				
				uCenteredTemperatureSlopes (i,j) = (temperature (i+1,j) - temperature (i-1,j)) / (2.0 * h)
			}
		
			// If both velocities move to the left (negative), use left neighbor temperature minus right neighbor temperature.
			else if (cellCenteredUVelocityWindow (i, j) <= 0 &&
					 cellCenteredUVelocityWindow (i + 1, j) <= 0){
				
				uCenteredTemperatureSlopes (i , j) = (temperatureWindow (i-1,j)) - temperatureWindow (i+1,j) / (2.0 * h)

			}
			
			// If velocities are in opposite directions, set U-direction temperature slope to zero
			else {
				uCenteredTemperatureSlopes (i , j) = 0.0
			}
		}
	}
	
	
	
	// Compute the V-direction cell centered slopes for Temperature.
	
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j <  N; ++j) {
			// If both velocities move to the right (positive), use right neighbor temperature minus left neighbor temperature.
			if (cellCenteredVVelocityWindow (i, j) >= 0 &&
				cellCenteredVVelocityWindow (i + 1, j) >= 0){
				
				vCenteredTemperatureSlopes (i,j) = (temperature (i,j+1) - temperature (i,j-1)) / (2.0 * h)
			}
			
			// If both velocities move to the left (negative), use left neighbor temperature minus right neighbor temperature.
			else if (cellCenteredvVelocityWindow (i, j) <= 0 &&
					 cellCenteredvVelocityWindow (i + 1, j) <= 0){
				
				vCenteredTemperatureSlopes (i , j) = (temperatureWindow (i,j-1) - temperatureWindow (i,j+1)) / (2.0 * h)
				
			}
			
			// If velocities are in opposite directions, set V-direction temperature slope to the average of the two.
			else {
				vCenteredTemperatureSlopes (i , j) = (temperatureWindow (i,j-1) + temperatureWindow (i,j+1)) / (2.0 * h)
			}
		}
	// return uCenteredTemperatureSlopes, vCenteredTemperatureSlopes
	}
	

				
				
				
				
				
			
