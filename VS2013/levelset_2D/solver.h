/*
 *  solver.h
 *  smoke
 *
 */

#include "common.h"

namespace solver {
	// Solve Ax = b
	float solve( float **A, float **x, float **b, int n, char subcell, char solver );
}
