/*
 *  solver.cpp
 *  smoke
 *
 */
#include <minmax.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"
#include "utility.h"

static char subcell = 0;
static char solver_mode = 0;


// Clamped Fetch
static float x_ref( float **A, float **x, int fi, int fj, int i, int j, int n ) {
	i = min(max(0,i),n-1);
	j = min(max(0,j),n-1);
	if( A[i][j] < 0.0 ) return x[i][j];
	return subcell ? A[i][j]/fmin(0.00001f,A[fi][fj])*x[fi][fj] : 0.0;
}

// Ans = Ax
static void compute_Ax( float **A, float **x, float **ans, int n ) {
	float h2 = 1.0/(n*n);
	FOR_EVERY_CELL(n) {
		if( A[i][j] < 0.0 ) {
			ans[i][j] = (4.0*x[i][j]-x_ref(A,x,i,j,i+1,j,n)-x_ref(A,x,i,j,i-1,j,n)-x_ref(A,x,i,j,i,j+1,n)-x_ref(A,x,i,j,i,j-1,n))/h2;
		} else {
			ans[i][j] = 0.0;
		}
	} END_FOR
}

// ans = x^T * x
static float product( float **A, float **x, float **y, int n ) {
	float ans = 0.0;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) ans += x[i][j]*y[i][j];
		}
	}
	return ans;
}

// x = 0
static void clear( float **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = 0.0;
		}
	}
}

static void flip( float **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = -x[i][j];
		}
	}
}

// x <= y
static void copy( float **x, float **y, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = y[i][j];
		}
	}
}
				 
// Ans = x + a*y
static void op( float **A, float **x, float **y, float **ans, float a, int n ) {
	static float **tmp = alloc2D<float>(n);
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) tmp[i][j] = x[i][j]+a*y[i][j];
		}
	}
	copy(ans,tmp,n);
}

// r = b - Ax
static void residual( float **A, float **x, float **b, float **r, int n ) {
	compute_Ax(A,x,r,n);
	op( A, b, r, r, -1.0, n );
}

static inline float square( float a ) {
	return a*a;
}

static float A_ref( float **A, int i, int j, int qi, int qj, int n ) {
	if( i<0 || i>n-1 || j<0 || j>n-1 || A[i][j]>0.0 ) return 0.0;
	if( qi<0 || qi>n-1 || qj<0 || qj>n-1 || A[qi][qj]>0.0 ) return 0.0;
	return -1.0;
}

static float A_diag( float **A, int i, int j, int n ) {
	float diag = 4.0;
	if( A[i][j] > 0.0 ) return diag;
	int q[][2] = { {i-1,j}, {i+1,j}, {i,j-1}, {i,j+1} };
	for( int m=0; m<4; m++ ) {
		int qi = q[m][0];
		int qj = q[m][1];
		if( qi<0 || qi>n-1 || qj<0 || qj>n-1 ) diag -= 1.0;
		else if( A[qi][qj] > 0.0 && subcell ) {
			diag -= A[qi][qj]/fmin(0.00001f,A[i][j]);
		}
	}
	return diag;
}

static float P_ref( float **P, int i, int j, int n ) {
	if( i<0 || i>n-1 || j<0 || j>n-1 ) return 0.0;
	return P[i][j];
}

static void buildPreconditioner( float **P, float **A, int n ) {
	clear(P,n);
	float t = solver_mode == 2 ? 0.97 : 0.0;
	float a = 0.25;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) {
				float left = A_ref(A,i-1,j,i,j,n)*P_ref(P,i-1,j,n);
				float bottom = A_ref(A,i,j-1,i,j,n)*P_ref(P,i,j-1,n);
				float mleft = A_ref(A,i-1,j,i,j,n)*A_ref(A,i,j-1,i,j,n)*square(P_ref(P,i-1,j,n));
				float mbottom = A_ref(A,i,j-1,i,j,n)*A_ref(A,i-1,j,i,j,n)*square(P_ref(P,i,j-1,n));
				
				float diag = A_diag( A, i, j, n );
				float e = diag - square(left) - square(bottom) - t*( mleft + mbottom );
				if( e < a*diag ) e = diag;
				P[i][j] = 1.0/sqrtf(e);
			}
		}
	}
}

static void applyPreconditioner( float **z, float **r, float **P, float **A, int n ) {
	if( solver_mode == 0 ) {
		copy(z,r,n);
		return;
	}
	
	static float **q = alloc2D<float>(n);
	clear(q,n);
	
	// Lq = r
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) {
				float left = A_ref(A,i-1,j,i,j,n)*P_ref(P,i-1,j,n)*P_ref(q,i-1,j,n);
				float bottom = A_ref(A,i,j-1,i,j,n)*P_ref(P,i,j-1,n)*P_ref(q,i,j-1,n);
				
				float t = r[i][j] - left - bottom;
				q[i][j] = t*P[i][j];
			}
		}
	}
	
	// L^T z = q
	for( int i=n-1; i>=0; i-- ) {
		for( int j=n-1; j>=0; j-- ) {
			if( A[i][j] < 0.0 ) {
				float right = A_ref(A,i,j,i+1,j,n)*P_ref(P,i,j,n)*P_ref(z,i+1,j,n);
				float top = A_ref(A,i,j,i,j+1,n)*P_ref(P,i,j,n)*P_ref(z,i,j+1,n);
				
				float t = q[i][j] - right - top;
				z[i][j] = t*P[i][j];
			}
		}
	}
}

static void conjGrad( float **A, float **P, float **x, float **b, int n ) {
	// Pre-allocate Memory
	static float **r = alloc2D<float>(n);
	static float **z = alloc2D<float>(n);
	static float **s = alloc2D<float>(n);
	
	clear(x,n);									// p = 0
	copy(r,b,n);								// r = b
	applyPreconditioner(z,r,P,A,n);				// Apply Conditioner z = f(r)
	copy(s,z,n);								// s = z
	
	float a = product( A, z, r, n );			// a = z . r
	for( int k=0; k<n*n; k++ ) {
		compute_Ax( A, s, z, n );				// z = applyA(s)
		float alpha = a/product( A, z, s, n );	// alpha = a/(z . s)
		op( A, x, s, x, alpha, n );				// p = p + alpha*s
		op( A, r, z, r, -alpha, n );			// r = r - alpha*z;
		float error2 = product( A, r, r, n );	// error2 = r . r
		if( error2/(n*n) < 1.0e-6 ) break;
		applyPreconditioner(z,r,P,A,n);			// Apply Conditioner z = f(r)
		float a2 = product( A, z, r, n );		// a2 = z . r
		float beta = a2/a;
		op( A, z, s, s, beta, n );				// s = z + beta*s
		a = a2;
	}
}

float solver::solve( float **A, float **x, float **b, int n, char subcell_aware, char solver_type ) {
	static float **r = alloc2D<float>(n);
	static float **P = alloc2D<float>(n);
	clear(r,n);
	
	// Save Mode
	subcell = subcell_aware;
	solver_mode = solver_type;
	
	// Flip Divergence
	flip(b,n);
	
	// Build Modified Incomplete Cholesky Precondioner Matrix
	if( solver_mode >= 1 ) buildPreconditioner(P,A,n);
	
	// Conjugate Gradient Method
	conjGrad(A,P,x,b,n);

	residual(A,x,b,r,n);
	float res = sqrt(product( A, r, r, n ))/(n*n);
	// printf( "Residual = %e\n", res );
	return res;
}