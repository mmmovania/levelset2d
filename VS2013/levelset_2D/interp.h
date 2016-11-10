/*
 *  interp.h
 */

#include "common.h"
namespace interp {
	float spline( float **d, float x, float y, int w, int h );
	float linear ( float **d, float x, float y, int w, int h );
	float interp ( float **d, float x, float y, int w, int h );
	///
	void setInterpMethod( int num );
}