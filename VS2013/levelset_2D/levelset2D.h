/*
 *  levelset2D.h
 */

#include "common.h"

namespace levelset2D {
	void init( int n );
	void buildLevelset( bool (*func)(float x, float y), float tolerance );
	void advect( void (*func)( float x, float y, float &u, float &v, float &dt ) );
	void redistance( float tolerance );
	void extrapolate( float **q, char **region );
	void display( bool cell_centered );
	void keyDown( char key );
	///
	float getLevelSet( int i, int j );
	void getLevelSet( float **dists );
	void setLevelSet( int i, int j, float d );
	float getVolume();
	///
	void setVisibility( bool show_grid, bool show_dist, bool show_region );
}
