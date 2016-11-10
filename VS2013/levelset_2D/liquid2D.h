/*
 *  liquid2D.h
 */

#include "common.h"

namespace liquid2D {
	void init( int n );
	void display();
	void keyDown( char key );
	void mouse( float x, float y, int state );
	void motion( float x, float y, float dx, float dy );
}