/*
 *  profile.cpp
 */

#include "profile.h"
#include "levelset2D.h"
#include <cmath>
namespace {
	int gn;
	float maxdist = 0.15f;
}

static bool sphere( float x, float y ) {
	return hypot(x-0.5f,y-0.75f) < 0.15f;
}

void profile::init( int n ) {
	gn = n;
	levelset2D::init(n);
	
	// Build Simple LevelSet
	levelset2D::buildLevelset(sphere,maxdist);
}

static void flow( float x, float y, float &u, float &v, float &dt ) {
	u = -(y-0.5f);
	v = x-0.5f;
	dt = 0.01f;
}

void profile::display() {
	// Advect
	levelset2D::advect(flow);
	levelset2D::redistance(maxdist);	
	levelset2D::display(true);
}

void profile::keyDown( unsigned char key ) {
	if( key == 'r' ) {
		levelset2D::redistance(maxdist);
	}
}