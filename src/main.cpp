#include "controller.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <sys/time.h>
#elif defined(WIN32)
#include "glut.h"
#include <windows.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#include <sys/time.h>
#endif

int win_x = 512;
int win_y = 512;
int prev_x, prev_y;
int mstat = 1;

static void display(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	controller::display();
	glutSwapBuffers();
}

static void init( int gn ) {
	glClearColor(0.0, 0.0, 0.0, 1.0);
	controller::init(gn);
}

static void reshape(int w, int h) {
	win_x = w;
	win_y = h;
	controller::reshape(w,h);
}

static void idle ( void ) {
	glutPostRedisplay ();
}

static void keyboard( unsigned char key, int x, int y ) {
	if( key == '\e' ) exit(0);
	controller::keyDown(key);
	glutPostRedisplay ();
}

static void mouse ( int button, int state, int x, int y ) {
	prev_x = x;
	prev_y = y;
	mstat = state;
	controller::mouse( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y, ! state );
}

static void motion ( int x, int y ) {
	if( mstat == 0 ) {
		controller::motion( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y,
						 (x-prev_x)/(GLdouble)win_x, -(y-prev_y)/(GLdouble)win_y );
	}
	prev_x = x;
	prev_y = y;
}


int main (int argc, char * argv[]) {
	
	int grid_size = 64;
	if( argc == 2  ) {
		sscanf( argv[1], "%d", &grid_size );
	}
	
	glutInit(&argc, argv);
#if _OPENMP
	printf( "Number of threads: %d\n", omp_get_num_procs() );
#endif
	glutInitDisplayMode(GLUT_RGBA | GL_DOUBLE);
	glutInitWindowPosition ( 100, 100 );
	glutInitWindowSize ( win_x, win_y );
	glutCreateWindow(argv[0]);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc (motion);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	init(grid_size);
	glutMainLoop();
	return 0;
}
