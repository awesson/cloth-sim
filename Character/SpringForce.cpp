#include "SpringForce.h"
//#include <GLUT/glut.h>

SpringForce::SpringForce(Particle* i_p1, Particle* i_p2, double i_dist, double i_ks, double i_kd) :
  p1(i_p1), p2(i_p2), dist(i_dist), ks(i_ks), kd(i_kd) {}

void SpringForce::draw()
{
	glBegin( GL_LINES );
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f( p1->Position[0], p1->Position[1], p1->Position[2]);
	glColor3f(0.0, 0.5, 1.0);
	glVertex3f( p2->Position[0], p2->Position[1], p2->Position[2]);
	glEnd();
}

void SpringForce::add_force(){
	Vec3 dx = p1->Position - p2->Position;
	double norm_dx = length(dx);
	double v_dx;
	Vec3 force1;
	// so that if the particles are on top of each other the program does not blow up
	if(norm_dx < EPS){
		force1 = make_vector(INF, 0.0, 0.0);
	} else{
		double inv_norm_dx = 1.0/norm_dx;
		v_dx = (p1->Velocity - p2->Velocity)*dx*inv_norm_dx;
		force1 = -(ks*(norm_dx - dist) + kd*v_dx)*dx*inv_norm_dx;
	}
	p1->forces += force1;
	p2->forces -= force1;
}
