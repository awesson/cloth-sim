#include "Particle.h"
//#include <GLUT/glut.h>

const GLfloat green[] = {0, 1.0, 0, 1};
GLfloat white_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat shininess[] = { 100.0 };

Particle::Particle(const Vec3 & i_ConstructPos, double i_inv_mass, int i_id) :
	ConstructPos(i_ConstructPos), Position(i_ConstructPos),
        Velocity(make_vector(0.0, 0.0, 0.0)), forces(make_vector(0.0, 0.0, 0.0)),
        deriv_position(make_vector(0.0, 0.0, 0.0)),
        deriv_velocity(make_vector(0.0, 0.0, 0.0)),
 		normal(make_vector(0.0, 0.0, 0.0)), id(i_id), inv_mass(i_inv_mass)
{
}

Particle::Particle(Particle const & b) :
	ConstructPos(b.ConstructPos), Position(b.Position),
        Velocity(b.Velocity), forces(b.forces),
        deriv_position(b.deriv_position),
        deriv_velocity(b.deriv_velocity),
 		normal(b.normal), id(b.id), inv_mass(b.inv_mass)
{
}

Particle::~Particle(void)
{
}

void Particle::reset()
{
	Position = ConstructPos;
	Velocity = make_vector(0.0, 0.0, 0.0);
    forces = make_vector(0.0, 0.0, 0.0);
}
void Particle::draw()
{
	const double h = 0.03;
	glColor3f(0.f, 1.f, 0.f);
	glBegin(GL_QUADS);
	glVertex3f(Position[0]-h/2.0, Position[1]-h/2.0, Position[2]);
	glVertex3f(Position[0]+h/2.0, Position[1]-h/2.0, Position[2]);
	glVertex3f(Position[0]+h/2.0, Position[1]+h/2.0, Position[2]);
	glVertex3f(Position[0]-h/2.0, Position[1]+h/2.0, Position[2]);
	glEnd();
}
