#pragma once

#include <Vector/Vector.hpp>
#include <Graphics/Graphics.hpp>

class Particle
{
public:

	Particle(const Vec3 & ConstructPos, double mass, int i_id);
	Particle(Particle const & b);
	virtual ~Particle(void);

	Particle &operator=(Particle const & b) {
		ConstructPos = b.ConstructPos;
		Position = b.Position;
		Velocity = b.Velocity;
		forces = b.forces;
		deriv_position = b.deriv_position;
		deriv_velocity = b.deriv_velocity;
		normal = b.normal;
		id = b.id;
		inv_mass = b.inv_mass;
		return *this;
	}

	void reset();
	void draw();

	Vec3 ConstructPos;
	Vec3 Position;
	Vec3 Velocity;
	Vec3 forces;
	// the derivative of position and velocity respectively
	Vec3 deriv_position;
	Vec3 deriv_velocity;
	// normal of the particle used for rendering
	Vec3 normal;
	// particle's numbering used in calculating lambda
	int id;
	double inv_mass;
	
};