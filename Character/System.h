#pragma once

#include "SpringForce.h"
#include <vector>
#include <stdlib.h>

#define G 0.0004
#define EPSILON 1.0e-30
#define del_x .01
#define Ks 0.5
#define Kd 0.8
#define restitution 0.2

class System
{
public:
	
	System(int N);
	~System(void);

	void deriv_eval(std::vector<Particle*> &);
	void get_state(std::vector<Particle*> &);
	void get_forces(std::vector<SpringForce*> &);
	bool resolve_collisions(std::vector<Particle*>& prev_state, Vec3 pos, Vec3 vel, double radius);

	void set_state(std::vector<Particle*>);
	// allow for adding spring forces after initialization
	void add_springForce(SpringForce*);

	void pop_springForce();
	int size();
	void draw();

private:

	std::vector<Particle*> pVector;
	std::vector<SpringForce*> forceVector;
	int N;

};
