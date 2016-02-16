#include "System.h"
#include <algorithm>

System::System(int i_N)
{
	N = i_N;
	const double dist = 0.25;
    const Vec3 center = make_vector(0.0, 10.5, 20.0);
    const Vec3 x_offset = make_vector(dist, 0.0, 0.0);
    const Vec3 y_offset = make_vector(0.0, -dist, 0.0);

    // Create an array of NxN particles connected by warp, weft, shear, and bend springs.

    // add particles
    for(int i = 0; i < N; ++i){
        for(int k = 0; k < N; ++k){
            pVector.push_back(new Particle(center + k*x_offset + i*y_offset, 1.f, i*N + k));
        }
    }

	// make it so the top corners are stationary
	pVector[0]->inv_mass = 0.0;
	pVector[N-1]->inv_mass = 0.0;
 
    for(int i = 0; i < N; ++i){
        for(int k = 0; k < N; ++k){
            if(k != N-1){
                // weft springs
                forceVector.push_back(new SpringForce(pVector[N*i + k], pVector[N*i + k + 1], dist, Ks, Kd));

                // shear springs
				if(i != 0)
                	forceVector.push_back(new SpringForce(pVector[N*i + k], pVector[N*i + k - N + 1], sqrt(2)*dist, Ks, Kd));
				if(i != N - 1)
	                forceVector.push_back(new SpringForce(pVector[N*i + k], pVector[N*i + k + N + 1], sqrt(2)*dist, Ks, Kd));
            }

			// warp springs
			if(i != N - 1)
            	forceVector.push_back(new SpringForce(pVector[N*i + k], pVector[N*i + k + N], dist, Ks, Kd));

			// bend springs
			if(k < N - 2)
				forceVector.push_back(new SpringForce(pVector[N*i + k], pVector[N*i + k + 2], 2*dist, Ks, Kd));
			if(i < N - 2)
				forceVector.push_back(new SpringForce(pVector[N*i + k], pVector[N*i + k + 2*N], 2*dist, Ks, Kd));
        }
    }

	// disturb the cloth initially so it is not "stuck" in 2d
	// for(int k = N/2 -1; k < N/2 + 2; ++k){
	// 	pVector[k]->Velocity = make_vector(0.0,0.0,0.0001);
	// }
}

System::~System(void)
{
	pVector.clear();
    forceVector.clear();
}

void System::deriv_eval(std::vector<Particle*>& o_pVector){
	int size = pVector.size();
	int num_f = forceVector.size();

	// reset forces to just gravity
	for(int i = 0; i < size; ++i){
		pVector[i]->forces = make_vector(0.0, -G, 0.0);
	}

	// add spring forces
	for(int i = 0; i < num_f; ++i){
		forceVector[i]->add_force();
	}

	// set the derivative of position to the velocity and
	// the derivative of the velocity to the total force divided by the mass
	for(int i = 0; i < size; ++i){
		pVector[i]->deriv_position = pVector[i]->Velocity;
		pVector[i]->deriv_velocity = pVector[i]->forces*pVector[i]->inv_mass;
	}

	// return particle vector with the derivatives
	o_pVector = pVector;
}

/**
 * Tests and resolves collisions between the cloth and a spherical obstacle. 
 * prev_state will be updated to have a velocity which will no longer cause
 * a collision on the next time step assuming no external forces.
 **/
bool System::resolve_collisions(std::vector<Particle*>& prev_state, Vec3 obstacle_pos, Vec3 obstacle_vel, double radius)
{
	Vec3 rel_pos, rel_vel, normal;
	double rel_dist;
	bool is_collision = false;
	for(unsigned int i = 0; i < pVector.size(); ++i){
		rel_pos = pVector[i]->Position - obstacle_pos;
		rel_dist = length(rel_pos);
		rel_vel = pVector[i]->Velocity - obstacle_vel;
		normal = rel_pos / rel_dist;
		if(rel_dist <= radius && rel_vel*normal < 0){ // collision has occurred
			// reflect the particle off the obstacle in the direction of the normal
			prev_state[i]->Velocity -= (1+restitution)*(rel_vel*normal)*normal;
			is_collision = true;
		}
	}
	return is_collision;
}

void System::get_state(std::vector<Particle*>& o_pVector)
{
	o_pVector.resize(pVector.size());
	for(unsigned int i = 0; i < pVector.size(); ++i)
		*o_pVector[i] = *pVector[i];
}

void System::set_state(std::vector<Particle*> i_pVector)
{
	pVector.resize(i_pVector.size());
	for(unsigned int i = 0; i < pVector.size(); ++i)
		*pVector[i] = *i_pVector[i];
}

void System::get_forces(std::vector<SpringForce*>& o_forceVector)
{
        o_forceVector = forceVector;
}

void System::add_springForce(SpringForce* f)
{
    forceVector.push_back(f);
}

void System::pop_springForce()
{
    forceVector.pop_back();
}

int System::size()
{
        return pVector.size();
}

void System::draw()
{	
	// calculate the normals
	for(int i = 0; i < N; ++i){
		for(int k = 0; k < N; ++k){
			pVector[N*i + k]->normal = make_vector(0.0, 0.0, 0.0);
			if(i != N-1 && k != N-1)
				pVector[N*i + k]->normal += normalize(cross_product(pVector[N*i + k]->Position - pVector[N*i + k + 1]->Position,
												  					pVector[N*i + k]->Position - pVector[N*i + k + N]->Position));
			if(i != N-1 && k != 0)
				pVector[N*i + k]->normal += normalize(cross_product(pVector[N*i + k]->Position - pVector[N*i + k + N]->Position,
												  					pVector[N*i + k]->Position - pVector[N*i + k - 1]->Position));
			if(i != 0 && k != 0)
				pVector[N*i + k]->normal += normalize(cross_product(pVector[N*i + k]->Position - pVector[N*i + k - 1]->Position,
												  					pVector[N*i + k]->Position - pVector[N*i + k - N]->Position));
			if(i != 0 && k != N-1)
				pVector[N*i + k]->normal += normalize(cross_product(pVector[N*i + k]->Position - pVector[N*i + k - N]->Position,
												  					pVector[N*i + k]->Position - pVector[N*i + k + 1]->Position));
			pVector[N*i + k]->normal = normalize(pVector[N*i + k]->normal);
		}
	}	
	
	glEnable(GL_NORMALIZE);
	glColor3d(0.0, 0.6, 0.0);
	glBegin(GL_TRIANGLES);
	// render a front and back side of the cloth separated by del_x to get normals correct
	for(int i = 0; i < N; ++i){
		for(int k = 0; k < N; ++k){
			// top left triangles
			if(i != N-1 && k != N-1){
				// front side
				glNormal3d(pVector[N*i + k]->normal[0],
						   pVector[N*i + k]->normal[1], 
						   pVector[N*i + k]->normal[2]);
				glVertex3d(pVector[N*i + k]->Position[0],
				 		   pVector[N*i + k]->Position[1], 
						   pVector[N*i + k]->Position[2]);
				glNormal3d(pVector[N*i + k + 1]->normal[0],
						   pVector[N*i + k + 1]->normal[1], 
					       pVector[N*i + k + 1]->normal[2]);
				glVertex3d(pVector[N*i + k + 1]->Position[0],
				 		   pVector[N*i + k + 1]->Position[1],
				 		   pVector[N*i + k + 1]->Position[2]);
				glNormal3d(pVector[N*i + k + N]->normal[0],
						   pVector[N*i + k + N]->normal[1], 
						   pVector[N*i + k + N]->normal[2]);
				glVertex3d(pVector[N*i + k + N]->Position[0],
				 		   pVector[N*i + k + N]->Position[1],
				 		   pVector[N*i + k + N]->Position[2]);
				
				// back side
				glNormal3d(-pVector[N*i + k]->normal[0],
						   -pVector[N*i + k]->normal[1], 
						   -pVector[N*i + k]->normal[2]);
				glVertex3d(pVector[N*i + k]->Position[0] - del_x*pVector[N*i + k]->normal[0],
				 		   pVector[N*i + k]->Position[1] - del_x*pVector[N*i + k]->normal[1],
				 		   pVector[N*i + k]->Position[2] - del_x*pVector[N*i + k]->normal[2]);
				glNormal3d(-pVector[N*i + k + 1]->normal[0],
						   -pVector[N*i + k + 1]->normal[1], 
						   -pVector[N*i + k + 1]->normal[2]);
				glVertex3d(pVector[N*i + k + 1]->Position[0] - del_x*pVector[N*i + k + 1]->normal[0],
				 		   pVector[N*i + k + 1]->Position[1] - del_x*pVector[N*i + k + 1]->normal[1],
				 		   pVector[N*i + k + 1]->Position[2] - del_x*pVector[N*i + k + 1]->normal[2]);
				glNormal3d(-pVector[N*i + k + N]->normal[0],
						   -pVector[N*i + k + N]->normal[1], 
						   -pVector[N*i + k + N]->normal[2]);
				glVertex3d(pVector[N*i + k + N]->Position[0] - del_x*pVector[N*i + k + N]->normal[0],
				 		   pVector[N*i + k + N]->Position[1] - del_x*pVector[N*i + k + N]->normal[1],
				 		   pVector[N*i + k + N]->Position[2] - del_x*pVector[N*i + k + N]->normal[2]);
			}
			
			// bottom right triangles
			if(i != N-1 && k != 0){
				// front side
				glNormal3d(pVector[N*i + k]->normal[0],
						   pVector[N*i + k]->normal[1], 
						   pVector[N*i + k]->normal[2]);
				glVertex3d(pVector[N*i + k]->Position[0],
				 		   pVector[N*i + k]->Position[1],
				 		   pVector[N*i + k]->Position[2]);
				glNormal3d(pVector[N*i + k + N]->normal[0],
						   pVector[N*i + k + N]->normal[1], 
						   pVector[N*i + k + N]->normal[2]);
				glVertex3d(pVector[N*i + k + N]->Position[0],
				 		   pVector[N*i + k + N]->Position[1],
				 		   pVector[N*i + k + N]->Position[2]);
				glNormal3d(pVector[N*i + k + N-1]->normal[0],
						   pVector[N*i + k + N-1]->normal[1], 
						   pVector[N*i + k + N-1]->normal[2]);
				glVertex3d(pVector[N*i + k + N-1]->Position[0],
				 		   pVector[N*i + k + N-1]->Position[1],
				 		   pVector[N*i + k + N-1]->Position[2]);
				
				// back side
				glNormal3d(-pVector[N*i + k]->normal[0],
						   -pVector[N*i + k]->normal[1], 
						   -pVector[N*i + k]->normal[2]);
				glVertex3d(pVector[N*i + k]->Position[0] - del_x*pVector[N*i + k]->normal[0],
				 		   pVector[N*i + k]->Position[1] - del_x*pVector[N*i + k]->normal[1],
				 		   pVector[N*i + k]->Position[2] - del_x*pVector[N*i + k]->normal[2]);
				glNormal3d(-pVector[N*i + k + N]->normal[0],
						   -pVector[N*i + k + N]->normal[1], 
						   -pVector[N*i + k + N]->normal[2]);
				glVertex3d(pVector[N*i + k + N]->Position[0] - del_x*pVector[N*i + k + N]->normal[0],
				 		   pVector[N*i + k + N]->Position[1] - del_x*pVector[N*i + k + N]->normal[1],
				 		   pVector[N*i + k + N]->Position[2] - del_x*pVector[N*i + k + N]->normal[2]);
				glNormal3d(-pVector[N*i + k + N-1]->normal[0],
						   -pVector[N*i + k + N-1]->normal[1], 
						   -pVector[N*i + k + N-1]->normal[2]);
				glVertex3d(pVector[N*i + k + N-1]->Position[0] - del_x*pVector[N*i + k + N-1]->normal[0],
				 		   pVector[N*i + k + N-1]->Position[1] - del_x*pVector[N*i + k + N-1]->normal[1],
				 		   pVector[N*i + k + N-1]->Position[2] - del_x*pVector[N*i + k + N-1]->normal[2]);
			}
		}
	}
	glEnd();
	glDisable(GL_NORMALIZE);
}

