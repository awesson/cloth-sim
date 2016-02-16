#include "integrator.h"

EulerIntegrator::EulerIntegrator(System sys, double dt)
{
	int size = sys.size();
	
	state.resize( size );
	
	for(int i = 0; i < size; ++i){
		state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
	}
}

EulerIntegrator::~EulerIntegrator()
{
	state.clear();
	deriv_state.clear();
}


/**
 * Uses the basic Euler integration method, x' = x + dx/dt * dt.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void EulerIntegrator::integrate( System sys, double dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize( size );
    deriv_state.resize( size );

    // get the current state
    sys.get_state( state );

    // compute the current derivative
    sys.deriv_eval( deriv_state );

    // update the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += deriv_state[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt;
        }
    }

    // set the updated state
    sys.set_state( state );
}

RK2Integrator::RK2Integrator(System sys, double dt)
{
	int size = sys.size();
	
	state.resize( size );
	
	for(int i = 0; i < size; ++i)
		state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
}

RK2Integrator::~RK2Integrator()
{
	state.clear();
	mid_deriv.clear();
	deriv_state.clear();
}

/**
 * Uses the midpoint integration method.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void RK2Integrator::integrate( System sys, double dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize( size );
    mid_deriv.resize( size );
    deriv_state.resize( size );

    // get the current state
    sys.get_state( state );

    // compute the current derivative
    sys.deriv_eval( deriv_state );

    // get midpoint state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += deriv_state[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get derivative at the midpoint
    sys.set_state( state );
    sys.deriv_eval( mid_deriv );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] -= deriv_state[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] -= deriv_state[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get full point state using midpoint derivative
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += mid_deriv[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += mid_deriv[i]->deriv_velocity[j] * dt;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( state );
}


RK4Integrator::RK4Integrator(System sys, double dt)
{
	int size = sys.size();
	
	state.resize( size );

	for(int i = 0; i < size; ++i)
		state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
}

RK4Integrator::~RK4Integrator()
{
	state.clear();
    k1.clear();
    k2.clear();
    k3.clear();
    k4.clear();
}


/**
 * Uses the 4th order Runge-Kutta integration method.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void RK4Integrator::integrate( System sys, double dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize( size );
    k1.resize( size );
    k2.resize( size );
    k3.resize( size );
    k4.resize( size );

    // get the current state
    sys.get_state( state );

    // compute the current derivative
    sys.deriv_eval( k1 );

    // get midpoint state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += k1[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] += k1[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get derivative at the midpoint
    sys.set_state( state );
    sys.deriv_eval( k2 );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] -= k1[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] -= k1[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get second midpoint state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += k2[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] += k2[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get derivative at the new midpoint
    sys.set_state( state );
    sys.deriv_eval( k3 );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] -= k2[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] -= k2[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get final state using new midpoint derivative
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += k3[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += k3[i]->deriv_velocity[j] * dt;
        }
    }

    // get derivative at the guess for the final state
    sys.set_state( state );
    sys.deriv_eval( k4 );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] -= k3[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] -= k3[i]->deriv_velocity[j] * dt;
        }
    }

    // get final state using all 4 derivatives
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += k1[i]->deriv_position[j] * dt/6.f + k2[i]->deriv_position[j] * dt/3.f
                                   + k3[i]->deriv_position[j] * dt/3.f + k4[i]->deriv_position[j] * dt/6.f;
            state[i]->Velocity[j] += k1[i]->deriv_velocity[j] * dt/6.f + k2[i]->deriv_velocity[j] * dt/3.f
                                   + k3[i]->deriv_velocity[j] * dt/3.f + k4[i]->deriv_velocity[j] * dt/6.f;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( state );
}

SymplecticEulerIntegrator::SymplecticEulerIntegrator(System sys, double dt)
{
	int size = sys.size();
	
	state.resize( size );
	
	for(int i = 0; i < size; ++i){
		state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
	}
}

SymplecticEulerIntegrator::~SymplecticEulerIntegrator()
{
	state.clear();
    deriv_state.clear();
}


/**
 * Uses a symplectic Euler integration method. First the position is
 * calculated explicitly and the the velocity is calculated implicitly.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void SymplecticEulerIntegrator::integrate( System sys, double dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize(size);
    deriv_state.resize(size);

    // get the current state (pos, t)
    sys.get_state( state );

    // compute the current derivative (pos')
    sys.deriv_eval( deriv_state );

    // update the velocity explicitly.
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt;
        }
    }

    // update the position implicitly with the velocity at t + dt
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += state[i]->Velocity[j] * dt;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( state );
}

BeemanIntegrator::BeemanIntegrator(System sys, double dt)
{
	int size = sys.size();
	
    state.resize(size);
	next_state.resize(size);
	prev_deriv_state.resize(size);

	sys.deriv_eval( prev_deriv_state );
	
	for(int i = 0; i < size; ++i){
		state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
		next_state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
	}
}

BeemanIntegrator::~BeemanIntegrator()
{
	state.clear();
	next_state.clear();
	prev_deriv_state.clear();
	next_deriv_state.clear();
    deriv_state.clear();
}

/**
 * Uses Beeman's integration method. First the position is
 * calculated explicitly and the the velocity is calculated implicitly.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void BeemanIntegrator::integrate( System sys, double dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize(size);
	next_state.resize(size);
	prev_deriv_state.resize(size);
	next_deriv_state.resize(size);
    deriv_state.resize(size);

    // get the current state (pos, t)
    sys.get_state( state );
	sys.get_state( next_state );

    // compute the current derivative (pos')
    sys.deriv_eval( deriv_state );

    // update the position with the position at t - dt
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            next_state[i]->Position[j] = state[i]->Position[j] + state[i]->Velocity[j]*dt
 									   + 2.0/3.0*deriv_state[i]->deriv_velocity[j]*dt*dt
									   - 1.0/6.0*prev_deriv_state[i]->deriv_velocity[j]*dt*dt;
        }
    }

    // update the velocity implicitly from the position.
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            next_state[i]->Velocity[j] = state[i]->Velocity[j] + 1.5*deriv_state[i]->deriv_velocity[j]*dt
 									   - .5*prev_deriv_state[i]->deriv_velocity[j]*dt;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( next_state );

	sys.deriv_eval( next_deriv_state );
	// update the velocity implicitly from the position.
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            next_state[i]->Velocity[j] = state[i]->Velocity[j] + 1.0/3.0*next_deriv_state[i]->deriv_velocity[j]*dt
 									   + 5.0/6.0*deriv_state[i]->deriv_velocity[j]*dt - 1.0/6.0*prev_deriv_state[i]->deriv_velocity[j]*dt;
        }
    }

	sys.deriv_eval( prev_deriv_state );

	// set the state to (pos + pos' * dt, t + dt)
    sys.set_state( next_state );
	
}

VerletIntegrator::VerletIntegrator(System sys, double dt)
{
	int size = sys.size();

    if (size == 0)
        return;

    state.resize(size);
	prev_state.resize(size);
    deriv_state.resize(size);

	for(int i = 0; i < size; ++i){
		state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
		prev_state[i] = new Particle(make_vector(0.0, 0.0, 0.0), 1.0, 1);
	}

    // get the current state (pos, t)
    sys.get_state( state );
	// get the current state as the previous state for the next iteration
	sys.get_state( prev_state );

    // compute the current derivative (pos')
    sys.deriv_eval( deriv_state );

    // update the position
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] += state[i]->Velocity[j]*dt + 0.5*deriv_state[i]->deriv_velocity[j]*dt*dt;
		}
    }

    // update the velocity
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j]*dt;
        }
    }

    // update the system
    sys.set_state( state );
}

VerletIntegrator::~VerletIntegrator()
{
	state.clear();
	prev_state.clear();
    deriv_state.clear();
}

/**
 * Uses Verlet integration. First the position is
 * calculated explicitly and the the velocity is calculated implicitly.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void VerletIntegrator::integrate( System sys, double dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize(size);
	prev_state.resize(size);
    deriv_state.resize(size);

    // get the current state (pos, t)
    sys.get_state( state );

    // compute the current derivative (pos')
    sys.deriv_eval( deriv_state );

    // update the position using the position at t - dt
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 3; ++j){
            state[i]->Position[j] = 2*state[i]->Position[j] - prev_state[i]->Position[j] + deriv_state[i]->deriv_velocity[j]*dt*dt;
       	}
    }

	// store the current position as the previous position for the next iteration
	sys.get_state( prev_state );

    // update the velocity implicitly from the positions
    for(int i = 0; i < size; ++i){
            for(int j = 0; j < 3; ++j){
                state[i]->Velocity[j] = (state[i]->Position[j] - prev_state[i]->Position[j])/dt;
            }
        }

    // update the system
    sys.set_state( state );
}

