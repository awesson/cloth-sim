#pragma once

#include <vector>
#include "System.h"

/**
 * Interface for integrators that can step the simulation of a system.
 */
class Integrator
{
public:
    Integrator() { }
    virtual ~Integrator() { }

    /**
     * Step the simulation of the given system by the given timestep.
     * @param sys The system to integrate. It should be integrated starting
     *   from the system's current step.
     * @param dt The length of the time step to integrate.
     */
    virtual void integrate( System sys, double dt ) const = 0;
    // used for storing state vectors locally
    // without allocating memory every time.
    typedef std::vector<Particle*> StateList;
};

/**
 * Uses the basic Euler integration method, x' = x + dx/dt * dt.
 */
class EulerIntegrator : public Integrator
{
public:
	EulerIntegrator(System sys, double dt);
	virtual ~EulerIntegrator();
    virtual void integrate( System sys, double dt ) const;
private:
    mutable StateList state;
    mutable StateList deriv_state;
};

/**
 * Uses the midpoint integration method.
 */
class RK2Integrator : public Integrator
{
public:
	RK2Integrator(System sys, double dt);
	virtual ~RK2Integrator();
	virtual void integrate( System sys, double dt ) const;
private:
        mutable StateList state;
        mutable StateList mid_deriv;
        mutable StateList deriv_state;
};

/**
 * Uses the 4th order Runge-Kutta integration method.
 */
class RK4Integrator : public Integrator
{
public:
	RK4Integrator(System sys, double dt);
	virtual ~RK4Integrator();
    virtual void integrate( System sys, double dt ) const;
private:
    mutable StateList state;
    // the derivatives at each guess
    mutable StateList k1;
    mutable StateList k2;
    mutable StateList k3;
    mutable StateList k4;
};

/**
 * Uses a sympletic Euler integration method, calculating
 * the position explicitly and then the velocity implicitly
 */
class SymplecticEulerIntegrator : public Integrator
{
public:
	SymplecticEulerIntegrator(System sys, double dt);
	virtual ~SymplecticEulerIntegrator();
    virtual void integrate( System sys, double dt ) const;
private:
	mutable StateList state;
	mutable StateList deriv_state;
};

/**
 * Uses Beeman's integration method, calculating
 * the position explicitly and then the velocity implicitly
 */
class BeemanIntegrator : public Integrator
{
public:
	BeemanIntegrator(System sys, double dt);
	virtual ~BeemanIntegrator();
    virtual void integrate( System sys, double dt ) const;
private:
	mutable StateList state;
	mutable StateList next_state;
	mutable StateList deriv_state;
	mutable StateList next_deriv_state;
	mutable StateList prev_deriv_state;
};

/**
 * Uses a verlet integration method, calculating
 * the position explicitly and then the velocity implicitly
 */
class VerletIntegrator : public Integrator
{
public:	
	// sets the previous state for the first time step
	VerletIntegrator(System sys, double dt);
	virtual ~VerletIntegrator();
    virtual void integrate( System sys, double dt ) const;
private:
	mutable StateList state;
	mutable StateList prev_state;
	mutable StateList deriv_state;
};

