#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------
class Compute {
public:
	/// Creates a compute instance with given geometry and parameter
	Compute (const Geometry* geom, const Parameter* param, const Communicator* comm=0);
	/// Deletes all grids
	~Compute ();

	/// Execute one time step of the fluid simulation (with or without debug info)
        // @ param printInfo print information about current solver state (residual etc.)
	void TimeStep (bool printInfo);

	/// Returns the simulated time in total
	const real_t& GetTime () const;

	/// Returns the pointer to U
	const Grid* GetU () const;
	/// Returns the pointer to V
	const Grid* GetV () const;
	/// Returns the pointer to P
	const Grid* GetP () const;
	/// Returns the pointer to RHS
	const Grid* GetRHS () const;

	/// Computes and returns the absolute velocity
	const Grid* GetVelocity ();
	/// Computes and returns the vorticity
	const Grid* GetVorticity ();
	/// Computes and returns the stream line values
	const Grid* GetStream ();
private:
        // current timestep
	real_t		_t;

        // donor-cell diffusion condition (p. 27)
	real_t		_dtlimit;

        // limit for residual
	real_t		_epslimit;

        // velocities
	Grid*		_u;
	Grid*		_v;

        // pressure
	Grid*		_p;

        // prel. vel
	Grid*		_F;
	Grid*		_G;

        // right-hand side
	Grid*		_rhs;

        // container for interpolating whichever values
	Grid*		_tmp;

	Solver*	   _solver;

	const Geometry*		_geom;
	const Parameter*	_param;
	const Communicator* _comm;

	/// Compute the new velocites u,v
	void NewVelocities	(const real_t& dt);
	/// Compute the temporary velocites F,G
	void MomentumEqu	(const real_t& dt);
	/// Compute the RHS of the poisson equation
	void RHS			(const real_t& dt);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
