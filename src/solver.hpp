#include "typedef.hpp"
#include "comm.hpp"
//------------------------------------------------------------------------------
#ifndef __SOLVER_HPP
#define __SOLVER_HPP
//------------------------------------------------------------------------------

/** abstract base class for an iterative solver
 */
class Solver {
public:
	/// Constructor of the abstract Solver class
	Solver (const Geometry* geom);
	/// Destructor of the Solver Class
	virtual ~Solver();

	/// This function must be implemented in a child class
        // @param [in][out] grid current values
        // @param [in]      rhs  right hand side values
        // @returns accumulated residual
	virtual real_t Cycle (Grid *grid, const Grid* rhs) const = 0;
  virtual void prepare( Grid* grid, const Grid* rhs );
protected:
	const Geometry* _geom;

	/// Returns the residual at [it] for the pressure-Poisson equation
	real_t localRes (const Iterator& it, const Grid *grid, const Grid* rhs) const;
};

//------------------------------------------------------------------------------

/** concrete SOR solver
 */
class SOR : public Solver {
public:
	/// Constructs an actual SOR solver
	SOR (const Geometry* geom, const real_t& omega);
	/// Destructor
	~SOR();

	/// Returns the total residual and executes a solver cycle
        // @param grid current pressure values
        // @param rhs right hand side
	real_t Cycle (Grid* grid, const Grid* rhs) const;
protected:
	real_t _omega;
};
//------------------------------------------------------------------------------

/** concrete Red or Balck SOR solver
 */
class RedOrBlackSOR : public SOR {
public:
	/// Constructs an actual SOR solver
	RedOrBlackSOR (const Geometry* geom, const real_t& omega);
	/// Destructor
	~RedOrBlackSOR();
	real_t RedCycle (Grid* grid, const Grid* rhs) const;
	real_t BlackCycle (Grid* grid, const Grid* rhs) const;
};

//------------------------------------------------------------------------------


class CG : public Solver
{
public:
  CG( const Geometry* geom, const Communicator* comm );
  ~CG( );

  void prepare( Grid* grid, const Grid* rhs );
  real_t Cycle( Grid* grid, const Grid *rhs ) const;
private:
  Grid* _residuum;
  Grid* _direction;
  const Communicator* _comm;
};

//------------------------------------------------------------------------------

class MG_Solver : public Solver {
public:
  MG_Solver(const Geometry* geom);
  ~MG_Solver();

  real_t Cycle(Grid *grid, const Grid *rhs) const;
  real_t collectResidual(Grid* grid, const Grid* rhs) const;
  void Iteration(Grid * grid, const Grid *rhs) const;
  real_t Smooth(Grid *grid, const Grid *rhs) const;
  // returns a pointer to the restricted residual
  Grid* Restrict(const Grid *fine_grid, const Grid *fine_rhs, Grid* coarse_residual) const;
  // returns a pointer to the interpolated grid
  void ProlongAndAdd(Grid* fine_grid, const Grid *coarse_grid) const;

private:
  Solver* _smoother;
};
//------------------------------------------------------------------------------
#endif // __SOLVER_HPP
