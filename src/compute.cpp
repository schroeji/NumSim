#include "compute.hpp"
#include "geometry.hpp"
#include "solver.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "parameter.hpp"

#include <math.h>
#include <iostream>
#include <stdio.h>
#include "assert.h"

Compute::Compute (const Geometry *geom, const Parameter *param) {
  _geom = geom;
  _param = param;
  _epslimit = _param->Eps();
  _t = 0.0;
  // Berechnung gemäß Skript-Abschnitt SOR
  _solver = new SOR(_geom, _param->Omega());

  // Erzeugen der Gitter evtl fehlen offsets
  real_t dx = _geom->Mesh()[0];
  real_t dy = _geom->Mesh()[1];
  _u = new Grid(_geom, {dx, dy/2.0} );
  _u->Initialize(0.0);
  _v = new Grid(_geom, {dx/2.0, dy} );
  _v->Initialize(0.0);
  _p = new Grid(_geom, {dx/2.0, dy/2.0});
  _p->Initialize(0.0);
  _F = new Grid(_geom);
  _F->Initialize(0.0);
  _G = new Grid(_geom);
  _G->Initialize(0.0);
  _rhs = new Grid(_geom);
  _rhs->Initialize(0.0);

  // Randwerte von F
}


void Compute::TimeStep(bool printinfo) {
  _geom->Update_U(_u);
  _geom->Update_V(_v);
  real_t dt = _param->Tau()*std::fmax(_geom->Mesh()[0],_geom->Mesh()[1])/std::fmax(_u->AbsMax(),_v->AbsMax());
	real_t dt2 = _param->Tau()*_param->Re()/2* (_geom->Mesh()[1]*_geom->Mesh()[1]*_geom->Mesh()[0]*_geom->Mesh()[0]);
	dt2 = dt2/(_geom->Mesh()[1]*_geom->Mesh()[1]+_geom->Mesh()[0]*_geom->Mesh()[0]);
	dt = std::min(dt2,std::min(dt,_param->Dt()));

  // const real_t dt = 0.001;
  _t += dt;
  if(printinfo) printf("Performing timestep t = %f\n", _t);
  // Randwerte setzen
  if(printinfo) printf("Setting boundary values for u,v...\n");
  // _geom->Update_U(_u);
  // _geom->Update_V(_v);
  if(printinfo) printf("calculating F and G for inner nodes...\n");
  MomentumEqu(dt);
  // Eigentlich nur einmal nötig
  _geom->Update_U(_F);
  _geom->Update_V(_G);
  if(printinfo) printf("done\n");

  if(printinfo) printf("calculating right-hand-sides...\n");
  RHS(dt);
  if(printinfo) printf("done\n");

  // Lösen der Poissongleichung
  if(printinfo) printf("solving with eps = %f \n", _epslimit);

  index_t counter = 0;
  real_t sum_of_squares;
  // Iterator it (_geom);
  // for (it.First();it.Valid();it.Next()){
  //   if (it.Pos()[1] == 128) {
  //     std::cout << "uvp" << std::endl;
  //     std::cout << _u->Cell(it) << std::endl;
  //     std::cout << _v->Cell(it) << std::endl;
  //     std::cout << _p->Cell(it) << std::endl;
  //     // std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  //     // std::cout << _rhs->Cell(it) << std::endl;
  //     // std::cout << "F_dl: " << _F->dx_l(it) << std::endl;
  //   }
  // }
  do {
    _geom->Update_P(_p);
    sum_of_squares = _solver->Cycle(_p, _rhs);
    // std::cout << sum_of_squares << std::endl;
    // std::cout << sqrt( sum_of_squares/(_geom->Size()[0] * _geom->Size()[1]) ) << std::endl;
    counter++;
  } while (  sum_of_squares > _epslimit  && counter < _param->IterMax());

  if(printinfo) printf("last residual = %f \n", sum_of_squares/(_geom->Size()[0] * _geom->Size()[1]) );
  // InteriorIterator ir (_geom);
  // for (ir.First();ir.Valid();ir.Next()){
    // if (ir.Pos()[1] == 128)
      // std::cout << _p->Cell(ir) << std::endl;
  // }
  if(printinfo) printf("Convergence after %i iterations\n", counter);
  // Update u,v
  NewVelocities(dt);
}


// void Compute::TimeStep(bool printInfo) {
// 	// Compute boundary Values
// 	// eigentlich erst nach dt, aber dann geht der erste Zeitschritt zu lang.
// 	_geom->Update_U(_u);
// 	_geom->Update_V(_v);
// 	_geom->Update_P(_p);

// 	// Compute dt
// 	//real_t dt = _param->Dt();
// 	// real_t dt = _param->Tau()*std::fmax(_geom->Mesh()[0],_geom->Mesh()[1])/std::fmax(_u->AbsMax(),_v->AbsMax());
// 	// real_t dt2 = _param->Tau()*_param->Re()/2* (_geom->Mesh()[1]*_geom->Mesh()[1]*_geom->Mesh()[0]*_geom->Mesh()[0]);
// 	// dt2 = dt2/(_geom->Mesh()[1]*_geom->Mesh()[1]+_geom->Mesh()[0]*_geom->Mesh()[0]);
// 	// dt = std::min(dt2,std::min(dt,_param->Dt()));
//   real_t dt = 0.004;

// 	// Compute F, G
// 	MomentumEqu(dt);
// 	//std::cin.ignore();

// 	// Compute RHS
// 	RHS(dt);

// 	// Compute p
// 	real_t res = 10000000;
// 	index_t i = 0;
// 	while(res >_param->Eps() && i < _param->IterMax() ) {
// 		res = _solver->Cycle(_p,_rhs);
// 		_geom->Update_P(_p);
// 		i++;
// 	}

// 	// Compute u,v
// 	NewVelocities(dt);

// 	// Next timestep
// 	_t += dt;

// 	// Print info
// 	if (printInfo) {
//     std::cout << "_t: " << _t << "  \tres: " << std::scientific << res << "\t progress: " << std::fixed << _t/_param->Tend()*100 << "%" << std::endl;
// 	}
// }


const real_t&
Compute::GetTime
(
   void
) const
{
   return _t;
}


Compute::~Compute () {
  // Loeschen der erzeugten Gitter etc
}

const Grid* Compute::GetU() const {
  return _u;
}

const Grid* Compute::GetV() const {
  return _v;
}

const Grid* Compute::GetP() const {
  return _p;
}

const Grid* Compute::GetRHS() const {
  return _rhs;
}


const Grid* Compute::GetVelocity() {
	// Initialize
  return _u;
	Iterator it = Iterator(_geom);

	// Cycle through all cells
	for(it.First(); it.Valid(); it.Next() ) {
		_tmp->Cell(it) = sqrt(( _u->Cell(it)*_u->Cell(it) + _v->Cell(it)*_v->Cell(it)));
	}

	return _tmp;
}


/// Computes and returns the vorticity
const Grid*
Compute::GetVorticity
(
   void
)
{
  return _tmp;
}


/// Computes and returns the stream line values
const Grid*
Compute::GetStream
(
   void
)
{
  return _tmp;
}

/// Compute the new velocites u,v
void
Compute::NewVelocities
(
   const real_t &dt
)
{
  InteriorIterator it(_geom);
  for (it.First(); it.Valid(); it.Next()){
    // _u->Cell(it) = _F->Cell(it) - dt* _p->Cell(it);
    _u->Cell(it) = _F->Cell(it) - dt* _p->dx_r(it);
    // _v->Cell(it) = _G->Cell(it) - dt* _p->Cell(it);
    _v->Cell(it) = _G->Cell(it) - dt* _p->dy_r(it);
  }
}


/// Compute the temporary velocites F,G
void
Compute::MomentumEqu
(
   const real_t &dt
)
{
  const real_t alpha = _param->Alpha();
  const real_t re = _param->Re();
  InteriorIterator it(_geom);
  for (it.First(); it.Valid(); it.Next()) {
    real_t A = (1/re) * (_u->dxx(it) + _u->dyy(it)) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
    real_t B = (1/re) * (_v->dxx(it) + _v->dyy(it)) - _v->DC_udv_x(it, alpha, _u) - _v->DC_vdv_y(it, alpha);
    // std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
    // std::cout << "A: " << A << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "v_udvx: " << _v->DC_udv_x(it, alpha, _u) << std::endl;
    // std::cout << "v_vdvy: " << _v->DC_vdv_y(it, alpha) << std::endl;
    // std::cout << "v_dxx: " << _v->dxx(it) << std::endl;
    // std::cout << "u_udux: " << _u->DC_udu_x(it, alpha) << std::endl;
    // std::cout << "u_vduy: " << _u->DC_vdu_y(it, alpha, _v) << std::endl;
    // std::cout << "u_dxx: " << _u->dxx(it) << std::endl;
    // std::cout << "u_dyy: " << _u->dyy(it) << std::endl;
    // std::cout << "u_dy_r: " << _u->dy_r(it) << std::endl;
    // std::cout << "u_dy_l: " << _u->dy_l(it) << std::endl;
    // std::cout << "u_top: " << _u->Cell(it.Top()) << std::endl;
    // std::cout << "u_it: " << _u->Cell(it) << std::endl;
    _F->Cell(it) = _u->Cell(it) + dt * A;
    _G->Cell(it) = _v->Cell(it) + dt * B;
  }
}
// void Compute::MomentumEqu(const real_t& dt) {
// 	// Initialize interior iterator
// 	InteriorIterator it = InteriorIterator(_geom);

// 	// Cycle through all interior cells
// 	while ( it.Valid() ) {
// 		// Get current velocities
// 		const real_t u = _u->Cell(it);
// 		const real_t v = _v->Cell(it);

// 		// Get parameter
// 		const real_t RE_inv = 1.0/_param->Re();
// 		const real_t alpha = _param->Alpha();

// 		// Calculate temporary velocietes
// 		real_t A = RE_inv * ( _u->dxx(it) + _u->dyy(it) ) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
// 		real_t B = RE_inv * ( _v->dxx(it) + _v->dyy(it) ) - _v->DC_vdv_y(it, alpha) - _v->DC_udv_x(it, alpha, _u);

// 		// Update new temporary velocities
// 		_F->Cell(it) = u + dt*A;
// 		_G->Cell(it) = v + dt*B;

// 		// Next cell
// 		it.Next();
// 	}

// 	_geom->Update_U(_F);
// 	_geom->Update_V(_G);
// }


/// Compute the RHS of the poisson equation
void
Compute::RHS
(
   const real_t &dt
)
{
  InteriorIterator it(_geom);
  for (it.First(); it.Valid(); it.Next()) {
    _rhs->Cell(it) = (1.0/dt) * (_F->dx_l(it) + _G->dy_l(it));
    assert(!std::isnan(_rhs->Cell(it)));
  }
}
