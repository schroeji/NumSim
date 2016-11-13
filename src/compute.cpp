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
  _geom->Update_U(_F);
  _geom->Update_V(_G);
}
void Compute::TimeStep(bool printinfo) {
  if(printinfo) printf("Performing timestep\n");
  const real_t dt = _param->Dt();
  _t += dt;
  // Randwerte setzen
  if(printinfo) printf("Setting boundary values for u,v...\n");
  _geom->Update_U(_u);
  _geom->Update_V(_v);
  if(printinfo) printf("calculating F and G for inner nodes...\n");
  MomentumEqu(dt);
  if(printinfo) printf("done\n");

  if(printinfo) printf("calculating right-hand-sides...\n");
  RHS(dt);
  if(printinfo) printf("done\n");

  // Lösen der Poissongleichung
  if(printinfo) printf("solving with eps = %f \n", _epslimit);

  index_t counter = 0;
  real_t sum_of_squares;
  do {
    _geom->Update_P(_p);
    sum_of_squares = _solver->Cycle(_p, _rhs);
    std::cout << sum_of_squares << std::endl;
    counter++;
  // } while (counter < 5);
  } while ( sqrt( sum_of_squares/(_geom->Size()[0] * _geom->Size()[1]) ) > _epslimit  && counter < _param->IterMax());
  _geom->Update_P(_p);
  if(printinfo) printf("Convergence after %i iterations\n", counter);
  // Update u,v
  NewVelocities(dt);
}



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


/// Computes and returns the absolute velocity
const Grid*
Compute::GetVelocity
(
   void
)
{
  _tmp = new Grid(_geom);
  Iterator it(_geom);
  for(it.First(); it.Valid(); it.Next()) {
    multi_real_t pos = {(real_t) 1/it.Pos()[0], (real_t) 1/it.Pos()[1]};
    // std::cout << pos[0] << std::endl;
    _tmp->Cell(it) = sqrt(_u->Interpolate(pos)*_u->Interpolate(pos) + _v->Interpolate(pos)* _v->Interpolate(pos));
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

}


/// Computes and returns the stream line values
const Grid*
Compute::GetStream
(
   void
)
{

}

/// Compute the new velocites u,v
void
Compute::NewVelocities
(
   const real_t &dt
)
{
  Iterator it(_geom);
  for (it.First(); it.Valid(); it.Next()){
    _u->Cell(it) = _F->Cell(it) - dt* _p->Cell(it);
    _v->Cell(it) = _G->Cell(it) - dt* _p->Cell(it);
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
    // std::cout << "u_udux: " << _u->DC_udu_x(it, alpha) << std::endl;
    // std::cout << "u_vduy: " << _u->DC_vdu_y(it, alpha, _v) << std::endl;
    // std::cout << "u_dxx: " << _u->dxx(it) << std::endl;
    // std::cout << "u_dyy: " << _u->dyy(it) << std::endl;
    _F->Cell(it) = _u->Cell(it) + dt * A;
    _G->Cell(it) = _v->Cell(it) + dt * B;
  }
}


/// Compute the RHS of the poisson equation
void
Compute::RHS
(
   const real_t &dt
)
{
  Iterator it(_geom);
  for (it.First(); it.Valid(); it.Next()) {
    _rhs->Cell(it) = (1/dt) * (_F->dx_l(it) + _G->dy_l(it));
    assert(!std::isnan(_rhs->Cell(it)));
  }
}
