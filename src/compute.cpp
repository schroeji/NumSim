#include "compute.hpp"
#include "geometry.hpp"
#include "solver.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "parameter.hpp"

#include <math.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include "assert.h"

Compute::Compute
(
   const Geometry* geom,
   const Parameter* param,
   const Communicator* communicator
)
{
  std::cout << "creating compute for " << communicator->ThreadNum() << std::endl;
  _geom = geom;
  _param = param;
  _epslimit = _param->Eps();
  _comm = communicator;
  _t = 0.0;
  // Berechnung gemäß Skript-Abschnitt SOR
  // _solver = new RedOrBlackSOR(_geom, _param->Omega());

  switch(param->Solver()) {

  case CG_SOLVER:
    _solver = new CG(_geom, _comm);
    std::cout << "CG solver" << std::endl;
    break;

  case MULTIGRID:
    _solver = new MG_Solver(_geom);
    std::cout << "MG solver" << std::endl;
    break;

  case SOR_SOLVER:
    _solver = new SOR(_geom, _param->Omega());
    std::cout << "SOR solver" << std::endl;
    break;

  case REDBLACK:
    _solver = new RedOrBlackSOR(_geom, _param->Omega(), communicator);
    std::cout << "RedOrBlackSOR solver" << std::endl;
    break;

}

  // _solver = new GS_Solver(_geom);
  // _cg_solver = new CG( _geom, _comm );

  // Erzeugen der Gitter evtl fehlen offsets
  real_t dx = _geom->Mesh()[0];
  real_t dy = _geom->Mesh()[1];
  _u = new Grid(_geom, {dx, dy/2.0}, communicator );
  _u->Initialize(0.0);
  _v = new Grid(_geom, {dx/2.0, dy}, communicator );
  _v->Initialize(0.0);
  _p = new Grid(_geom, {dx/2.0, dy/2.0}, communicator );
  _p->Initialize(0.0);
  _F = new Grid(_geom, communicator );
  _F->Initialize(0.0);
  _G = new Grid(_geom, communicator );
  _G->Initialize(0.0);
  _rhs = new Grid(_geom, communicator );
  _rhs->Initialize(0.0);
  _vort = new Grid(_geom, {dx, dy}, communicator);
  _vort->Initialize(0.0);
  _stream = new Grid(_geom, {dx, dy}, communicator);
  _stream->Initialize(0.0);
  std::cout << "created grids for " << communicator->ThreadNum() << std::endl;
  // initial randwerte
  communicator->wait();
  _geom->Update_U(_u);

  _geom->Update_V(_v);
  _geom->Update_U(_F);
  _geom->Update_V(_G);
}



void Compute::TimeStep(bool printinfo) {

  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  const real_t diff_cond = (dx*dx * dy*dy* _param->Re())/(2*dx*dx + 2*dy*dy);

  real_t absMaxValueU = _u->AbsMax();
  //real_t absMaxOverAllProcessesU = _comm->geatherMax( absMaxValueU );
  real_t absMaxValueV = _v->AbsMax();
//  real_t absMaxOverAllProcessesV = _comm->geatherMax( absMaxValueV );
  // const real_t conv_cond = std::min(dx/_u->AbsMax(), dy/_v->AbsMax() );
  const real_t conv_cond = std::min(dx/absMaxValueU, dy/absMaxValueV);
	// real_t dt = std::min(diff_cond, std::min(conv_cond, _param->Dt()));
	real_t dt = std::min(diff_cond*_param->Tau(), std::min(conv_cond * _param->Tau(),_param->Dt()));

   dt = _comm->geatherMin(dt);
  _t += dt;
  if(printinfo) printf("Performing timestep t = %f\n", _t);
  // Randwerte setzen
  if(printinfo) printf("Setting boundary values for u,v...\n");
  _geom->Update_U( _u );
  _geom->Update_V( _v );

  if(printinfo) printf("calculating F and G for inner nodes...\n");
  MomentumEqu(dt);
  // Eigentlich nur einmal nötig
  // aber doppelte Randwerte rechts(_F) und oben (_G) werden vom InteriorIterator verändert
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

 _solver->prepare( _p, _rhs );
  do {
    sum_of_squares = _solver->Cycle(_p, _rhs);
    counter++;
    sum_of_squares = _comm->geatherSum( sum_of_squares );
    _geom->Update_P(_p);
  } while (  std::sqrt(sum_of_squares) > _epslimit  && counter < _param->IterMax());

  if(printinfo) printf("last residual = %f \n", std::sqrt(sum_of_squares));

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


const Grid* Compute::GetVelocity() {
   // Initialize
  _tmp = new Grid(_geom);
   Iterator it = Iterator(_geom);
   for(it.First(); it.Valid(); it.Next() ) {
    multi_real_t pos = {((real_t)it.Pos()[0]) * _geom->Mesh()[0], ((real_t)it.Pos()[1]) * _geom->Mesh()[1]};
    real_t u_val =  _u->Interpolate(pos);
    real_t v_val =  _v->Interpolate(pos);
      _tmp->Cell(it) = sqrt( u_val*u_val + v_val*v_val);
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
  // real_t dx = _geom->Mesh()[0];
  // real_t dy = _geom->Mesh()[1];
  // multi_real_t offset = {dx, dy};
  // _tmp = new Grid(_geom, offset);
  Iterator it = Iterator(_geom);
  for(it.First(); it.Valid(); it.Next()){
    _vort->Cell(it) = _u->dy_r(it) - _v->dx_r(it);
  }
  return _vort;
}


/// Computes and returns the stream line values
const Grid*
Compute::GetStream
(
   void
)
{
  real_t dx = _geom->Mesh()[0];
  real_t dy = _geom->Mesh()[1];
  // multi_real_t offset = {dx, dy};
  // _tmp = new Grid(_geom, offset);
  // _tmp->Initialize(0.0);
  Iterator it(_geom);
  //unten links anfangen und null setzen
  it.First();
  _stream->Cell(it) = 0.0;
  // inkrementelle Berechnung
  for(it.First(); it.Valid(); it.Next()){
    if(it.Pos()[1] == 0) {
      _stream->Cell(it) = _stream->Cell(it.Left()) - dx * _v->Cell(it);
    } else {
      _stream->Cell(it) = _stream->Cell(it.Down()) + dy * _u->Cell(it);
    }
  }

  // Werden nur benutzt um Werte oben links und unten rechts auszulesen
  BoundaryIterator top_left(_geom);
  BoundaryIterator bottom_right(_geom);
  top_left.SetBoundary(3);
  bottom_right.SetBoundary(2);
  top_left.First();
  bottom_right.First();
  real_t add_value = _comm->send_rcv_streamoffset(_stream->Cell(bottom_right),  _stream->Cell(top_left));
  if(! (_comm->isLeft() && _comm->isBottom()) ) {
    for (it.First(); it.Valid(); it.Next()) {
      _stream->Cell(it) += add_value;
    }
  }
  return _stream;
}

/// Compute the new velocites u,v
void
Compute::NewVelocities
(
   const real_t &dt
)
{
  InteriorIterator it(_geom);
  for (it.First(); it.Valid(); it.Next())
  {
    _u->Cell(it) = _F->Cell(it) - dt* _p->dx_r(it);
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
  auto reInverse = 1/re;

  for (it.First(); it.Valid(); it.Next()) {
    real_t A = reInverse* (_u->dxx(it) + _u->dyy(it)) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
    real_t B = reInverse * (_v->dxx(it) + _v->dyy(it)) - _v->DC_udv_x(it, alpha, _u) - _v->DC_vdv_y(it, alpha);
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
  InteriorIterator it(_geom);
  for (it.First(); it.Valid(); it.Next()) {
    _rhs->Cell(it) = (1.0/dt) * (_F->dx_l(it) + _G->dy_l(it));
//     if(std::isnan( _rhs->Cell(it) ) ) {
//       printf(" Nan! Pos: %i;%i Value: %d", it.Pos()[0], it.Pos()[1], _rhs->Cell( it ) );
//     }
    assert(!std::isnan(_rhs->Cell(it)));
  }
}
