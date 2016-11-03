#include "compute.hpp"
#include "geometry.hpp"
#include "solver.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "parameter.hpp"

#include <math.h>
#include <stdio.h>

Compute::Compute (const Geometry *geom, const Parameter *param) {
  _geom = geom;
  _param = param;
  _t = 0.0;
  // Berechnung gemäß Skript-Abschnitt SOR
  _solver = new SOR(_geom, _param->Omega());

  // Erzeugen der Gitter evtl fehlen offsets
  _u = new Grid(_geom);
  _v = new Grid(_geom);
  _p = new Grid(_geom);

  _F = new Grid(_geom);
  _G = new Grid(_geom);
  _rhs = new Grid(_geom);

}
void Compute::TimeStep(bool printinfo) {
  if(printinfo) printf("Performing timestep\n");
  const real_t re = _param->Re();
  const real_t omega = _param->Omega();
  const real_t alpha = _param->Alpha();
  const real_t dt = _param->Dt();
  const real_t eps = _param->Eps();

  // Randwerte setzen
  if(printinfo) printf("WARNING: no boundary values\n");

  if(printinfo) printf("calculating F and G for inner nodes...\n");
  InteriorIterator it(_geom);
  for (it.First(); it.Valid(); it.Next()) {
    real_t A = (1/re) * (_u->dxx(it) + _u->dyy(it)) - _u->DC_udu_x(it, alpha) - _u->DC_vdu_y(it, alpha, _v);
    real_t B = (1/re) * (_v->dxx(it) + _v->dyy(it)) - _u->DC_udv_x(it, alpha, _u) - _v->DC_vdv_y(it, alpha);
    _F->Cell(it) = _u->Cell(it) + dt * A;
    _G->Cell(it) = _v->Cell(it) + dt * B;
  }
  if(printinfo) printf("done\n");

  if(printinfo) printf("calculating right-hand-sides...\n");
  for (it.First(); it.Valid(); it.Next()) {
    _rhs->Cell(it) = (1/dt) * (_F->dx_l(it) + _G->dy_l(it));
  }
  // Lösen der Poissongleichung
  // Update u,v

  _t += 1.0;
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
