#include "solver.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "grid.hpp"

Solver::Solver(const Geometry *geom){
  _geom = geom;
}

real_t Solver::localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const {
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  const real_t center = grid->Cell(it);
  const real_t left = grid->Cell(it.Left());
  const real_t top = grid->Cell(it.Top());
  const real_t down = grid->Cell(it.Down());
  const real_t right = grid->Cell(it.Right());

  return rhs->Cell(it) - (left + down - 4*center + right + top)/(dx*dy);
}

SOR::SOR(const Geometry *geom, const real_t &omega)
  : Solver(geom) {
  _omega = omega;
}

real_t SOR::Cycle(Grid *grid, const Grid *rhs) const {
  Iterator it(_geom);
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  real_t sum_of_squares = 0;

  for (it.First(); it.Valid(); it.Next()) {
    const real_t res = localRes(it, grid, rhs);
    sum_of_squares += res*res;
    grid->Cell(it) = grid->Cell(it) - _omega * ((dx*dx * dy*dy) / (2 * (dx*dx + dy*dy))) * res;
  }
  return sum_of_squares;
}
