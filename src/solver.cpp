#include "solver.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iostream"
#include "math.h"


Solver::Solver(const Geometry *geom){
  _geom = geom;
}

Solver::~Solver() {

}


real_t Solver::localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const {
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  // const real_t center = grid->Cell(it);
  const real_t left = grid->Cell(it.Left());
  const real_t top = grid->Cell(it.Top());
  const real_t down = grid->Cell(it.Down());
  const real_t right = grid->Cell(it.Right());
  return (left + right)/(dx*dx) + (down + top)/(dy*dy) - rhs->Cell(it);
}



SOR::SOR(const Geometry *geom, const real_t &omega)
  : Solver(geom) {
  _omega = omega;
}

SOR::~SOR(){

}

real_t SOR::Cycle(Grid *grid, const Grid *rhs) const {
  InteriorIterator it(_geom);
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  real_t sum_of_squares = 0.0;
  real_t res;
  // for (it.First(); it.Valid(); it.Next()) {
  for( auto it : _geom->getFLUID() ) {
    const real_t center = grid->Cell(it);
    const real_t factor = (dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
    res = localRes(it, grid, rhs);
    // assert(!std::isnan(res));
    sum_of_squares += std::fabs(res - center/factor);
    // sum_of_squares += (res - center/factor)*(res - center/factor);
    grid->Cell(it) = (1-_omega) * center + _omega * factor * res;
  }
  // return sum_of_squares*dx*dy;
  return sum_of_squares/_geom->getFLUID().size();
}
