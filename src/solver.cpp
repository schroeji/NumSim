#include "solver.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iostream"
#include "math.h"
#include "assert.h"
Solver::Solver(const Geometry *geom){
  _geom = geom;
}

Solver::~Solver() {

}


real_t Solver::localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const {
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  const real_t center = grid->Cell(it);
  const real_t left = grid->Cell(it.Left());
  const real_t top = grid->Cell(it.Top());
  const real_t down = grid->Cell(it.Down());
  const real_t right = grid->Cell(it.Right());
  // std::cout << "it: " << it << std::endl;
  // std::cout << "it Left: " << grid->Cell(it.Left()) << std::endl;
  // std::cout << it.Pos()[0] << ";" << it.Pos()[1]<< std::endl;
  // std::cout << "center: " << center << std::endl;
  // std::cout << "left: " << left << std::endl;
  // std::cout << "top: " << top << std::endl;
  // std::cout << "down: " << down << std::endl;
  // std::cout << "right: " << right << std::endl;
  // std::cout << "rhs: " <<  rhs->Cell(it) << std::endl;
  real_t res1 = rhs->Cell(it) - (left + down - 4*center + right + top)/(dx*dy);
  // real_t res1 = (left + right)/(dx*dx) + (down + top)/(dy*dy) - rhs->Cell(it);
  return res1;
}

SOR::SOR(const Geometry *geom, const real_t &omega)
  : Solver(geom) {
  _omega = omega;
}

SOR::~SOR(){

}

real_t SOR::Cycle(Grid *grid, const Grid *rhs) const {
  Iterator it(_geom);
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  real_t sum_of_squares = 0.0;
  real_t res;
  for (it.First(); it.Valid(); it.Next()) {
    // res = localRes(it, grid, rhs) *  ((dx*dx * dy*dy) / (2 * (dx*dx + dy*dy)));
    res = localRes(it, grid, rhs);
    assert(!std::isnan(res));
    sum_of_squares += res*res ;
    // if(it.Value() > 4)
      // printf("res %u: %f \n", it.Value(), res);
      // std::cout << it << std::endl;
    // std::cout << "sos: " << sum_of_squares << std::endl;
    grid->Cell(it) = grid->Cell(it) - _omega*0.25*(dx*dy) * res;
    // grid->Cell(it) = (1-_omega)*grid->Cell(it) + _omega * ((dx*dx * dy*dy) / (2 * (dx*dx + dy*dy))) * res;
    // grid->Cell(it) = (1-_omega)*grid->Cell(it) + _omega * res;
  }
  return sum_of_squares;
}
