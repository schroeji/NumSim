#include "solver.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iostream"


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
  real_t res1;

  const bool leftBoundaryNeighbor = it.Pos()[0] == 1;
  const bool rightBoundaryNeighbor = it.Pos()[0] == _geom->Size()[0];
  const bool underBoundaryNeighbor = it.Pos()[1] == 1;
  const bool topBoundaryNeighbor = it.Pos()[1] == _geom->Size()[1];

  //std::cout << it.Pos()[0] << " " << it.Pos()[1] << std::endl;
  if( leftBoundaryNeighbor )
  {
//     std::cout << it.Pos()[0] << " " << it.Pos()[1] << std::endl;
     if( topBoundaryNeighbor )
     {
        res1 = rhs->Cell(it) - ( -center + right )/(dx*dx) - ( down - center )/( dy * dy );
     }
     else if( underBoundaryNeighbor )
     {
        res1 = rhs->Cell(it) - ( -center + right )/(dx*dx) - ( -center + top )/( dy * dy );
     }
     else
     {
        res1 = rhs->Cell(it) - ( -center + right )/(dx*dx) - ( down -2.0*center + top )/( dy * dy );
     }
  }
  else if( rightBoundaryNeighbor )
  {
     //std::cout << it.Pos()[0] << " " << it.Pos()[1] << std::endl;
     if( topBoundaryNeighbor )
     {
        res1 = rhs->Cell(it) - ( left - center )/(dx*dx) - ( down - center )/( dy * dy );
     }
     else if( underBoundaryNeighbor )
     {
        res1 = rhs->Cell(it) - ( left - center )/(dx*dx) - ( -center + top )/( dy * dy );
     }
     else
     {
        res1 = rhs->Cell(it) - ( left - center )/(dx*dx) - ( down -2.0*center + top )/( dy * dy );
     }
  }
  else if( topBoundaryNeighbor )
  {
     //std::cout << it.Pos()[0] << " " << it.Pos()[1] << std::endl;
     res1 = rhs->Cell(it) - ( left  - 2.0*center + right )/(dx*dx) - ( down - center )/( dy * dy );
  }
  else if( underBoundaryNeighbor )
  {
     //std::cout << it.Pos()[0] << " " << it.Pos()[1] << std::endl;
     res1 = rhs->Cell(it) - ( left  - 2.0*center + right )/(dx*dx) - ( top - center )/( dy * dy );
  }
  else
  {
    res1 = rhs->Cell(it) - ( left  - 2*center + right )/(dx*dx) - ( -2*center + top + down )/( dy * dy );
  }
  //real_t res2 = (left + right)/(dx*dx) + (down - top)/(dy*dy) - rhs->Cell(it);
  return res1;
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
  for (it.First(); it.Valid(); it.Next())
  {
      res = localRes( it, grid, rhs );
      //   const real_t res = localRes(it, grid, rhs);
      //   // std::cout << "cycle2:" << grid->Data()[0] << std::endl;
      //std::cout << "res in cycle: " << res <<std::endl;
      sum_of_squares += res*res;
//      std::cout << "sum_of_suares " << sum_of_squares <<std::endl;
      grid->Cell(it) = grid->Cell(it) - _omega * ((dx*dx * dy*dy) / (2 * (dx*dx + dy*dy))) * res;
  }
  std::cout << "sum_of_suares " << sum_of_squares <<std::endl;
  return sum_of_squares;
}
