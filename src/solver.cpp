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
  // const real_t center = grid->Cell(it);
  const real_t left = grid->Cell(it.Left());
  const real_t top = grid->Cell(it.Top());
  const real_t down = grid->Cell(it.Down());
  const real_t right = grid->Cell(it.Right());
  // return (left + right)/(dx*dx) + (down + top)/(dy*dy) - rhs->Cell(it);
  return rhs->Cell(it) - grid->dxx(it) - grid->dyy(it);
}

//------------------------------------------------------------------------------

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
  for (it.First(); it.Valid(); it.Next()) {
    const real_t center = grid->Cell(it);
    const real_t factor = (dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
    res = localRes(it, grid, rhs);
    // assert(!std::isnan(res));
    sum_of_squares += fabs(res - center/factor);
    // sum_of_squares += (res - center/factor)*(res - center/factor);
    grid->Cell(it) = (1-_omega) * center + _omega * factor * res;
    assert(!std::isnan(grid->Cell(it)));
  }
  return sum_of_squares*dx*dy;
}

//------------------------------------------------------------------------------
GS_Solver::GS_Solver(const Geometry *geom)
  : Solver(geom) {
}

GS_Solver::~GS_Solver() {
}

real_t GS_Solver::Cycle(Grid *grid, const Grid *rhs) const {
  InteriorIterator it(_geom);
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  real_t sum_of_squares = 0.0;
  real_t res;
  for (it.First(); it.Valid(); it.Next()) {
    const real_t center = grid->Cell(it);
    const real_t factor = (dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
    res = -localRes(it, grid, rhs);
    // assert(!std::isnan(res));
    sum_of_squares += fabs(res - center/factor);
    // sum_of_squares += fabs(res);
    // sum_of_squares += (res - center/factor)*(res - center/factor);
    grid->Cell(it) = grid->Cell(it) + factor * res;
    assert(!std::isnan(grid->Cell(it)));
  }
  return sum_of_squares*dx*dy;
}

//------------------------------------------------------------------------------

RedOrBlackSOR::RedOrBlackSOR
(
   const Geometry* geom,
   const real_t& omega
) : SOR( geom, omega )
{

}


RedOrBlackSOR::~RedOrBlackSOR
(
   void
)
{

}



real_t RedOrBlackSOR::RedCycle
(
   Grid* grid,
   const Grid* rhs
) const
{
   InteriorIterator it( _geom );
   const real_t dx = _geom->Mesh()[0];
   const real_t dy = _geom->Mesh()[1];
   real_t sum_of_squares = 0.0;
   real_t res;
	it.First();
   for( it.Next(); it.Valid(); it.Next() )
   {
     const real_t center = grid->Cell(it);
     const real_t factor = (dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
     res = localRes(it, grid, rhs);
     // assert(!std::isnan(res));
     sum_of_squares += fabs(res - center/factor);
     // sum_of_squares += (res - center/factor)*(res - center/factor);
     grid->Cell(it) = (1-_omega) * center + _omega * factor * res;
     it.Next();
   }
   return sum_of_squares*dx*dy;
}


real_t RedOrBlackSOR::BlackCycle
(
   Grid* grid,
   const Grid* rhs
) const
{
   InteriorIterator it( _geom );
   const real_t dx = _geom->Mesh()[0];
   const real_t dy = _geom->Mesh()[1];
   real_t sum_of_squares = 0.0;
   real_t res;
   for( it.First(); it.Valid(); it.Next() )
   {
     const real_t center = grid->Cell(it);
     const real_t factor = (dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
     res = localRes(it, grid, rhs);
     // assert(!std::isnan(res));
     sum_of_squares += fabs(res - center/factor);
     // sum_of_squares += (res - center/factor)*(res - center/factor);
     grid->Cell(it) = (1-_omega) * center + _omega * factor * res;
     it.Next();
   }
   return sum_of_squares*dx*dy;
}

//------------------------------------------------------------------------------

MG_Solver::MG_Solver(const Geometry *geom)
  : Solver(geom) {
  _smoother = new GS_Solver(geom);
}

MG_Solver::~MG_Solver() {

}

real_t MG_Solver::Cycle(Grid* grid, const Grid *rhs) const {
  Iteration(grid, rhs);

  // calc residual
  return collectResidual(grid, rhs);
}

real_t MG_Solver::collectResidual(Grid* grid, const Grid* rhs) const {
  InteriorIterator it(_geom);
  const real_t dx = _geom->Mesh()[0];
  const real_t dy = _geom->Mesh()[1];
  real_t sum_of_squares = 0.0;
  real_t res;
  for (it.First(); it.Valid(); it.Next()) {
    const real_t center = grid->Cell(it);
    const real_t factor = (dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
    res = localRes(it, grid, rhs);
    // assert(!std::isnan(res));
    sum_of_squares += fabs(res - center/factor);
    // sum_of_squares += fabs(res);
    // sum_of_squares += (res - center/factor)*(res - center/factor);
    // grid->Cell(it) = factor * res;
    assert(!std::isnan(grid->Cell(it)));
  }
  return sum_of_squares*dx*dy;
}

void MG_Solver::Iteration(Grid* grid, const Grid *rhs) const {
  // smooth
  // std::cout << "smoothing" << std::endl;
  real_t res = 0.0;
  res = Smooth(grid, rhs);
  // std::cout << "after smooth1:" << res << std::endl;
  if(_geom->Size()[0] <= 8) {
    // solve
    // SOR *solver = new SOR(_geom, 1.7);
    GS_Solver *solver = new GS_Solver(_geom);
    index_t counter = 0;
    // std::cout << "solving... res:" << collectResidual(grid,rhs) << std::endl;
    do {
      // _geom->Update_P(grid);
      res = solver->Cycle(grid, rhs);
      counter++;
      // std::cout << "Res: " << res<< std::endl;
      // std::cout << "collected res:" << collectResidual(grid,rhs) << std::endl;
    } while (counter < 20);
    // std::cout << "solved last res:" << collectResidual(grid,rhs) << std::endl;
    return;
  }
  else {
    // residual and restrict
    // std::cout << "restricting" << std::endl;
    Geometry* temp_geom = new Geometry(_geom);
    temp_geom->setSize( {_geom->Size()[0]/2, _geom->Size()[1]/2} );
    const Geometry* coarse_geom = const_cast<const Geometry*>(temp_geom);
    const multi_real_t coarse_offset = {coarse_geom->Mesh()[0]/2, coarse_geom->Mesh()[1]/2};
    Grid* coarse_residual = new Grid(coarse_geom, coarse_offset, coarse_geom->comm());
    coarse_residual->Initialize(0.0);
    Restrict(grid, rhs, coarse_residual);

    // recursive call
    // std::cout << "recursive" << std::endl;
    Grid* coarse_grid = new Grid(coarse_geom, coarse_offset, coarse_geom->comm());
    coarse_grid->Initialize(0.0);
    MG_Solver* coarse = new MG_Solver(coarse_geom);
    coarse->Iteration(coarse_grid, coarse_residual);

    // interpolate and added
    // std::cout << "prolong" << std::endl;
    // std:: cout << "before Prolong:" << collectResidual(grid, rhs) << std::endl;
    ProlongAndAdd(grid, coarse_grid);
    // std:: cout << "after Prolong:" << collectResidual(grid, rhs) << std::endl;
    // std::cout << "smoothing2" << std::endl;
    res = Smooth(grid, rhs);
    // std::cout << "after smooth2:" << res << std::endl;
  }
}


Grid* MG_Solver::Restrict(const Grid *fine_grid, const Grid* fine_rhs, Grid* coarse_residual) const {

  // bestimmen der Residuen
  InteriorIterator it(coarse_residual->getGeometry());
  // std::cout << "coarse_geom:"  << coarse_residual->getGeometry()->Size()[0] << "x" << coarse_residual->getGeometry()->Size()[1] << std::endl;
  // std::cout << "fine_geom:"  << _geom->Size()[0] << "x" << _geom->Size()[1] << std::endl;
  for (it.First(); it.Valid(); it.Next()) {
    Iterator fine_it(_geom, 2*it.Pos()[0],  2*it.Pos()[1]);
    // coarse_residual->Cell(it) = -0.25* ( localRes(fine_it, fine_grid, fine_rhs)
    //                                     + localRes(fine_it.Right(), fine_grid, fine_rhs)
    //                                     + localRes(fine_it.Top(), fine_grid, fine_rhs)
    //                                     + localRes(fine_it.Right().Top(), fine_grid, fine_rhs)
    //                                     );
    coarse_residual->Cell(it) = 0.25* ( localRes(fine_it, fine_grid, fine_rhs)
                                        + localRes(fine_it.Left(), fine_grid, fine_rhs)
                                        + localRes(fine_it.Down(), fine_grid, fine_rhs)
                                        + localRes(fine_it.Down().Left(), fine_grid, fine_rhs)
                                         );
  }
  return coarse_residual;
}

real_t MG_Solver::Smooth(Grid* grid, const Grid *rhs) const {
  real_t res = 0.0;
  for(size_t i = 0; i < 2 ; i++ ){
    res = _smoother->Cycle(grid, rhs);
  }
  return res;
}

void  MG_Solver::ProlongAndAdd(Grid* fine_grid, const Grid* coarse_grid) const {
  // std::cout << "Prolong fine_geom:"  << _geom->Size()[0] << "x" << _geom->Size()[1] << std::endl;
  InteriorIterator it(coarse_grid->getGeometry());
  for (it.First(); it.Valid(); it.Next()) {
    Iterator fine_it(_geom, 2*it.Pos()[0],  2*it.Pos()[1]);
    // std::cout << it.Pos()[0] << "x" << it.Pos()[1] <<std::endl;
    assert(!std::isnan(coarse_grid->Cell(it)));

    // eventuell interpolieren
    // fine_grid->Cell(fine_it) += coarse_grid->Cell(it);
    // fine_grid->Cell(fine_it.Top()) += coarse_grid->Cell(it);
    // fine_grid->Cell(fine_it.Top().Right()) += coarse_grid->Cell(it);
    // fine_grid->Cell(fine_it.Right()) += coarse_grid->Cell(it);

    fine_grid->Cell(fine_it) += coarse_grid->Cell(it);
    fine_grid->Cell(fine_it.Left()) += coarse_grid->Cell(it);
    fine_grid->Cell(fine_it.Down().Left()) += coarse_grid->Cell(it);
    fine_grid->Cell(fine_it.Down()) += coarse_grid->Cell(it);
  }
}
