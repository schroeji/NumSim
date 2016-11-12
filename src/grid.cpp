#include <algorithm>
#include "stdlib.h"
#include "math.h"
#include <assert.h>
#include <iostream>

#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

Grid::Grid(const Geometry *geom)
{
  multi_index_t geom_size = geom->Size();
  _geom = geom;
  _data = (real_t*) malloc( (geom_size[0] + 2) * (geom_size[1] + 2) * sizeof(real_t));
  _offset = {0.0, 0.0};
}



Grid::Grid(const Geometry *geom, const multi_real_t &offset)
{
  multi_index_t geom_size = geom->Size();
  _geom = geom;
  _data = (real_t*) malloc(geom_size[0] * geom_size[1] * sizeof(real_t));
  _offset = offset;
}



Grid::~Grid()
{
  free (_data);
}



void Grid::Initialize(const real_t &value)
{
  Iterator it(_geom);
  for ( it.First(); it.Valid(); it.Next() )
  {
    _data[ it.Value( ) ] = value;
  }
}




real_t
Grid::Interpolate
(
   const multi_real_t &pos
) const
{
  std ::cout << pos[0] << ";" << pos[1] << std::endl;
  // std ::cout << _offset[0] << ";" << _offset[1] << std::endl;
  index_t xPosition = static_cast< index_t >( (pos[0] - _offset[0]) / _geom->Mesh()[0] + 1); // round down
  index_t yPosition = static_cast< index_t >( (pos[1] - _offset[1]) / _geom->Mesh()[1] + 1); // round down

  // std ::cout << (pos[0] - _offset[0]) / _geom->Mesh()[0] << ";" << (pos[1] - _offset[1])  / _geom->Mesh()[1]<< std::endl;
  std ::cout << xPosition << ";" << yPosition << std::endl;
   Iterator asseccIterator( _geom, xPosition + ( _geom->Size()[0] ) * yPosition );
   auto weightForX =  fmod( xPosition, _geom->Mesh()[0] ) /_geom->Mesh()[0];
   auto weightForY =  fmod( yPosition, _geom->Mesh()[1] ) / _geom->Mesh()[1];
   // std ::cout << asseccIterator << ";" << asseccIterator.Top() << std::endl;
   real_t r_interpolate =     ( 1 - weightForY )
                            * ( ( 1 - weightForX ) * Cell( asseccIterator )
                              + weightForX * ( asseccIterator.Right() ) )
                          + weightForY
                            * ( ( 1 - weightForX ) * Cell( asseccIterator.Top() )
                              + weightForX * Cell( asseccIterator.Top().Right() ) );

    return r_interpolate;
}



real_t& Grid::Cell(const Iterator &it)
{
  return _data[ it.Value() ];
}



const real_t& Grid::Cell(const Iterator &it) const
{
  // std ::cout << it << std::endl;
   assert( it.Valid() );
   return _data[ it.Value() ];
}



real_t Grid::Max() const {
  Iterator it(_geom);
  real_t max_val = _data[0];
  for (it.First(); it.Valid(); it.Next())
  {
    max_val = std::max(max_val, _data[it.Value()]);
  }
  return max_val;
}



real_t Grid::Min() const
{
  Iterator it(_geom);
  real_t min_val = _data[0];
  for (it.First(); it.Valid(); it.Next())
  {
    min_val = std::min(min_val, _data[it.Value()]);
  }

  return min_val;
}



real_t Grid::AbsMax() const
{
  Iterator it(_geom);
  real_t max_val = fabs(_data[0]);
  for (it.First(); it.Valid(); it.Next())
  {
    max_val = std::max(max_val, fabs(_data[it.Value()]));
  }

  return max_val;
}



real_t* Grid::Data()
{
  return _data;
}



real_t
Grid::dx_l
(
   const Iterator &it
) const
{
   assert( it.Left().Valid() && it.Valid() );
   real_t r_diff = ( Cell( it.Left( ) ) -  Cell( it ) )/_geom->Mesh()[0];
   return r_diff;
}



real_t
Grid::dx_r
(
   const Iterator &it
) const
{
   assert( it.Valid() && it.Right().Valid() );
   real_t r_diff = ( Cell( it ) - Cell( it.Right( ) ) )/_geom->Mesh()[0];
    return r_diff;
}



real_t
Grid::dy_r
(
      const Iterator &it
) const
{
   assert( it.Valid() && it.Top().Valid() );
   real_t r_diff = ( Cell( it.Top() ) - Cell( it ) )/_geom->Mesh()[1];
   return r_diff;
}



real_t
Grid::dy_l
(
      const Iterator &it
) const
{
   assert( it.Valid() && it.Down().Valid() );
   auto valueOfCell = Cell(it );
   auto valueOfCellDown = Cell( it.Down() );
   real_t r_diff = ( Cell( it ) - Cell( it.Down() ) )/_geom->Mesh()[1];
   return r_diff;
}



real_t
Grid::dxx
(
      const Iterator &it
) const
{
   real_t r_ddiff = ( dx_r(it) - dx_l(it) )/_geom->Mesh()[0];
   return r_ddiff;
}



real_t
Grid::dyy
(
      const Iterator &it
) const
{
   real_t r_ddiff = ( dy_r(it) - dy_l(it) )/_geom->Mesh()[1];
   return r_ddiff;
}



/// Computes u*du/dx with the donor cell method
real_t
Grid::DC_udu_x
(
   const Iterator &it,
   const real_t &alpha
) const
{
   // only if grid are u grid.
   real_t firstTerm =  0.25 *( ( Cell( it ) + Cell( it.Right() ) ) * ( Cell( it ) + Cell( it.Right() ) )
                             - ( Cell( it.Left() ) + Cell( it ) ) * ( Cell( it.Left() ) + Cell( it ) ) ) ;

   real_t secondTerm = alpha * 0.25 * ( std::abs( Cell( it ) + Cell( it.Right() ) ) * ( Cell( it ) - Cell( it.Right() ) )
                                         - std::abs( Cell( it.Left() ) + Cell( it ) ) * ( Cell( it.Left() ) - Cell( it ) )  );
   real_t r_DC_udu_x = 1.0/_geom->Mesh()[0] * firstTerm + 1.0/_geom->Mesh()[0]*secondTerm;
   return r_DC_udu_x;
}



/// Computes v*du/dy with the donor cell method
real_t
Grid::DC_vdu_y
(
   const Iterator &it,
   const real_t &alpha,
   const Grid *v
) const
{
   real_t firstTerm =  0.25 *( ( v->Cell( it ) + v->Cell( it.Right() ) )
                               * ( Cell( it ) + Cell( it.Top() ) )
                             - ( v->Cell( it.Down() ) + v->Cell( it.Down().Right() ) )
                               * ( Cell( it.Down() ) + Cell( it ) ) ) ;

   real_t secondTerm =  alpha
                      * 0.25
                      * ( std::abs( v->Cell( it ) + v->Cell( it.Right() ) )
                        * ( Cell( it ) - Cell( it.Top() ) )
                      - std::abs( v->Cell( it.Down() ) + v->Cell( it.Down().Right() ) )
                        * ( Cell( it.Down() ) - Cell( it ) )  );
   real_t r_DC_vdu_y = 1.0/_geom->Mesh()[1] * firstTerm + 1.0/_geom->Mesh()[1] * secondTerm;
   return r_DC_vdu_y;
}



/// Computes u*dv/dx with the donor cell method
real_t
Grid::DC_udv_x
(
   const Iterator &it,
   const real_t &alpha,
   const Grid *u
) const
{
   real_t firstTerm =  0.25 *( ( u->Cell( it ) + u->Cell( it.Right() ) )
                               * ( Cell( it ) + Cell( it.Right() ) )
                             - ( u->Cell( it.Left() ) + u->Cell( it ) )
                               * ( Cell( it.Left() ) + Cell( it ) ) ) ;

   real_t secondTerm =  alpha
                      * 0.25
                      * ( std::abs( u->Cell( it ) + u->Cell( it.Right() ) )
                        * ( Cell( it ) - Cell( it.Right() ) )
                      - std::abs( u->Cell( it.Left() ) + u->Cell( it ) )
                        * ( Cell( it.Left() ) - Cell( it ) )  );
   real_t r_DC_udv_x = 1.0/_geom->Mesh()[0] * firstTerm + 1.0/_geom->Mesh()[0] * secondTerm;
   return r_DC_udv_x;
}



/// Computes v*dv/dy with the donor cell method
real_t
Grid::DC_vdv_y
(
   const Iterator &it,
   const real_t &alpha
) const
{
   real_t firstTerm =  0.25 *( ( Cell( it ) + Cell( it.Down() ) ) * ( Cell( it ) + Cell( it.Down() ) )
                             - ( Cell( it.Top() ) + Cell( it ) ) * ( Cell( it.Top() ) + Cell( it ) ) ) ;

   real_t secondTerm = alpha * 0.25 * ( std::abs( Cell( it ) + Cell( it.Top() ) ) * ( Cell( it ) - Cell( it.Down() ) )
                                         - std::abs( Cell( it.Top() ) + Cell( it ) ) * ( Cell( it.Top() ) - Cell( it ) )  );
   real_t r_DC_udu_y = 1.0/_geom->Mesh()[1] * firstTerm + 1.0/_geom->Mesh()[1]*secondTerm;
   return r_DC_udu_y;
}
