#include "iterator.hpp"
#include "geometry.hpp"

Iterator::Iterator( const Geometry *geom ){
  _geom = geom;
  _value = 0;
  const multi_index_t& geom_size = geom->Size();
  _valid = (geom_size[0] > 0) && (geom_size[1] > 0);
}



Iterator::Iterator( const Geometry *geom, const index_t &value )
{
  _geom = geom;
  _value = value;
  _valid = Valid();
}



const index_t&
Iterator::Value
(
   void
) const
{
  return _value;
}



Iterator::operator const index_t& () const
{
  return _value;
}



multi_index_t Iterator::Pos() const
{
  index_t row_len = _geom->Size()[0];
  index_t ind1 = _value % row_len;
  index_t ind2 = _value / row_len;
  multi_index_t ret = {ind1, ind2};
  return ret;
}



void Iterator::First()
{
  _value = 0;
}



void Iterator::Next()
{
  ++_value;
  Valid();
}



bool
Iterator::Valid
(
   void
) const
{
  const multi_index_t& geom_size = _geom->Size();
  return ( _value + 1 ) < geom_size[0]*geom_size[1];
}



bool
Iterator::isInteriorIterator
(
   void
) const
{
   bool r_isInterior  =    _value > 0                         // not the bottom boundary
                        && ( ( _value + 1 ) % _geom->Size()[1] != 0 ) // not on right boundary
                        && ( ( _value + 1 ) % _geom->Size()[1] != 1 ) // not on left boundary
                        && ( ( _value + 1 ) < ( _geom->Size()[0] * ( _geom->Size()[1] - 1 ) ) ); // not the upper boundary
   return r_isInterior;
}



Iterator
Iterator::Left
(
   void
) const
{
  multi_index_t pos = Pos();
  if( pos[1] == 0 )
  {
    return Iterator( _geom, _geom->Size()[0]*_geom->Size()[1] ); // invalid Iterator
  }
  else
  {
    return Iterator(_geom, _value - 1);
  }
}



Iterator
Iterator::Right
(
   void
) const
{
  const multi_index_t& geom_size = _geom->Size();
  multi_index_t pos = Pos();
  if( pos[0] == geom_size[1] - 1)
  {
    return Iterator(_geom, _geom->Size()[0]*_geom->Size()[1] ); // invalid Iterator
  }
  else
  {
    return Iterator(_geom, _value + 1);
  }
}



Iterator
Iterator::Top
(
   void
) const
{
  multi_index_t pos = Pos();
  if( pos[1] == 0 )
  {
    return Iterator(_geom, _geom->Size()[0]*_geom->Size()[1] ); // invalid Iterator!
  }
  else
  {
    return Iterator(_geom, _value - _geom->Size()[0]);
  }
}



Iterator
Iterator::Down
(
   void
) const
{
  const multi_index_t& geom_size = _geom->Size();
  multi_index_t pos = Pos();
  if(pos[1] == geom_size[1])
  {
    return Iterator(_geom, _geom->Size()[0]*_geom->Size()[1]); // invalid Iterator
  }
  else
  {
    return Iterator(_geom, _value - geom_size[0]);
  }
}




/// Constructs a new BoundaryIterator
BoundaryIterator::BoundaryIterator(const Geometry *geom) : Iterator( geom )
{
   _value = _geom->Size()[0]*_geom->Size()[1];
   _boundary = -1;
   _valid = false; // create invalid BoundaryIterator;
}



/// Sets the boundary to iterate
void
BoundaryIterator::SetBoundary
(
   const index_t &boundary
)
{
   // _boundary == 1 under boundary
   // _boundary == 2 right boundary
   // _boundary == 3 upper boundary
   // _boundary == 4 left boundary
   assert( ( 1 <= boundary ) && ( 4 >= boundary ) );

   _boundary = boundary;
}



bool
BoundaryIterator::Valid
(
   void
) const
{
   bool r_isValid = Iterator::Valid() && !Iterator::isInteriorIterator();
   return r_isValid;
}



/// Sets the iterator to the first element
void BoundaryIterator::First()
{
   if( _boundary == 1 )
   {
      _value = 0;
   }
   else if( _boundary == 2 )
   {
      _value = _geom->Size()[0] - 1;
   }
   else if( _boundary == 3 )
   {
      _value = _geom->Size()[0]*( _geom->Size()[1] - 1 );
   }
   else if( _boundary == 4 )
   {
      _value = 0;
   }
}



/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next()
{
   if( _boundary == 1 ) // under boundary
   {
      _value++;
   }
   else if( _boundary == 2 ) // right boundary
   {
      _value += _geom->Size()[1];
   }
   else if( _boundary == 3  ) // upper boundary
   {
      _value++;
   }
   else if( _boundary == 4 ) // left boundary
   {
      _value += _geom->Size()[0];
   }
}



InteriorIterator::InteriorIterator(const Geometry *geom) : Iterator( geom )
{
   _valid = Valid();
}



bool
InteriorIterator::Valid
(
   void
) const
{
   return Iterator::isInteriorIterator( );
}



void
InteriorIterator::First
(
   void
) {
  // Element (1,1)
  _value = _geom->Size()[0] + 1;
  _valid = true;
}

void
InteriorIterator::Next
(
   void
)
{
  const multi_index_t& geom_size = _geom->Size();
  ++_value;
 // Wenn auf dem rechten Rand nochmal 2 Schritte
  if(_value % geom_size[0] == geom_size[0] - 1)
  {
    _value += 2;
  }
  _valid = Valid();
}