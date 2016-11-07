#include "iterator.hpp"
#include "geometry.hpp"

Iterator::Iterator( const Geometry *geom ){
  _geom = geom;
  _value = 0;
  const multi_index_t& geom_size = geom->Size();
  _valid = (geom_size[0] > 0) && (geom_size[1] > 0);
}



Iterator::Iterator(const Geometry *geom, const index_t &value)
{
  _geom = geom;
  _value = value;
  const multi_index_t& geom_size = geom->Size();
  _valid = ( geom_size[0] * geom_size[1] > value );
}


Iterator::~Iterator( )
{

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
  return _value < geom_size[0]*geom_size[1];
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
    return Iterator( _geom, - 1 ); // invalid Iterator
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
    return Iterator(_geom, -1 ); // invalid Iterator
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
    return Iterator(_geom, -1 ); // invalid Iterator!
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
    return Iterator(_geom, -1); // invalid Iterator
  }
  else
  {
    return Iterator(_geom, _value - geom_size[0]);
  }
}





/// Constructs a new BoundaryIterator
//BoundaryIterator::BoundaryIterator(const Geometry *geom)
//{
//
//
//}

BoundaryIterator::~BoundaryIterator( )
{

}

/// Sets the boundary to iterate
void BoundaryIterator::SetBoundary(const index_t &boundary)
{


}

/// Sets the iterator to the first element
void BoundaryIterator::First()
{


}
/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next()
{


}
