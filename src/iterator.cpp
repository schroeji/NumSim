#include "iterator.hpp"
#include "geometry.hpp"

Iterator::Iterator(const Geometry *geom){
  _geom = geom;
  _value = 0;
  const multi_index_t& geom_size = geom->Size();
  _valid = (geom_size[0] > 0) && (geom_size[1] > 0);
}

Iterator::Iterator(const Geometry *geom, const index_t &value){
  _geom = geom;
  _value = value;
  const multi_index_t& geom_size = geom->Size();
  _valid = (geom_size[0] > value) && (geom_size[1] > value);
}

const index_t& Iterator::Value() const{
  return _value;
}

Iterator::operator const index_t& () const {
  return _value;
}

multi_index_t Iterator::Pos() const{
  index_t row_len = _geom->Size()[0];
  index_t ind1 = _value % row_len;
  index_t ind2 = _value / row_len;
  multi_index_t ret = {ind1, ind2};
  return ret;
}

void Iterator::First() {
  _value = 0;
}

void Iterator::Next() {
  ++_value;
  const multi_index_t& geom_size = _geom->Size();
  _valid = _value < geom_size[0]*geom_size[1];
}

bool Iterator::Valid() const {
  return _valid;
}

Iterator Iterator::Left() const {
  multi_index_t pos = Pos();
  if(pos[0] == 0){
    return Iterator(_geom, _value);
  }
  else{
    return Iterator(_geom, _value - 1);
  }
}

Iterator Iterator::Right() const {
  const multi_index_t& geom_size = _geom->Size();
  multi_index_t pos = Pos();
  if(pos[0] == geom_size[1] - 1){
    return Iterator(_geom, _value);
  }
  else{
    return Iterator(_geom, _value + 1);
  }
}

Iterator Iterator::Top() const {
  multi_index_t pos = Pos();
  if(pos[1] == 0){
    return Iterator(_geom, _value);
  }
  else{
    return Iterator(_geom, _value - _geom->Size()[0]);
  }
}

Iterator Iterator::Down() const {
  const multi_index_t& geom_size = _geom->Size();
  multi_index_t pos = Pos();
  if(pos[1] == geom_size[1]){
    return Iterator(_geom, _value);
  }
  else{
    return Iterator(_geom, _value - geom_size[0]);
  }
}

void InteriorIterator::First() {
  // Element (1,1)
  _value = _geom->Size()[0] + 1;
}

void InteriorIterator::Next() {
  const multi_index_t& geom_size = _geom->Size();
  ++_value;
  // Wenn auf dem rechten Rand nochmal 2 Schritte
  if(_value % geom_size[0] == geom_size[0] - 1){
    ++_value;
    ++_value;
  }
  _valid = _value < geom_size[0]*geom_size[1] - 1 - geom_size[0];
}
