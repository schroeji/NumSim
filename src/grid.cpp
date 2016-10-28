#include <algorithm>
#include "stdlib.h"
#include "math.h"

#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

Grid::Grid(const Geometry *geom){
  multi_index_t geom_size = geom->Size();
  _geom = geom;
  _data = (real_t*) malloc(geom_size[0] * geom_size[1] * sizeof(real_t));
  _offset = {0.0, 0.0};
}

Grid::Grid(const Geometry *geom, const multi_real_t &offset) {
  multi_index_t geom_size = geom->Size();
  _geom = geom;
  _data = (real_t*) malloc(geom_size[0] * geom_size[1] * sizeof(real_t));
  _offset = offset;
}

void Grid::Initialize(const real_t &value) {
  Iterator it(_geom);
  for (it.First(); it.Valid(); it.Next()) {
    _data[it] = value;
  }
}

real_t& Grid::Cell(const Iterator &it) {
  return _data[it];
}

const real_t& Grid::Cell(const Iterator &it) const {
  return _data[it];
}

real_t Grid::Max() const {
  Iterator it(_geom);
  real_t max_val = _data[0];
  for (it.First(); it.Valid(); it.Next()) {
    max_val = std::max(max_val, _data[it]);
  }
  return max_val;
}

real_t Grid::Min() const {
  Iterator it(_geom);
  real_t min_val = _data[0];
  for (it.First(); it.Valid(); it.Next()) {
    min_val = std::min(min_val, _data[it]);
  }
  return min_val;
}

real_t Grid::AbsMax() const {
  Iterator it(_geom);
  real_t max_val = fabs(_data[0]);
  for (it.First(); it.Valid(); it.Next()) {
    max_val = std::max(max_val, fabs(_data[it]));
  }
  return max_val;
}

real_t* Grid::Data(){
  return _data;
}
Grid::~Grid(){
  free (_data);
}
