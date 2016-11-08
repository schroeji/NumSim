#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "geometry.hpp"

Geometry::Geometry() {
  _size = {128, 128};
  _length = {1.0, 1.0};
  _h = {_length[0] / _size[0], _length[1] / _size[1]};
  _velocity = {0.0, 0.0};
  _pressure = 0.1;
}

void Geometry::Load(const char *file){
  FILE* handle = fopen(file,"r");
  double inval[2] = {0.0, 0.0};
  char name[20];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s =", name)) continue;
    if (strcmp(name,"size") == 0) {
      if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _size[0] = inval[0];
        _size[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name,"length") == 0) {
      if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _length[0] = inval[0];
        _length[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name,"velocity") == 0) {
      if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _velocity[0] = inval[0];
        _velocity[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name,"pressure") == 0) {
      if (fscanf(handle," %lf\n",&inval[0]))
        _pressure = inval[0];
      continue;
    }
  }
  fclose(handle);

}

const multi_index_t& Geometry::Size() const {
  return _size;
}

const multi_real_t& Geometry::Length() const {
  return _length;
}

const multi_real_t& Geometry::Mesh() const {
  return _h;
}
