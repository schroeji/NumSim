#include "parameter.hpp"
#include <cstdio>
#include <cstring>
Parameter::Parameter () {
  _re = 1000;
  _omega = 1.7;
  _alpha = 0.9;
  _dt = 0.005;
  _tend = 50;
  _itermax = 1000;
  _eps = 0.001;
  _tau = 0.5;
}

void Parameter::Load(const char* file){
  FILE* handle = fopen(file, "r");
  double val;
  char name[20];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s = %lf\n", name, &val))
      continue;
    if (strcmp(name, "re") == 0) _re = val;
    else if (strcmp(name, "omega") == 0) _omega = val;
    else if (strcmp(name, "alpha") == 0) _alpha = val;
    else if (strcmp(name, "dt") == 0) _dt = val;
    else if (strcmp(name, "tend") == 0) _tend = val;
    else if (strcmp(name, "itermax") == 0) _itermax = val;
    else if (strcmp(name, "eps") == 0) _eps = val;
    else if (strcmp(name, "tau") == 0) _tau = val;
    else printf("unknown parameter %s\n", name);
  }
  fclose(handle);
}

const real_t& Parameter::Re() const {
  return _re;
}

const real_t& Parameter::Omega() const {
  return _omega;
}

const real_t& Parameter::Alpha() const{
  return _alpha;
}

const real_t& Parameter::Dt() const {
  return _dt;
}

const real_t& Parameter::Tend() const {
  return _tend;
}

const index_t& Parameter::IterMax() const {
  return _itermax;
}

const real_t& Parameter::Eps() const {
  return _eps;
}

const real_t& Parameter::Tau() const {
  return _tau;
}
