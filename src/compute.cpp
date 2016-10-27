#include "compute.hpp"

Compute::Compute (const Geometry *geom, const Parameter *param) {
  _geom = geom;
  _param = param;
  _t = 0.0;

  //Erzeugen der Fehlenden Objekte und zuweisen von Größen
}

Compute::~Compute () {
  // Loeschen der erzeugten Gitter etc
}
