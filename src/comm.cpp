#include "stdlib.h"
#inlcude "mpi/mpi.h"

#include "comm.hpp"
#include "grid.hpp"
#include "geometry.hpp"

bool Communicator::copyLeftBoundary	(Grid* grid) const {
  // tue nichts wenn am linken Rand
  if(this->isLeft())
    return false;

  const index_t height = grid->Size()[1];
  real_t* buffer = (real_t*) malloc(height * sizeof(real_t));
  Geometry* geom = grid->getGeometry();

  // linker Rand
  BoundaryIterator it(geom);
  it.SetBoundary(4);

  // Auslesen der Werte
  index_t i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    buffer[it] = grid->Cell(it.Right());
    ++i;
  }
  const int tag = 0;
  //Adresse einfügen
  const int dest = 0;
  // senden
  MPI_Sendrecv_replace(buffer, height, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, &stat);
  // zurück kopieren
  for(it.First(); it.Valid(); it.Next()){
    grid->Cell(it) = buffer[it];
    ++i;
  }
}

bool Communicator::copyRightBoundary	(Grid* grid) const {

}

bool Communicator::copyTopBoundary	(Grid* grid) const {
}

bool Communicator::copyBottomBoundary	(Grid* grid) const {

}
