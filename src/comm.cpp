#include "comm.hpp"
#include "grid.hpp"
#include "geometry.hpp"

#include <mpi.h>
#include <assert.h>
#include "stdlib.h"
#include "iterator.hpp"

Communicator::Communicator
(
   int* argc,
   char*** argv
)
{
   //todo evenodd
   MPI_Init( argc, argv );

   MPI_Comm_rank( MPI_COMM_WORLD, &_rank );

   MPI_Comm_size( MPI_COMM_WORLD, &_size );

   assert( DIM == 2);
   int dims[DIM] = {0, 0};
   MPI_Dims_create( _size, DIM, dims );

   for( index_t i = 0; i < DIM; i++)
   {
      _tdim[i] = dims[i];
   }

   int periodic = 0;
   int reorder = 0;
   MPI_Comm mpi_communicator;
   MPI_Cart_create( MPI_COMM_WORLD, (int)DIM, dims, &periodic, reorder, &mpi_communicator );

   int pos[DIM];
   MPI_Cart_coords( mpi_communicator, _rank, (int)DIM, pos );

   for( int i = 0; i < DIM; i++)
   _tidx[i] = pos[i];
}



Communicator::~Communicator
(
   void
)
{
   MPI_Finalize();
}



const multi_index_t&
Communicator::ThreadIdx
(
   void
) const
{
   return _tidx;
}



const multi_index_t&
Communicator::ThreadDim
(
   void
) const
{
   return _tdim;
}



const int& Communicator::ThreadNum() const
{
   return _rank;
}



const int& Communicator::ThreadCnt() const
{
   return _size;
}



const bool& Communicator::EvenOdd() const
{
   return _evenodd;
}



real_t Communicator::geatherSum(const real_t& val) const
{
   real_t r_value;
   MPI_Allreduce( &val, &r_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return r_value;
}



real_t Communicator::geatherMin( const real_t& val ) const
{
   real_t r_value;
   MPI_Allreduce( &val, &r_value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   return r_value;
}



real_t Communicator::geatherMax( const real_t& val ) const
{
   real_t r_value;
   MPI_Allreduce( &val, &r_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   return r_value;
}



void  Communicator::copyBoundary( Grid* grid ) const
{
   copyBottomBoundary( grid ); // dont copy if bottom, is impemented in method
   copyRightBoundary( grid ); // dont copy if right, is impemented in method
   copyTopBoundary( grid ); // dont copy if Top, is impemented in method
   copyLeftBoundary( grid ); // dont copy if Left, is impemented in method
}



const bool Communicator::isLeft() const
{
   const bool r_isLeft = _tidx[1] == 0 && _tdim[1] != 1;
   return r_isLeft;
}



const bool Communicator::isRight() const
{
  return _tidx[0] == _tdim[1] - 1 && _tdim[1] != 1;
}



const bool Communicator::isTop() const
{
   return _tidx[1] == _tdim[0] - 1 && _tdim[0] != 1;
}



const bool Communicator::isBottom() const
{
   return _tidx[1] == 0 && _tdim[0] != 1;
}



bool Communicator::copyLeftBoundary (Grid* grid) const {
  // tue nichts wenn am linken Rand
  if(this->isLeft())
    return false;

  const index_t height = grid->Size()[1];
  real_t* buffer = (real_t*) malloc(height * sizeof(real_t));
  const Geometry* geom = grid->getGeometry();

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
  MPI_Status stat;
  MPI_Sendrecv_replace( buffer, height, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, &stat );
  // zurück kopieren
  for(it.First(); it.Valid(); it.Next()){
    grid->Cell(it) = buffer[it];
    ++i;
  }
  return true;
}


bool Communicator::copyRightBoundary(Grid* grid) const
{

   if(this->isRight())
     return false;

   const index_t height = grid->Size()[1];
   real_t* buffer = (real_t*) malloc(height * sizeof(real_t));
   const Geometry* geom = grid->getGeometry();

   // linker Rand
   BoundaryIterator it(geom);
   it.SetBoundary(2);

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
   MPI_Status stat;
   MPI_Sendrecv_replace( buffer, height, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, &stat );
   // zurück kopieren
   for(it.First(); it.Valid(); it.Next())
   {
     grid->Cell(it) = buffer[it];
     ++i;
   }
   return true;
 }



bool Communicator::copyTopBoundary(Grid* grid) const
{
   if(this->isTop())
     return false;

   const index_t width = grid->Size()[0];
   real_t* buffer = (real_t*) malloc(width * sizeof(real_t));
   const Geometry* geom = grid->getGeometry();


   BoundaryIterator it(geom);
   it.SetBoundary(3);

   // Auslesen der Werte
   index_t i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     buffer[it] = grid->Cell(it.Right());
     ++i;
   }
   const int tag = 0;
   //Adresse einfügen
   const int dest = 0;
   // senden
   MPI_Status stat;
   MPI_Sendrecv_replace( buffer, width, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, &stat );
   // zurück kopieren
   for(it.First(); it.Valid(); it.Next())
   {
     grid->Cell(it) = buffer[it];
     ++i;
   }
   return true;
 }



bool Communicator::copyBottomBoundary(Grid* grid) const
{
   if(this->isTop())
     return false;

   const index_t width = grid->Size()[0];
   real_t* buffer = (real_t*) malloc(width * sizeof(real_t));
   const Geometry* geom = grid->getGeometry();


   BoundaryIterator it(geom);
   it.SetBoundary(1);

   // Auslesen der Werte
   index_t i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     buffer[it] = grid->Cell(it.Right());
     ++i;
   }
   const int tag = 0;
   //Adresse einfügen
   const int dest = 0;
   // senden
   MPI_Status stat;
   MPI_Sendrecv_replace( buffer, width, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, &stat );
   // zurück kopieren
   for(it.First(); it.Valid(); it.Next())
   {
     grid->Cell(it) = buffer[it];
     ++i;
   }
   return true;
 }
