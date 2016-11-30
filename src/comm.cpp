#include "comm.hpp"
#include "geometry.hpp"

#include <assert.h>
#include "stdlib.h"
#include "iterator.hpp"
#include "iostream"

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

   MPI_Cart_create( MPI_COMM_WORLD, (int)DIM, dims, &periodic, reorder, &_mpi_communicator );

   int pos[DIM];
   MPI_Cart_coords( _mpi_communicator, _rank, (int)DIM, pos );
   for( int i = 0; i < DIM; i++) {
       _tidx[i] = pos[i];
   }
   _evenodd = (_tidx[0] % 2 == 1) ^ (_tidx[1] % 2 == 1);
   std::cout << "Rank: " << _rank << " Pos: " << pos[0] << ";" << pos[1] << " even: " << _evenodd << std::endl;
   std::cout << "tdim: " << _tdim[0] << ";" << _tdim[1] << std::endl;
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
   MPI_Allreduce( &val, &r_value, 1, MPI_DOUBLE, MPI_SUM, _mpi_communicator);
   return r_value;
}



real_t Communicator::geatherMin( const real_t& val ) const
{
   real_t r_value;
   MPI_Allreduce( &val, &r_value, 1, MPI_DOUBLE, MPI_MIN, _mpi_communicator);
   return r_value;
}



real_t Communicator::geatherMax( const real_t& val ) const
{
   real_t r_value;
   MPI_Allreduce( &val, &r_value, 1, MPI_DOUBLE, MPI_MAX, _mpi_communicator);
   return r_value;
}



void  Communicator::copyBoundary( Grid* grid ) const
{
  std::cout << "entering copyboundary Rank:" << _rank << std::endl;
  if(_evenodd) {
    copyBottomBoundary( grid ); // dont copy if bottom, is impemented in method
    copyTopBoundary( grid ); // dont copy if Top, is impemented in method
    copyRightBoundary( grid ); // dont copy if right, is impemented in method
    copyLeftBoundary( grid ); // dont copy if Left, is impemented in method
  }
  else {
    copyTopBoundary( grid ); // dont copy if Top, is impemented in method
    copyBottomBoundary( grid ); // dont copy if bottom, is impemented in method
    copyLeftBoundary( grid ); // dont copy if Left, is impemented in method
    copyRightBoundary( grid ); // dont copy if right, is impemented in method
  }
  std::cout << "exiting copyboundary Rank:" << _rank << std::endl;
}



const bool Communicator::isLeft() const
{
   // const bool r_isLeft = _tidx[1] == 0 && _tdim[1] != 1;
   // return r_isLeft;
  return _tidx[0] == 0;
}



const bool Communicator::isRight() const
{
  // return _tidx[0] == _tdim[1] - 1 && _tdim[1] != 1;
  return _tidx[0] == _tdim[0] - 1;
}



const bool Communicator::isTop() const
{
   // return _tidx[1] == _tdim[0] - 1 && _tdim[0] != 1;
  return _tidx[1] == _tdim[1] - 1;
}



const bool Communicator::isBottom() const
{
   // return _tidx[1] == 0 && _tdim[0] != 1;
  return _tidx[1] == 0;
}



bool Communicator::copyLeftBoundary (Grid* grid) const {
  // tue nichts wenn am linken Rand
  if(this->isLeft())
    return false;

  const index_t height = grid->Size()[1] + 2;
  real_t* buffer = (real_t*) malloc(height * sizeof(real_t));
  const Geometry* geom = grid->getGeometry();

  // linker Rand
  BoundaryIterator it(geom);
  it.SetBoundary(4);

  // Auslesen der Werte
  index_t i = 0;
  for(it.First(); it.Valid(); it.Next()) {
    buffer[i] = grid->Cell(it.Right());
    ++i;
  }
  // std::cout << "read into buffer rank:" << _rank << std::endl;
  const int tag = 0;
  //Adresse einfügen
  const int dest = _rank - _tdim[1];
  // senden
  MPI_Status stat;
  std::cout << "send_rcv from: " << _rank << " to: " << dest << std::endl;
  MPI_Sendrecv_replace( buffer, height, MPI_DOUBLE, dest, tag, dest, tag, _mpi_communicator, &stat );
  std::cout << "send_rcv from: " << _rank << " to: " << dest <<  " done !" << std::endl;
  // zurück kopieren
  i = 0;
  for(it.First(); it.Valid(); it.Next()){
    grid->Cell(it) = buffer[i];
    ++i;
  }
  return true;
}


bool Communicator::copyRightBoundary(Grid* grid) const
{

   if(this->isRight())
     return false;

   const index_t height = grid->Size()[1] + 2;
   printf("height: %i", height);
   real_t* buffer = (real_t*) malloc(height * sizeof(real_t));
   const Geometry* geom = grid->getGeometry();

   // linker Rand
   BoundaryIterator it(geom);
   it.SetBoundary(2);

   // Auslesen der Werte
   index_t i = 0;
   for(it.First(); it.Valid(); it.Next()) {
     buffer[i] = grid->Cell(it.Left().Left());
     ++i;
   }
   // std::cout << "read into buffer rank:" << _rank << std::endl;
   const int tag = 0;
   //Adresse einfügen
   const int dest = _rank + _tdim[1];
   // senden
   MPI_Status stat;
   std::cout << "send_rcv from: " << _rank << " to: " << dest << std::endl;
   MPI_Sendrecv_replace( buffer, height, MPI_DOUBLE, dest, tag, dest, tag, _mpi_communicator, &stat );
   std::cout << "send_rcv from: " << _rank << " to: " << dest <<  " done !" << std::endl;
   // zurück kopieren
   i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     grid->Cell(it) = buffer[i];
     ++i;
   }
   return true;
 }


bool Communicator::copyTopBoundary(Grid* grid) const
{
   if(this->isTop())
     return false;
   const index_t width = grid->Size()[0] + 2;
   real_t* buffer = (real_t*) malloc(width * sizeof(real_t));
   const Geometry* geom = grid->getGeometry();


   BoundaryIterator it(geom);
   it.SetBoundary(3);

   // Auslesen der Werte
   index_t i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     buffer[i] = grid->Cell(it.Down());
     ++i;
   }
   // std::cout << "read into buffer rank:" << _rank << std::endl;
   const int tag = 0;
   //Adresse einfügen
   const int dest = _rank + 1;
   // senden
   MPI_Status stat;
   std::cout << "send_rcv from: " << _rank << " to: " << dest << std::endl;
   MPI_Sendrecv_replace( buffer, width, MPI_DOUBLE, dest, tag, dest, tag, _mpi_communicator, &stat);
   std::cout << "send_rcv from: " << _rank << " to: " << dest <<  " done !" << std::endl;
   // zurück kopieren
   i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     grid->Cell(it) = buffer[i];
     ++i;
   }

   return true;
 }



bool Communicator::copyBottomBoundary(Grid* grid) const
{
   if(this->isBottom())
     return false;

   const index_t width = grid->Size()[0] + 2;
   real_t* buffer = (real_t*) malloc(width * sizeof(real_t));
   const Geometry* geom = grid->getGeometry();


   BoundaryIterator it(geom);
   it.SetBoundary(1);

   // Auslesen der Werte
   index_t i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     buffer[i] = grid->Cell(it.Top());
     ++i;
   }
   // std::cout << "read into buffer rank:" << _rank << std::endl;
   const int tag = 0;
   //Adresse einfügen
   // const int dest = 0;
   const int dest = _rank - 1;
   // senden
   MPI_Status stat;
   std::cout << "send_rcv from: " << _rank << " to: " << dest << std::endl;
   MPI_Sendrecv_replace( buffer, width, MPI_DOUBLE, dest, tag, dest, tag, _mpi_communicator, &stat );
   std::cout << "send_rcv from: " << _rank << " to: " << dest <<  " done !" << std::endl;
   // zurück kopieren
   i = 0;
   for(it.First(); it.Valid(); it.Next())
   {
     grid->Cell(it) = buffer[i];
     ++i;
   }
   return true;
 }
