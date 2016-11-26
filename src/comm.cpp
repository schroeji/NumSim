#include "comm.hpp"
#include <mpi.h>
#include <assert.h>


Communicator::Communicator
(
   int* argc,
   char*** argv
)
{
   MPI_Init( argc, argv );

   MPI_Comm_rank( MPI_COMM_WORLD, &_rank );

   MPI_Comm_size( MPI_COMM_WORLD, &_size );

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
   assert( false && "ToDo: not implemented yet");
}



const bool Communicator::isLeft() const
{
   assert( false && "ToDo: not implemented yet");
}



const bool Communicator::isRight () const
{
   assert( false && "ToDo: not implemented yet");
}



const bool Communicator::isTop() const
{
   assert( false && "ToDo: not implemented yet");
}



const bool Communicator::isBottom() const
{
   assert( false && "ToDo: not implemented yet");
}




bool Communicator::copyLeftBoundary(Grid* grid) const
{
   assert( false && "ToDo: not implemented yet");
}



bool Communicator::copyRightBoundary(Grid* grid) const
{
   assert( false && "ToDo: not implemented yet");
}



bool Communicator::copyTopBoundary(Grid* grid) const
{
   assert( false && "ToDo: not implemented yet");
}



bool Communicator::copyBottomBoundary(Grid* grid) const
{
   assert( false && "ToDo: not implemented yet");
}



