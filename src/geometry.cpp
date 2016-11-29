#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <assert.h>

#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"

Geometry::Geometry() {
  _size = {128, 128};
  _length = {1.0, 1.0};
  _h = {_length[0] / _size[0], _length[1] / _size[1]};
  _velocity = {1.0, 0.0};
  _pressure = 0.0;
}

Geometry::Geometry(const Communicator* comm) {
  index_t x_dim = comm->ThreadDim()[0];
  index_t y_dim = comm->ThreadDim()[1];
  assert(_size[0] % x_dim == 0 && _size[1] % y_dim  == 0);

  _size = {128, 128};
  _bsize = {_size[0]/x_dim, _size[1]/y_dim};

  _length = {1.0, 1.0};
  _blength = {_length[0]/x_dim, _length[1]/y_dim};

  _h = {_blength[0] / _bsize[0], _blength[1] / _bsize[1]};
  _velocity = {1.0, 0.0};
  _pressure = 0.0;
  _comm = comm;
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
      } continue;}
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

  index_t x_dim = _comm->ThreadDim()[0];
  index_t y_dim = _comm->ThreadDim()[1];
  assert(_size[0] % x_dim == 0 && _size[1] % y_dim  == 0);

  _blength = {_length[0]/x_dim, _length[1]/y_dim};
  _bsize = {_size[0]/x_dim, _size[1]/y_dim};

  fclose(handle);
}

const multi_index_t& Geometry::Size() const {
  return _bsize;
}

const multi_index_t& Geometry::TotalSize() const {
  return _size;
}

const multi_real_t& Geometry::Length() const {
  return _blength;
}

const multi_real_t& Geometry::TotalLength() const {
  return _length;
}

const multi_real_t& Geometry::Mesh() const {
  return _h;
}


// Updates the velocity field u
void
Geometry::Update_U
(
   Grid *u
) const
{
  _comm->copyBoundary(u);
   BoundaryIterator it( this );
   if( u->isBottom() )
   {
      it.SetBoundary(1);
      for( it.First(); it.Valid(); it.Next() )
      {
        u->Cell( it ) = 0.0 - u->Cell(it.Top()) ;
      }
   }
   if(u->isRight())
   {
      it.SetBoundary(2);
      for( it.First(); it.Valid(); it.Next() )
      {
        u->Cell( it ) = 0.0;
        u->Cell( it.Left() ) = 0.0;
      }
   }

   if( u->isTop() )
   {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next()) {
        u->Cell(it) = 2*_velocity[0] - u->Cell(it.Down());
      }
   }

   if( u->isLeft() )
   {
      it.SetBoundary(4);
      for( it.First(); it.Valid(); it.Next() )
        {
          u->Cell( it ) = 0.0;
        }
   }
}


// Updates the velocity field v
void
Geometry::Update_V
(
   Grid *v
) const
{
  _comm->copyBoundary(v);
   BoundaryIterator it( this );
   if( v->isBottom() )
   {
      it.SetBoundary(1);
      for( it.First(); it.Valid(); it.Next() )
      {
         v->Cell( it ) = 0.0;
      }
   }

   if( v->isRight())
   {
      it.SetBoundary(2);
      for( it.First(); it.Valid(); it.Next() )
      {
         v->Cell( it ) = 0.0 - v->Cell( it.Left() );
      }
   }

   if( v->isTop() )
   {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next()) {
         v->Cell(it) = 0.0;
         v->Cell(it.Down()) = 0.0;
      }
   }

   if( v->isLeft() )
   {
      it.SetBoundary(4);
      for( it.First(); it.Valid(); it.Next() )
      {
        v->Cell( it ) = 0.0 - v->Cell(it.Right() );
      }
   }
}


/// Updates the pressure field p
void
Geometry::Update_P
(
   Grid *p
) const
{
  _comm->copyBoundary(p);
   BoundaryIterator it( this );
   if( p->isBottom() )
   {
      it.SetBoundary(1);
      for( it.First(); it.Valid(); it.Next() )
      {
        p->Cell( it ) = p->Cell(it.Top());
      }
   }

   if( p->isRight() )
   {
      it.SetBoundary(2);
      for( it.First(); it.Valid(); it.Next() )
      {
        p->Cell( it ) = p->Cell(it.Left());
      }
   }

   if( p->isTop() )
   {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next()) {
        p->Cell(it) = p->Cell(it.Down());
      }
   }

   if( p->isLeft() )
   {
      it.SetBoundary(4);
      for( it.First(); it.Valid(); it.Next() )
      {
        p->Cell( it ) = p->Cell(it.Right());
      }
   }
}
