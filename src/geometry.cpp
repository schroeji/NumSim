#include <cstdio>
#include <cstring>
#include <cstdlib>

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


/// Updates the velocity field u
void
Geometry::Update_U
(
   Grid *u
) const
{

   BoundaryIterator it( this );
   it.SetBoundary(1);
   for( it.First(); it.Valid(); it.Next() )
   {
     u->Cell( it ) = 0.0 - u->Cell(it.Top()) ;
   }

   it.SetBoundary(2);
   for( it.First(); it.Valid(); it.Next() )
   {
     u->Cell( it ) = 0.0;
   }

   it.SetBoundary(4);
   for( it.First(); it.Valid(); it.Next() )
   {
     u->Cell( it ) = 0.0;
   }

   // als letztes setzen wegen doppelter eckpunkte
   it.SetBoundary(3);
   for (it.First(); it.Valid(); it.Next()) {
     u->Cell(it) = 2*_velocity[0] - u->Cell(it.Down());
   }
}


/// Updates the velocity field v
void
Geometry::Update_V
(
   Grid *v
) const
{
   BoundaryIterator it( this );


   it.SetBoundary(1);
   for( it.First(); it.Valid(); it.Next() )
   {
      v->Cell( it ) = 0.0;
   }

   it.SetBoundary(2);
   for( it.First(); it.Valid(); it.Next() )
   {
     v->Cell( it ) = 0.0 - v->Cell( it.Left() );
   }
   it.SetBoundary(3);
   for (it.First(); it.Valid(); it.Next()) {
      v->Cell(it) = 0.0;
   }
   it.SetBoundary(4);
   for( it.First(); it.Valid(); it.Next() )
   {
     v->Cell( it ) = 0.0 - v->Cell(it.Right() );
   }
}


/// Updates the pressure field p
void
Geometry::Update_P
(
   Grid *p
) const
{
   BoundaryIterator it( this );


   it.SetBoundary(1);
   for( it.First(); it.Valid(); it.Next() )
   {
     p->Cell( it ) = p->Cell(it.Top());
   }

   it.SetBoundary(2);
   for( it.First(); it.Valid(); it.Next() )
   {
     p->Cell( it ) = p->Cell(it.Left());
   }

   it.SetBoundary(3);
   for (it.First(); it.Valid(); it.Next()) {
     p->Cell(it) = p->Cell(it.Down());
   }

   it.SetBoundary(4);
   for( it.First(); it.Valid(); it.Next() )
   {
     p->Cell( it ) = p->Cell(it.Right());
   }
}
