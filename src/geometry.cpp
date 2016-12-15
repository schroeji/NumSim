#include <cstdio>
#include <iostream>
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
  char line_buffer[2048];
  int geom_line = -2;
  // char
  while (!feof(handle)) {
    while(geom_line >= 0) {
      fscanf(handle, "%s\n", line_buffer);
      // printf("line %d: %s\n", geom_line, line_buffer);
      Iterator it(this, (_size[0]+2) * geom_line);
      for(index_t i = 0; i < _size[0] + 2; ++i){
        char flag = line_buffer[i] - '0';
        _flags[it] = (BoundaryType) flag;
        // printf("Stelle %d: %d\n", it.Value(), _flags[it]);
        // std::cout << "flag:" << (int) _flags[it] << std::endl;
        it.Next();
      }
      Iterator it40(this, 40);
      --geom_line;
      continue;
    }
    if (!fscanf(handle, "%s =", name)) continue;
    // printf("%s\n", name);
    if (strcmp(name,"size") == 0) {
      if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _size[0] = inval[0] - 2; // - 2 weil das erzeugte Gitter automatisch um 2 größer ist (siehe Grid)
        _size[1] = inval[1] - 2;
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
    if (strcmp(name,"geometry_start") == 0 && geom_line == -2) {
      fscanf(handle, " %s\n", line_buffer);
      _flags = (BoundaryType*) malloc( (_size[0] + 2) * (_size[1] + 2) * sizeof(BoundaryType));
      geom_line = _size[1] + 1;
      continue;
    }
  }
  fclose(handle);
  Iterator it40(this, 40);
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


// Updates the velocity field u
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
     u->Cell( it.Left() ) = 0.0;
   }
   it.SetBoundary(3);
   for (it.First(); it.Valid(); it.Next()) {
     u->Cell(it) = 2*_velocity[0] - u->Cell(it.Down());
   }

   it.SetBoundary(4);
   for( it.First(); it.Valid(); it.Next() )
     {
       u->Cell( it ) = 0.0;
     }

}


// Updates the velocity field v
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
      v->Cell(it.Down()) = 0.0;
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

BoundaryType Geometry::Flag(Iterator it){
  // std:: cout << "val:" << it.Value() << std::endl;
  // std:: cout << "flag:" << (int) _flags[it] << std::endl;
  return _flags[it];
}

void Geometry::print(){
  for (index_t line = _size[1]+1; line >=0; line--) {
    Iterator it(this, (_size[0]+2) * line);
    for(it; it.Value() < (_size[0]+2) * (line + 1);  it.Next()) {
      std::cout << (int) this->Flag(it);
      if(it.Value() % (this->Size()[0]+2) == this->Size()[0] + 1)
        std::cout << std::endl;
    }
  }
}
