#include <cstdio>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "assert.h"

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

void Geometry::drivenCavity() {
  int geom_line = _size[1] + 1;
  char line_buffer[2048];
  _flags = (BoundaryType*) malloc( (_size[0] + 2) * (_size[1] + 2) * sizeof(BoundaryType));

  for (index_t i = 0; i < _size[0]+2; ++i)
    line_buffer[i] = static_cast <char> (BoundaryType::INFLOW) + '0';
  addLine(line_buffer, geom_line);
  --geom_line;

  while(geom_line > 0) {
    line_buffer[0] = static_cast <char> (BoundaryType::NOSLIP) + '0';
    for (index_t i = 1; i < _size[0] + 1; ++i)
      line_buffer[i] = static_cast <char> (BoundaryType::FLUID) + '0';
    line_buffer[_size[0] + 1] = static_cast <char> (BoundaryType::NOSLIP) + '0';

    addLine(line_buffer, geom_line);
    --geom_line;
  }
  for (index_t i = 0; i < _size[0] + 2; ++i)
    line_buffer[i] = static_cast <char> (BoundaryType::NOSLIP) + '0';
  addLine(line_buffer, 0);
}

void Geometry::addLine(char* line_buffer, int geom_line ){
  Iterator it(this, (_size[0]+2) * geom_line);
  for(index_t i = 0; i < _size[0] + 2; ++i) {
    char flag = line_buffer[i] - '0';
    _flags[it] = static_cast< BoundaryType >( flag );
    switch( static_cast< BoundaryType >( flag ) )
		  {
      case BoundaryType::FLUID :
        _FLUID.push_back( it );
        break;
      case BoundaryType::OUTFLOW :
        _OUTFLOW.push_back( it );
        break;
      case BoundaryType::INFLOW :
        _INFLOW.push_back( it );
        break;
      case BoundaryType::NOSLIP :
        _NOSLIP.push_back( it );
        break;
      case BoundaryType::SLIP :
        _SLIP.push_back( it );
        break;
      case BoundaryType::OBSTACLE :
        _OBSTACLE.push_back( it );
        break;
      case BoundaryType::PRESSURE :
        _PRESSURE.push_back( it );
		  }
    // printf("Stelle %d: %d\n", it.Value(), _flags[it]);
    // std::cout << "flag:" << (int) _flags[it] << std::endl;
    it.Next();
  }
}

void Geometry::Load(const char *file){
  FILE* handle = fopen(file,"r");
  double inval[2] = {0.0, 0.0};
  char name[20];
  char line_buffer[2048];
  int geom_line = -2;
  while (!feof(handle)) {
    while(geom_line >= 0) {
      fscanf(handle, "%s\n", line_buffer);
      addLine(line_buffer, geom_line);
      --geom_line;
      // continue;
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
  _h = {_length[0] / _size[0], _length[1] / _size[1]};
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
     u->Cell( it ) = u->Cell(it.Top()) ;
   }
   it.SetBoundary(2);
   for( it.First(); it.Valid(); it.Next() )
   {
     u->Cell( it ) = 0.0;
     u->Cell( it.Left() ) = 0.0;
   }
   it.SetBoundary(3);
   for (it.First(); it.Valid(); it.Next()) {
     u->Cell(it) = u->Cell(it.Down());
   }

   it.SetBoundary(4);
   for( it.First(); it.Valid(); it.Next() )
  {
    u->Cell( it ) = 0.0;
  }

}




void
Geometry::Update_GF
(
   Grid* F,
   const Grid* u,
   Grid* G,
   const Grid* v
) const
{
	for( const auto it : _INFLOW )
	{
      assert( it.Right().Valid() );
		F->Cell( it ) = u->Cell( it );
		G->Cell( it ) = v->Cell( it );
	}

	for( const auto it : _OUTFLOW )
	{
      assert( it.Left().Valid() );
		F->Cell( it ) = u->Cell( it );
		G->Cell( it ) = v->Cell( it );
	}


	for( const auto it : _NOSLIP )
	{
		F->Cell( it ) = u->Cell( it );
		G->Cell( it ) = v->Cell( it );
	}
}



void
Geometry::Update_UVP
(
   Grid* u,
   Grid* v,
   Grid* p
) const
{
	for( const auto it : _INFLOW )
	{
      assert( it.Right().Valid() );
		u->Cell( it ) = _velocity[0];
		v->Cell( it ) = 2.0*_velocity[1] - v->Cell( it.Right() );
		p->Cell( it ) = 2.0*_pressure - p->Cell( it.Right() );
	}

	for( const auto it : _OUTFLOW )
	{
      assert( it.Left().Valid() );
		u->Cell( it ) = u->Cell( it.Left() );
		v->Cell( it ) = v->Cell( it.Left() );
		p->Cell( it ) = - p->Cell( it.Left() );
	}


  //NOSLIP haben genau 1 FLUID Zelle angrenzend
	for( const auto it : _NOSLIP )
	{
		const bool isTopFLUID = it.Top().Valid() && BoundaryType::FLUID == Flag( it.Top() );
		const bool isDownFLUID = it.Down().Valid() && BoundaryType::FLUID == Flag( it.Down() );
		const bool isLeftFLUID = it.Left().Valid() && BoundaryType::FLUID == Flag( it.Left() );
		const bool isRightFLUID = it.Right().Valid() && BoundaryType::FLUID == Flag( it.Right() );
		if( isTopFLUID ) // e.g. under boundary
		{
      u->Cell(it) = - u->Cell(it.Top());
      v->Cell(it) = 0.0;
      p->Cell( it ) = p->Cell(it.Top());
		}
		else if( isDownFLUID ) // e.g. upper boundary
		{
      u->Cell( it ) = - u->Cell( it.Down() );
      v->Cell( it ) = 0.0;
      v->Cell( it.Down() ) = 0.0;
      p->Cell( it ) = p->Cell( it.Down() );
    }
		else if( isLeftFLUID )
		{
      u->Cell( it ) = 0.0;
      u->Cell( it.Left() ) = 0.0;
      v->Cell( it ) = - v->Cell( it.Left() );
      p->Cell( it ) = p->Cell( it.Left() );

		}
		else if( isRightFLUID )
		{
      u->Cell( it ) = 0.0;
      v->Cell( it ) = -v->Cell( it.Right());
      p->Cell( it ) = p->Cell( it.Right() );
		}
	}


   for( const auto it : _OBSTACLE )
	{
		const bool isTopFLUID = it.Top().Valid() && BoundaryType::FLUID == Flag( it.Top() );
		const bool isDownFLUID = it.Down().Valid() && BoundaryType::FLUID == Flag( it.Down() );
		const bool isLeftFLUID = it.Left().Valid() && BoundaryType::FLUID == Flag( it.Left() );
		const bool isRightFLUID = it.Right().Valid() && BoundaryType::FLUID == Flag( it.Right() );
		if( isTopFLUID ) // e.g. under boundary
		{
      assert(!isDownFLUID);     //breite wäre nur 1 feld
			if( isLeftFLUID )
			{
				u->Cell( it ) = 0.0;
				v->Cell( it ) = 0.0;
				p->Cell( it ) = 0.5 * ( p->Cell( it.Top() ) + p->Cell( it.Left() ) );
			}
			else if( isRightFLUID )
			{
				u->Cell( it ) = 0.0;
				v->Cell( it ) = 0.0;
				p->Cell( it ) = 0.5 * ( p->Cell( it.Top() ) + p->Cell( it.Right() ) );
			}
			else // nur oben ist FLUID
			{
				u->Cell( it ) = -u->Cell( it.Top() );
				v->Cell( it ) = 0.0;
				p->Cell( it ) = p->Cell( it.Top() );
			}

		}
		else if( isDownFLUID ) // e.g. upper boundary
		{
      assert(!isTopFLUID);     //breite wäre nur 1 feld
			if( isLeftFLUID )
			{
				u->Cell( it ) = 0.0;
				v->Cell( it ) = 0.0;
				p->Cell( it ) = 0.5 * ( p->Cell( it.Down() ) + p->Cell( it.Left() ) );
			}
			else if( isRightFLUID )
			{
				u->Cell( it ) = 0.0;
				v->Cell( it ) = 0.0;
				p->Cell( it ) = 0.5 * ( p->Cell( it.Down() ) + p->Cell( it.Right() ) );
			}
			else
			{
				u->Cell( it ) = -u->Cell( it.Down() );
				v->Cell( it ) = 0.0;
				p->Cell( it ) = p->Cell( it.Down() );
			}
		}
		else if( isLeftFLUID )
		{
      assert(!isRightFLUID);     //breite wäre nur 1 feld
      u->Cell( it ) = 0.0;
			v->Cell( it ) = -v->Cell( it.Left() );
			p->Cell( it ) = p->Cell( it.Left() );
		}
		else if( isRightFLUID )
		{
      assert(!isLeftFLUID);     //breite wäre nur 1 feld
      u->Cell( it ) = 0.0;
			v->Cell( it ) = -v->Cell( it.Right() );
			p->Cell( it ) = p->Cell( it.Right() );
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
	for( const auto it : _NOSLIP )
	{
		const bool isTopFLUID = it.Top().Valid() && BoundaryType::FLUID == Flag( it.Top() );
		const bool isDownFLUID = it.Down().Valid() && BoundaryType::FLUID == Flag( it.Down() );
		const bool isLeftFLUID = it.Left().Valid() && BoundaryType::FLUID == Flag( it.Left() );
		const bool isRightFLUID = it.Right().Valid() && BoundaryType::FLUID == Flag( it.Right() );
		if( isTopFLUID ) // e.g. under boundary
		{
			if( isLeftFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Top() ) + p->Cell( it.Left() ) );
			}
			else if( isRightFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Top() ) + p->Cell( it.Right() ) );
			}
			else
			{
				p->Cell( it ) = p->Cell( it.Top() );
			}

		}
		else if( isDownFLUID ) // e.g. upper boundary
		{
			if( isLeftFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Down() ) + p->Cell( it.Left() ) );
			}
			else if( isRightFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Down() ) + p->Cell( it.Right() ) );
			}
			else
			{
				p->Cell( it ) = p->Cell( it.Down() );
			}
		}
		else if( isLeftFLUID )
		{
			p->Cell( it ) = p->Cell( it.Left() );
		}
		else if( isRightFLUID )
		{
			p->Cell( it ) = p->Cell( it.Right() );
		}
	}



	for( const auto it : _OBSTACLE )
	{
		const bool isTopFLUID = it.Top().Valid() && BoundaryType::FLUID == Flag( it.Top() );
		const bool isDownFLUID = it.Down().Valid() && BoundaryType::FLUID == Flag( it.Down() );
		const bool isLeftFLUID = it.Left().Valid() && BoundaryType::FLUID == Flag( it.Left() );
		const bool isRightFLUID = it.Right().Valid() && BoundaryType::FLUID == Flag( it.Right() );
		if( isTopFLUID ) // e.g. under boundary
		{
			if( isLeftFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Top() ) + p->Cell( it.Left() ) );
			}
			else if( isRightFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Top() ) + p->Cell( it.Right() ) );
			}
			else
			{
				p->Cell( it ) = p->Cell( it.Top() );
			}

		}
		else if( isDownFLUID ) // e.g. upper boundary
		{
			if( isLeftFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Down() ) + p->Cell( it.Left() ) );
			}
			else if( isRightFLUID )
			{
				p->Cell( it ) = 0.5 * ( p->Cell( it.Down() ) + p->Cell( it.Right() ) );
			}
			else
			{
				p->Cell( it ) = p->Cell( it.Down() );
			}
		}
		else if( isLeftFLUID )
		{
			p->Cell( it ) = p->Cell( it.Left() );
		}
		else if( isRightFLUID )
		{
			p->Cell( it ) = p->Cell( it.Right() );
		}
	}



	for( const auto it : _INFLOW )
	{
		assert( it.Right().Valid() );
		p->Cell( it ) = 2.0*_pressure - p->Cell(it.Right());
	}

	for( const auto it :_OUTFLOW )
	{
		p->Cell( it ) = -p->Cell(it.Left());
	}

/*
   BoundaryIterator it( this );


   it.SetBoundary(1);
   for( it.First(); it.Valid(); it.Next() )
   {
     p->Cell( it ) = p->Cell(it.Top());
   }

   it.SetBoundary(2);
   for( it.First(); it.Valid(); it.Next() )
   {
     p->Cell( it ) = -p->Cell(it.Left());
   }

   it.SetBoundary(3);
   for (it.First(); it.Valid(); it.Next()) {
     p->Cell(it) = p->Cell(it.Down());
   }

   it.SetBoundary(4);
   for( it.First(); it.Valid(); it.Next() )
   {

   }*/
}

BoundaryType Geometry::Flag(Iterator it) const {
//    std:: cout << "val:" << it.Value() << std::endl;
//    std:: cout << "flag:" << (int) _flags[it] << std::endl;
  return _flags[ it.Value() ];
}

void Geometry::print(){
	std::cout << std::endl;
  for (int line = _size[1]+1; line >=0; line--) {
    Iterator it(this, (_size[0]+2) * line);
    for(; it.Value() < (_size[0]+2) * (line + 1);  it.Next()) {
      std::cout << (int) this->Flag(it);
      if(it.Value() % (this->Size()[0]+2) == this->Size()[0] + 1)
        std::cout << std::endl;
    }
  }
}
