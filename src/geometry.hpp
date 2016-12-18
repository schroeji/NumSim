/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "typedef.hpp"
#include <vector>
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------


/// Geometry class that loads geometry definitons from a file and is responsible for setting boundary conditions.

class Geometry {
public:
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //      u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //      u=0, v=0
  Geometry();

  /// Loads a geometry from a file
  void Load(const char *file);

  /// Returns the number of cells in each dimension
  const multi_index_t &Size() const;
  /// Returns the length of the domain
  const multi_real_t &Length() const;
  /// Returns the meshwidth
  const multi_real_t &Mesh() const;

  /// Updates the velocity field u
  void Update_U(Grid *u) const;
  /// Updates the velocity field v
  void Update_V(Grid *v) const;
  /// Updates the pressure field p
  void Update_P(Grid *p) const;
  ///Prints Parameters
  void PrintVariables();
  void print();
  BoundaryType Flag(Iterator it) const;
  void Update_UVP( Grid* u, Grid* v, Grid* p ) const;
  void Update_GF( Grid* F,const Grid* u, Grid* G , const Grid* v ) const;
  inline const std::vector< Iterator >& getFLUID( void ) const { return _FLUID;};
  void addLine(char* line_buffer, int line);
  void drivenCavity();
private:
  multi_index_t _size;
  multi_real_t _length;
  multi_real_t _h;

  multi_real_t _velocity;
  real_t _pressure;
  std::vector< Iterator > _FLUID;
  std::vector< Iterator > _NOSLIP;
  std::vector< Iterator > _SLIP;
  std::vector< Iterator > _INFLOW;
  std::vector< Iterator > _OUTFLOW;
  std::vector< Iterator > _OBSTACLE;
  BoundaryType* _flags;
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
