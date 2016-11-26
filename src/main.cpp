/*
 *  Copyright (C) 2015   Malte Brunn
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include "typedef.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "iterator.hpp"
#include "visu.hpp"
#include "vtk.hpp"
#include "comm.hpp"

void testIterator(Geometry *geom) {
  Iterator it(geom);
  std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  it.Next();
  std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  it = it.Right();
  std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  it = it.Top();
  std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  it = it.Down();
  it = it.Down();
  std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  std::cout << it.Value() << std::endl;

  // for(it.First(); it.Valid(); it.Next()){
     // std::cout << it.Pos()[0] << ";" << it.Pos()[1] << std::endl;
  // }

  InteriorIterator intit(geom);
  intit.First();
  std::cout << intit.Pos()[0] << ";" << intit.Pos()[1] << std::endl;
  // for(intit.First(); intit.Valid(); intit.Next()){
    // std::cout << intit.Pos()[0] << ";" << intit.Pos()[1] << std::endl;
  // }

  BoundaryIterator boundit(geom);
  for (int boundary = 1; boundary <= 4; boundary++ ) {
    printf("------ Boundary %d ------\n", boundary);
    boundit.SetBoundary(boundary);
    for(boundit.First(); boundit.Valid(); boundit.Next()) {
      std::cout << boundit.Pos()[0] << ";" << boundit.Pos()[1] << std::endl;
    }
  }
}

int main(int argc, char **argv) {

  Communicator communicator( &argc, &argv);

  // Create parameter and geometry instances with default values
  Parameter param;
  Geometry geom;
  // Create the fluid solver
  Compute comp( &geom, &param, &communicator );
  // testIterator(&geom);
#ifdef USE_DEBUG_VISU
 // Create and initialize the visualization
 Renderer visu(geom.Length(), geom.Mesh());
 visu.Init(800, 800);
#endif // USE_DEBUG_VISU

  // Create a VTK generator
  VTK vtk(geom.Mesh(), geom.Size());
  const Grid *visugrid;
  bool run = true;
  visugrid = comp.GetVelocity();

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend() && run) {
#ifdef USE_DEBUG_VISU
   // Render and check if window is closed
   switch (visu.Render(visugrid)) {
   case -1:
     run = false;
     break;
   case 0:
     visugrid = comp.GetVelocity();
     break;
   case 1:
     visugrid = comp.GetU();
     break;
   case 2:
     visugrid = comp.GetV();
     break;
   case 3:
     visugrid = comp.GetP();
     break;
   default:
     break;
   };
#endif // DEBUG_VISU

    // Create a VTK File in the folder VTK (must exist)
    vtk.Init("/home/alex/git/NumSim/NumSim-master/VTK/field");
    vtk.AddField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddScalar("Pressure", comp.GetP());

    vtk.Finish();


    // std::cin.get();
    // Run a few steps
    for (uint32_t i = 0; i < 9; ++i)
      comp.TimeStep(false);
    comp.TimeStep(true);
    // comp.TimeStep(true);
  }
  return 0;
}
