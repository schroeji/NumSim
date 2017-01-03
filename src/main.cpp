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
#include <stdlib.h>
#include "fstream"

#include "typedef.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "iterator.hpp"
#include "visu.hpp"
#include "vtk.hpp"
#include "grid.hpp"

void write_points(const Grid* grid, std::string& path) {
  Iterator it1205(grid->getGeom(), 120, 5);
  // std::cout << it1205.Pos()[0] << "x" << it1205.Pos()[1] << std::endl;
  Iterator it6464(grid->getGeom(), 64, 64);
  // std::cout << it6464.Pos()[0] << "x" << it6464.Pos()[1] << std::endl;
  Iterator it5120(grid->getGeom(), 5, 120);
  // std::cout << it5120.Pos()[0] << "x" << it5120.Pos()[1] << std::endl;
  std::ofstream f;
  f.open(path, std::ofstream::app);
  f << grid->Cell(it1205) << "," << grid->Cell(it6464) << "," << grid->Cell(it5120) << std::endl;
  f.close();
}

void write_reynolds(real_t re, std::string& path) {
  std::ofstream f;
  f.open(path, std::ofstream::app);
  f << re << std::endl;
  f.close();
}

bool file_exists (const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}
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
  // Create parameter and geometry instances with default values
  Parameter param;
  param.Load("default.param");
  Geometry geom;
  if(param.useGeo()) {
    char* path  = argv[1];
    geom.Load(path);
    printf("loaded geometry %s\n", path);
  }
  else {
    // geom.drivenCavity();
    // geom.print();
  }
  // geom.print();
  // Create the fluid solver
  Compute comp(&geom, &param);

  std::string name_prefix = "uvalues/run_";
  uint file_num = 0;
  std::string name = name_prefix + std::to_string(file_num);
  while(file_exists(name)) {
    file_num++;
    name = name_prefix + std::to_string(file_num);
  }
  write_reynolds(param.Re(), name);
  // testIterator(&geom);
#ifdef USE_DEBUG_VISU
 // Create and initialize the visualization
 Renderer visu(geom.Length(), geom.Mesh());
 // visu.Init(800, 800);
 visu.Init(4*geom.Size()[0], 4*geom.Size()[1]);
 const Grid *visugrid;
 visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU
  // VTK vtk(geom.Mesh(), geom.Size());
  bool run = true;

  // visugrid = comp.GetP();
  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend() && run) {
	  // std::cout << "start: render " << std::endl;
    write_points(comp.GetU(), name);
#ifdef USE_DEBUG_VISU
   // Render and check if window is closed
   switch (visu.Render(visugrid)) {
   case -1:
     run = false;
     break;
   case 0:
		std::cout << "getVelocity " << std::endl;
     visugrid = comp.GetVelocity();
     // visugrid = comp.GetP();
     break;
   case 1:
		std::cout << "getU " << std::endl;
     visugrid = comp.GetU();
     break;
   case 2:
		std::cout << "start:getV " << std::endl;
     visugrid = comp.GetV();
     break;
   case 3:
		std::cout << "start:getP " << std::endl;
     visugrid = comp.GetP();
     break;
   default:
     break;
   };
#endif // DEBUG_VISU
	// std::cout << "end: render " << std::endl;
    // Create a VTK File in the folder VTK (must exist)

  // Create a VTK generator
	  // vtk.Init("VTK/field");
    // vtk.AddField("Velocity", comp.GetU(), comp.GetV());
    // vtk.AddScalar("Pressure", comp.GetP());
    // vtk.AddScalar("Vorticity", comp.GetVorticity());
    // vtk.AddScalar("Stream", comp.GetStream());
    // vtk.Finish();


    // std::cin.get();
    // Run a few steps
    // for (uint32_t i = 0; i < 9; ++i)
      comp.TimeStep(false);
    // comp.TimeStep(true);
  }
  return 0;
}
