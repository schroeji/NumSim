#include "../src/typedef.hpp"

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "iostream"
#include "fstream"
#include "random"
#include "chrono"

using namespace std;

multi_index_t size = {160, 32};
multi_real_t length = {5.0, 1.0};
real_t deltaP = 0.1;

void print_line() {
  printf("-----------------------------------------------\n");
}
void print_usage(char* name) {

  print_line();
  printf("Verwendung:\n");
  print_line();

  printf("Kanal-Geometrie:\n");
  printf("%s channel xLength yLength iMax jMax deltaP\n", name);
  print_line();

  printf("Stufen-Geometrie:\n");
  printf("%s step xLength yLength iMax jMax deltaP\n", name);
  printf("Für ungerade jMax wird jMax = jMax - 1 verwendet.\n");
  print_line();

  printf("Karman-Geometrie:\n");
  printf("%s Karman xLength yLength iMax jMax deltaP\n", name);
  printf("Der Winkel alpha ist auf 45° gesetzt.\n");
  printf("Für ungerade jMax wird jMax = jMax - 1 verwendet.\n");
  printf("Falls jMax dann nicht noch durch 4 teilbar ist, wird der Teil unter dem Balken ergänzt.\n");
  print_line();

  printf("Parameterdatei:\n");
  printf("%s param\n", name);
  print_line();

  printf("Alle Geometrien:\n");
  printf("%s all xLength yLength iMax jMax deltaP\n", name);
  printf("Damit werden alle drei Geometrien mit den angegebenen Parametern erzeugt.\n");
  printf("%s all\n", name);
  printf("Erzeugt alle Geometrien mit length=%f %f; size=%d %d; deltaP=%f und die Parameterdatei.\n", length[0], length[1], size[0], size[1], deltaP);

}

void write_parameters(real_t xLength, real_t yLength, int iMax, int jMax, real_t deltaP, ofstream& f){
  f << "size = " << iMax << " " << jMax << endl;
  f << "length = " << xLength << " " << yLength << endl;
  f << "velocity = " << "1.0 0.0" << endl;
  f << "pressure = " << deltaP << endl;
  f << "geometry_start = true" << endl;
}

void write_channel(real_t xLength, real_t yLength, int iMax, int jMax, real_t deltaP, string path) {
  ofstream f;
  f.open(path);
  write_parameters(xLength, yLength, iMax, jMax, deltaP, f);
  // Parameter Ausgabe
  //Geometry Ausgabe
  for (int i = 0; i < iMax; i++) {
    f << (int) BoundaryType::NOSLIP;
  }
  f << endl;

  for (int j = 0; j < jMax - 2; j++) {
    f << (int) BoundaryType::INFLOW;
    for (int i = 0; i < iMax - 2; i++) {
      f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OUTFLOW << endl;
  }

  for (int i = 0; i < iMax; i++) {
    f << (int) BoundaryType::NOSLIP;
  }
  f << endl;

  f.close();
}


void write_step(real_t xLength, real_t yLength, int iMax, int jMax, real_t deltaP, string path) {
  // falls jMax ungerade: jMax = jMax -1
  if(jMax % 2 == 1)
    jMax--;
  int stepsize = jMax/2 - 1;

  ofstream f;
  f.open(path);
  // Parameter Ausgabe
  write_parameters(xLength, yLength, iMax, jMax, deltaP, f);
  //Geometry Ausgabe
  for (int i = 0; i < iMax; i++) {
    f << (int) BoundaryType::NOSLIP;
  }
  f << endl;

  //Teil über der Stufe
  for (int j = 0; j < stepsize; j++) {
    f << (int) BoundaryType::INFLOW;
    for (int i = 0; i < iMax - 2; i++) {
      f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OUTFLOW << endl;
  }

  //Teil mit Stufe
  for (int j = 0; j < stepsize; j++) {
    f << (int) BoundaryType::NOSLIP;
    for (int i = 1; i < iMax - 1; i++) {
      if (i < stepsize + 1)
        f << (int) BoundaryType::OBSTACLE;
      else
        f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OUTFLOW << endl;
  }


  for (int i = 0; i < iMax; i++) {
    f << (int) BoundaryType::NOSLIP;
  }
  f << endl;

  f.close();
}

void write_karman(real_t xLength, real_t yLength, int iMax, int jMax, real_t deltaP, string path) {
  // falls jMax ungerade: jMax = jMax -1
  if(jMax % 2 == 1)
    jMax--;
  int stepsize = jMax/2;
  // -1 jeweils wegen der oberen und unteren noslip boundary
  int topspace = jMax/4 - 1;
  int botspace = jMax - stepsize - topspace - 2;
  ofstream f;
  f.open(path);
  // Parameter Ausgabe
  write_parameters(xLength, yLength, iMax, jMax, deltaP, f);
  //Geometry Ausgabe
  for (int i = 0; i < iMax; i++) {
    f << (int) BoundaryType::NOSLIP;
  }
  f << endl;

  //Teil über Balken
  for (int j = 0; j < topspace; j++) {
    f << (int) BoundaryType::INFLOW;
    for (int i = 0; i < iMax - 2; i++) {
      f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OUTFLOW << endl;
  }

  // Balken
  for (int j = 0; j < stepsize; j++) {
    int i;
    //Flüssigkeit vor dem Balken
    int fluid_len = 2*stepsize - j - 1;

    //Korrektur weil nur max 2 Fluidzellen nachbarn sein dürfen
    if(j==stepsize - 2)
      --fluid_len;

    f << (int) BoundaryType::INFLOW;

    for (i = 1; i < fluid_len; i++) {
      f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OBSTACLE << (int) BoundaryType::OBSTACLE;
    //Korrektur weil nur max 2 Fluidzellen nachbarn sein dürfen
    if(j == 1 || j == stepsize - 2){
      f << (int) BoundaryType::OBSTACLE;
      ++i;
    }
    for (i = i+2; i < iMax - 1; i++){
      f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OUTFLOW;
    f << endl;
  }

  //Teil unter Balken
  for (int j = 0; j < botspace; j++) {
    f << (int) BoundaryType::INFLOW;
    for (int i = 0; i < iMax - 2; i++) {
      f << (int) BoundaryType::FLUID;
    }
    f << (int) BoundaryType::OUTFLOW << endl;
  }


  for (int i = 0; i < iMax; i++) {
    f << (int) BoundaryType::NOSLIP;
  }
  f << endl;

  f.close();
}

void write_default(int size) {
  ofstream f;
  f.open("default.geom");
  write_parameters(1.0, 1.0, size, size, 0.0, f);
  f.close();
}

void write_parameterfile(real_t re, real_t tend, real_t dt, int solver, int itermax,  string path) {
  // real_t re = 1000;
  real_t omega = 1.7;
  real_t alpha = 0.9;
  // real_t dt = 0.5;
  real_t eps = 0.001;
  real_t tau = 0.9;

  ofstream f;
  f.open(path);

  f << "re = " << re << endl;
  f << "omega = " << omega << endl;
  f << "alpha = " << alpha << endl;
  f << "dt = " << dt << endl;
  f << "tend = " << tend << endl;
  f << "itermax = " << itermax << endl;
  f << "eps = " << eps << endl;
  f << "tau = " << tau << endl;
  f << "solver = " << solver << endl;
  f << "useGeometry = 0" << endl;
}

void run_monte_carlo() {
  real_t my = 1500;
  real_t delta = 1000.0/6.0;
  real_t dt = 0.004;            // Reynoladszahl sollte mit diesem dt=0.004 nicht kleiner 400 sein
  // real_t dx = 1/128.0;
  // real_t dy = 1/128.0;
  real_t re = 0;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  normal_distribution<real_t> distrib(my, delta);
  re = std::max(distrib(generator), 400.0);
  // const real_t conv_cond = std::min(dx/1, dy/1);
  // const real_t diff_cond = (dx*dx * dy*dy*re)/(2*dx*dx + 2*dy*dy);
	// const real_t dt_safe = std::min(diff_cond*0.8, conv_cond * 0.8);
  // std::cout << dt_safe << std::endl;
  write_parameterfile(re, 50, dt, 0, 100, "default.param");
  system("./numsim montecarlo");
}



void run_uniformly_distributed
(
   const real_t re
)
{
  real_t dt = 0.004;
  write_parameterfile(re, 50, dt, 0, 100, "default.param");
  system("./numsim uniformly");
}


void run_convergence(int solver) {
  index_t sizes[] = {16, 32, 64, 128};
  real_t dt = 0.001;
  write_parameterfile(1500, 0.099, dt, solver, 0, "default.param");
  for(index_t size : sizes) {
    write_default(size);
    system("mpirun -n 1 ./NumSim");
  }
}


int main (int argc, char** argv) {
  if(argc == 1){
    print_usage(argv[0]);
  }
  else if(!strcmp(argv[1], "channel")) {
    if(argc != 7) {
      printf("%d Parameter erhalten, aber 5 erwartet für Kanal-Geometrie.\n", argc - 2);
      print_usage(argv[0]);
    }
    else{
      write_channel(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "channel.geom");
    }
  }
  else if(!strcmp(argv[1], "step")) {
    if(argc != 7) {
      printf("%d Parameter erhalten, aber 5 erwartet für Stufen-Geometrie.\n", argc - 2);
      print_usage(argv[0]);
    }
    else{
      write_step(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "step.geom");
    }
  }
  else if(!strcmp(argv[1], "karman")) {
    if(argc != 7) {
      printf("%d Parameter erhalten, aber 5 erwartet für Karman-Geometrie.\n", argc - 2);
      print_usage(argv[0]);
    }
    else{
      write_karman(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "karman.geom");
    }
  }
  else if(!strcmp(argv[1], "param") ) {
    write_parameterfile(1000.0, 50, 0.5, 0, 100, "default.param");
  }
  else if(!strcmp(argv[1], "all")) {
    if(argc == 2) {
      write_channel(length[0], length[1], size[0], size[1], deltaP, "channel.geom");
      write_step(length[0], length[1], size[0], size[1], deltaP, "step.geom");
      write_karman(length[0], length[1], size[0], size[1], deltaP, "karman.geom");
      write_default(128);
      write_parameterfile(1000.0, 50, 0.05, 3, 1, "default.param");
    }
    else if(argc != 7) {
      printf("%d Parameter erhalten, aber 5 erwartet für alle Geometrien.\n", argc - 2);
      print_usage(argv[0]);
    }
    else{
      write_channel(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "channel.geom");
      write_step(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "stp.geom");
      write_karman(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "karman.geom");
      write_parameterfile(1000.0, 50, 0.5, 3, 1, "default.param");
    }
  }
  else if(!strcmp(argv[1], "montecarlo")) {
    int runs = atoi(argv[2]);
    for(int i = 0; i < runs; i++){
      run_monte_carlo();
    }
  }
  else if(!strcmp(argv[1], "uniformly")) {
    int runs = atoi(argv[2]);
	 const real_t mu = 1500.0;
	 const real_t sigma = 1000.0/6.0;
	 const real_t step = ( 6.0 * sigma ) / ( runs - 2 );
    for(int i = mu - 3.0 * sigma; i <= mu + 3.0 * sigma ; i += step ){
      run_uniformly_distributed( i );
    }
  }
  else if(!strcmp(argv[1], "konvergenz")) {
    run_convergence(atoi(argv[2]));
  }
  else {
    print_usage(argv[0]);
  }

  return 0;
}
