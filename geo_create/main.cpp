#include "../src/typedef.hpp"

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "iostream"
#include "fstream"

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
  f << "velocity = " << "0.0 0.0" << endl;
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

void write_parameterfile(string path) {
  real_t re = 1000;
  real_t omega = 1.7;
  real_t alpha = 0.9;
  real_t dt = 0.5;
  real_t tend = 50;
  real_t itermax = 200;
  real_t eps = 0.0001;
  real_t tau = 0.5;

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
  f << "useGeometry = 1" << alpha << endl;
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
    write_parameterfile("default.param");
  }
  else if(!strcmp(argv[1], "all")) {
    if(argc == 2) {
      write_channel(length[0], length[1], size[0], size[1], deltaP, "channel.geom");
      write_step(length[0], length[1], size[0], size[1], deltaP, "step.geom");
      write_karman(length[0], length[1], size[0], size[1], deltaP, "karman.geom");
      write_parameterfile("default.param");
    }
    else if(argc != 7) {
      printf("%d Parameter erhalten, aber 5 erwartet für alle Geometrien.\n", argc - 2);
      print_usage(argv[0]);
    }
    else{
      write_channel(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "channel.geom");
      write_step(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "stp.geom");
      write_karman(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), "karman.geom");
      write_parameterfile("default.param");
    }
  }
  else {
    print_usage(argv[0]);
  }

  return 0;
}
