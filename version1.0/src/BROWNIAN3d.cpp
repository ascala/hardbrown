/*
 
Copyright 2008, Antonio Scala
 
This file is part of HardBrown: 
"HARD particles BROWNian and granular simulation".
 
HardBrown is free software:
you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
HardBrown is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with HardBrown.  If not, see <http://www.gnu.org/licenses/>.
 
*/

/*! \file
 * \brief skeleton of a 3d event driven simulation of hard 
 * Brownian spheres
 * */


#include <iostream>
#include <fstream>  
      
#include <vector>      
#include <map>
#include <queue>

#include <cmath>
#include <cstdlib>

#include <gaussdev.h>

#include <Sphere.h>
#include <SimBox.h>
#include <EDsimul.h>

#include <Gofr.h>

using namespace std;

const int dim=3;

SimBox<Sphere,dim> box;
int Nside=6, N=216; 

double phi=0.01, radius=5; 
bool cubic_lattice=false, fcc_lattice=false, bcc_lattice=false;

int nsteps=100, ntherm=0; //BD steps, thermalization steps
double taubrownian; // brownian step size

EDsimul<Sphere,dim>  simulation;


void printHelp(char **argv) {
	cerr<< "\nUsage : "<<argv[0]<<" [ --option option_val ]\n"
		<< "Options: \n"
  	   	<< "   --particles :  Number of particles\n"
       	<< "   --phi :   packing fraction\n"
       	<< "   --radius :   radius of the particle\n"
       	<< "   --tau :   tau for brownian dynamics\n"
       	<< "   --cubic :  start from cubic lattice\n"
       	<< "   --fcc :  start from fcc lattice\n"
       	<< "   --bcc :  start from fcc lattice\n"
       	<< "   --nside :  number of lattice cells per side\n"
       	<< "   --steps :  number of BD steps\n"
       	<< "   --ntherm :  number of BD thermalization steps\n"
       	<< endl;
  	exit(1);
}

void checkArgs(int argc, char **argv) {
  int i=1; 
  if (argc>1) do {
	if ( !strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") )
 		{ printHelp(argv); exit(0);}
    else if (!strcmp(argv[i], "--particles"))
       { N = atoi(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--nside"))
       { Nside = atoi(argv[i + 1]); N=Nside*Nside*Nside; i+=2; }
    else if (!strcmp(argv[i], "--steps"))
       { nsteps = atoi(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--ntherm"))
       { ntherm = atoi(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--phi"))
       { phi = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--radius"))
       { radius = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--tau"))
       { taubrownian = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--cubic"))
       { cubic_lattice = true; i+=1; }
    else if (!strcmp(argv[i], "--fcc"))
       { fcc_lattice = true; i+=1; }
    else if (!strcmp(argv[i], "--bcc"))
       { bcc_lattice = true; i+=1; }
    else { printHelp(argv); exit(1);};
  } while(i<argc);
}

int main(int argc, char **argv){
	srand(time(NULL)); srand48(time(NULL));
	
  	checkArgs(argc,argv);
	if(cubic_lattice) box.CreateSC(Nside,phi,radius); 
	else if(fcc_lattice) box.CreateFCC(Nside,phi,radius); 
	else if(bcc_lattice) box.CreateBCC(Nside,phi,radius); 
	else box.CreateRND(N,phi,radius);
	box.ThermalizeVel();  

	simulation.init(box);
	simulation.timenow=0;
  	simulation.InitBrownianEvents(taubrownian);	
	
	Gofr<Sphere,Vec3d> gofr(100,box.L,box.Npart);
	
	Event currentevent; 
	cerr<<"Thermalization: ";
	for(int step=0;step<=ntherm;step++){
		do currentevent=simulation.BrownianStep();
		while (currentevent.kind!=BROWNIANRESET);
		cerr<<step<<" - ";
	}cerr<<"done\n";
	cerr<<"Steps: ";
	for(int step=0;step<=nsteps;step++){
		do currentevent=simulation.BrownianStep();
		while (currentevent.kind!=BROWNIANRESET);
		gofr.accumulate(box.particles);
		cerr<<step<<" - ";
	}cerr<<"done\n";
	
	gofr.normalize3d();
	cout << gofr << endl << flush;
	box.write_conf("edbd3d.cnf");
}


