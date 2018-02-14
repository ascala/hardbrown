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
 * \brief Animated event-driven simulation of 2d shaken inelastic disks
 * */

#include <iostream>
#include <fstream>  
      
#include <vector>      
#include <map>
#include <queue>

#include <cmath>
#include <cstdlib>

#include <gaussdev.h>

#include <Disk.h>
#include <SimBox.h>
#include <EDsimul.h>

using namespace std;

const int dim=2;

SimBox<Disk,dim> box;
int N=256, Nside=16;

double phi=0.40, radius=5; 
bool square_lattice=false;

double taugranular=1.0; // granular step size
double eps=0.99; // restitution coefficient

EDsimul<Disk,dim>  simulation;

#include "wrap_plplot.h"
wrap_plstream pl;

int step4frame=1; // numbero of EDBD steps per frame
void simulstep(){

	Event currentevent; 
	for(int step=0;step<=step4frame;step++)
		do currentevent=simulation.GranularStep();
			while (currentevent.kind!=GRANULARRESET);
}

void printHelp() {
	cout<< "\nUsage : ./EDBDplplot.x [ --option option_val ]\n"
		<< "Options: \n"
  	   	<< "   --particles : Number of particles\n"
       	<< "   --phi :   packing fraction\n"
       	<< "   --radius :   radius of the particle\n"
       	<< "   --tau :   tau for periodic drive\n"
       	<< "   --eps :   restitution coefficient\n"
      	<< "   --step4frame :  simulation steps per frame\n"
       	<< "   --square :  start from square lattice\n"
       	<< "   --nside :  number of particles per side\n"
       	<< endl;
  	exit(1);
}


void checkArgs(int argc, char **argv) {
  int i=1; 
  if(argc>1) do {
  	if ( !strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") )
 		{ printHelp(); exit(0);} 
    else if (!strcmp(argv[i], "--particles"))
       { N = atoi(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--nside"))
       { Nside = atoi(argv[i + 1]); N=Nside*Nside; i+=2; }
    else if (!strcmp(argv[i], "--phi"))
       { phi = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--radius"))
       { radius = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--tau"))
       { taugranular = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--eps"))
       { eps = atof(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--step4frame"))
       { step4frame = atoi(argv[i + 1]); i+=2; }
    else if (!strcmp(argv[i], "--square"))
       { square_lattice = true; i+=1; }
    else { printHelp(); exit(1);};
  } while(i<argc);
}



int main(int argc, char **argv){
	checkArgs(argc,argv);
box.tol=1.0e-1;
	if(square_lattice) box.CreateSQ(Nside,phi,radius); 
	else box.CreateRND(N,phi,radius);
box.tol=1.0e-20;
	box.ThermalizeVel();  
	
  	pl.draw_init(1.1*box.L);

	simulation.init(box);
	simulation.timenow=0;
  	simulation.InitGranularEvents(taugranular,eps);	
  	
  	while(!simulation.EventQueue.empty()){
  		
		simulstep();
    	Disk *c=box.particles;
		pl.clear();
		simulation.PropagateAll();
		double radius=5;
		for(int i=0;i<box.Npart;++i) {
			double ri=c[i].radius; int col;
			if(ri>radius) col=3; else if(ri<radius)col=5; else col=7;
			pl.draw(c[i],col);
		}

  	}

}


