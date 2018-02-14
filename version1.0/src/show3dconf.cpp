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
 * \brief shows 3d configuration using Open Inventor
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


using namespace std;

const int dim=3;

SimBox<Sphere,dim> box;
#include <InventorShow.h>
//needed from InventorShow
void simulstep(){}

int main(int argc, char **argv){

	if(argc<2) {cerr<<"Input file?\n"<<endl<<flush; exit(0);}
  	box.read_conf(argv[1]);

  	InventorShow();

}


