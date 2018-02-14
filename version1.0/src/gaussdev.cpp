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


#include<cctype>
#include<climits>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cmath>

#include<gaussdev.h>

double drand(){return (double)rand()/(1+(double)RAND_MAX);}
double gaussdev(void){return grand48();}

/* Approximation for gaussian distributed rnd numbers with zero mean and unit */
/* variance using algorithm from Knuth - "The art of computer programming"    */
double gsum(void){
  const double  a=3.949846138, b=0.252408784, c=0.076542912, 
    			d=0.008355968, e=0.029899776;

  double sum = 0.0;
  for(int i=0;i<12;++i) sum += drand();
  
  double r=(sum-6.0)/4.0, r2=r*r;
  return  ((((e*r2+d)*r2+c)*r2+b)*r2+a)*r;
}

double gsum48(void){
  const double  a=3.949846138, b=0.252408784, c=0.076542912, 
    			d=0.008355968, e=0.029899776;

  double sum = 0.0;
  for(int i=0;i<12;++i) sum += drand48();
  
  double r=(sum-6.0)/4.0, r2=r*r;
  return  ((((e*r2+d)*r2+c)*r2+b)*r2+a)*r;
}

/* Generates gaussian distributed rnd numbers with zero mean and unit variance */
/* using Box-Muller method with drand48() as the uniform rnd numbers generator */
double grand(void)  {  
  static bool have_grand=false; static double old_grand; 
  double BoxMuller,r2,x,y; 

  if (have_grand) {// Gaussian rnd number already in memory 
    have_grand=false; // will not have gaussian rnd number in memory 
    return old_grand; // return old gaussian rnd 
  } 
  else { /* No gaussian rnd number is in memory */ 

    do { 
      x=2.0*drand()-1.0; // two uniform rnd numbers in [-1,1]x[-1,1] 
      y=2.0*drand()-1.0; 
      r2=x*x+y*y;     
    } 
    while (r2 >= 1.0 || r2 == 0.0); // force x,y in the unit disk
 
    // Box-Muller transform factor
    BoxMuller=sqrt(-2.0*log(r2)/r2); 
    
    old_grand=x*BoxMuller; have_grand=true; // store one gaussian rnd number
    return y*BoxMuller; // return the other gaussian rnd number
  } 

}

double grand48(void)  {  
  static bool have_grand48=false; static double old_grand48; 
  double BoxMuller,r2,x,y; 

  if (have_grand48) {// Gaussian rnd number already in memory 
    have_grand48=false; // will not have gaussian rnd number in memory 
    return old_grand48; // return old gaussian rnd 
  } 
  else { /* No gaussian rnd number is in memory */ 

    do { 
      x=2.0*drand48()-1.0; // two uniform rnd numbers in [-1,1]x[-1,1] 
      y=2.0*drand48()-1.0; 
      r2=x*x+y*y;     
    } 
    while (r2 >= 1.0 || r2 == 0.0); // force x,y in the unit disk
 
    // Box-Muller transform factor
    BoxMuller=sqrt(-2.0*log(r2)/r2); 
    
    old_grand48=x*BoxMuller; have_grand48=true; // store one gaussian rnd number
    return y*BoxMuller; // return the other gaussian rnd number
  } 

}
