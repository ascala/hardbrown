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



#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <Sphere.h>

double Sphere::area(){return 4*M_PI*radius*radius;}

double Sphere::volume(){return 4*M_PI*radius*radius*radius/3;}
		
bool Sphere::overlap(Sphere const &B){
	double sigma2=radius+B.radius; sigma2*=sigma2;
	Vec3d dr=r-B.r;
	return (dr.norm2()<sigma2);
}

bool Sphere::overlap(Sphere const &B, double const &sigma2){
	Vec3d dr=r-B.r; 
	return (dr.norm2()<sigma2);
}

bool overlap(Sphere const &A, Sphere const &B, double const &sigma2){
	Vec3d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	return (dr.norm2()<sigma2);
}

bool overlap(Sphere const &A, Sphere const &B){
	double sigma2=A.radius+B.radius; sigma2*=sigma2;
	Vec3d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	return (dr.norm2()<sigma2);
}


bool Sphere::overlap(Sphere const &B, double const &sigma2, Vec3d const &L){
	Vec3d dr=r-B.r; 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

bool Sphere::overlap(Sphere const &B, Vec3d const &L){
	double sigma2=radius+B.radius; sigma2*=sigma2;
	Vec3d dr=r-B.r; 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

bool overlap(Sphere const &A, Sphere const &B, double const &sigma2, Vec3d const &L){
	Vec3d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

bool overlap(Sphere const &A, Sphere const &B, Vec3d const &L){
	double sigma2=A.radius+B.radius; sigma2*=sigma2;
	Vec3d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

