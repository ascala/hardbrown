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

#include <Disk.h>

double Disk::perimeter(){return 2*M_PI*radius;}

double Disk::area(){return M_PI*radius*radius;}
		
bool Disk::overlap(Disk const &B){
	double sigma2=radius+B.radius; sigma2*=sigma2;
	Vec2d dr(r[0]-B.r[0],r[1]-B.r[1]);
	return (dr.norm2()<sigma2);
}

bool Disk::overlap(Disk const &B, double const &sigma2){
	Vec2d dr(r[0]-B.r[0],r[1]-B.r[1]); 
	return (dr.norm2()<sigma2);
}

bool overlap(Disk const &A, Disk const &B, double const &sigma2){
	Vec2d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	return (dr.norm2()<sigma2);
}

bool overlap(Disk const &A, Disk const &B){
	double sigma2=A.radius+B.radius; sigma2*=sigma2;
	Vec2d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	return (dr.norm2()<sigma2);
}


bool Disk::overlap(Disk const &B, double const &sigma2, Vec2d const &L){
	Vec2d dr(r[0]-B.r[0],r[1]-B.r[1]); 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

bool Disk::overlap(Disk const &B, Vec2d const &L){
	double sigma2=radius+B.radius; sigma2*=sigma2;
	Vec2d dr(r[0]-B.r[0],r[1]-B.r[1]); 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

bool overlap(Disk const &A, Disk const &B, double const &sigma2, Vec2d const &L){
	Vec2d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

bool overlap(Disk const &A, Disk const &B, Vec2d const &L){
	double sigma2=A.radius+B.radius; sigma2*=sigma2;
	Vec2d dr(A.r[0]-B.r[0],A.r[1]-B.r[1]); 
	dr.in_box(L);
	return (dr.norm2()<sigma2);
}

