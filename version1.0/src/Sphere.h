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

/*! \file Sphere.h
 * \brief 3d spheres
 * */


#ifndef SPHERE_H_
#define SPHERE_H_

#include <vectorNdim.h>

class Sphere {
public:
	Vec3d r; // center of mass
	Vec3d v; // velocity

	double radius, radius2;
	double m; // mass
	double D; // diffusion constant

    double t;        // local time of the particle
    long int ncoll;  // number of collisions up to the local time

	// creator/distructor	
	Sphere(): radius(0.5), radius2(0.25){}
	Sphere(double const &Rmax): radius(Rmax), radius2(Rmax*Rmax) {}
	Sphere(Vec3d const &pos, double const &Rmax): 
		r(pos), radius(Rmax), radius2(Rmax*Rmax) {}
	Sphere(Vec3d const &pos, Vec3d const &vel): 
		r(pos), v(vel), radius(0.5), radius2(0.25) {}
	Sphere(Vec3d const &pos, Vec3d const &vel, double const &Rmax): 
		r(pos), v(vel), radius(Rmax), radius2(Rmax*Rmax) {}

	// input/output
	friend std::ostream& operator<<(std::ostream& output, const Sphere& s)
	{ return (output << s.r << " " << s.v << " " << s.radius); }
	friend std::istream& operator>>(std::istream& input, Sphere& s)
	{ input >> s.r  >> s.v >> s.radius; s.radius2=s.radius*s.radius; return input; }

// prototypes 

	// geometrical 
	double rmax(){return radius;}
	double area();
	double volume();
	double d_volume(){return this->volume();}
	
	// periodic conditions
	void in_box(void){ r.in_box();}
	void in_box(Vec3d const& L){ r.in_box(L);}
	
	// sphere-sphere overlap
	bool overlap(Sphere const &B);
	bool overlap(Sphere const &B, double const &sigma2);
	friend bool overlap(Sphere const &A, Sphere const &B);
	friend bool overlap(Sphere const &A, Sphere const &B, double const &sigma2);

	bool overlap(Sphere const &B, Vec3d const &L);
	bool overlap(Sphere const &B, double const &sigmaAB2, Vec3d const &L);
	friend bool overlap(Sphere const &A, Sphere const &B, Vec3d const &L);
	friend bool overlap(Sphere const &A, Sphere const &B, double const &sigma2, Vec3d const &L);

};


#endif /*SPHERE_H_*/

