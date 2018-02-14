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

/*! \file Disk.h
 * \brief 2d disks
 * */


#ifndef DISK_H_
#define DISK_H_

#include <vectorNdim.h>

class Disk {
public:
	Vec2d r; //! center of mass
	Vec2d v; //! velocity

	double radius, radius2; 
	double m; //! mass
	double D; //! diffusion constant

    double t;        // local time of the particle
    long int ncoll;  // number of collisions up to the local time

	// creator/distructor	
	Disk(): radius(0.5), radius2(0.25){}
	Disk(double const &Rmax): radius(Rmax), radius2(Rmax*Rmax) {}
	Disk(Vec2d const &pos, double const &Rmax): 
		r(pos), radius(Rmax), radius2(Rmax*Rmax) {}
	Disk(Vec2d const &pos, Vec2d const &vel): 
		r(pos), v(vel), radius(0.5), radius2(0.25) {}
	Disk(Vec2d const &pos, Vec2d const &vel, double const &Rmax): 
		r(pos), v(vel), radius(Rmax), radius2(Rmax*Rmax) {}

	// input/output
	friend std::ostream& operator<<(std::ostream& output, const Disk& s)
	{ return (output << s.r << " " << s.v << " " << s.radius); }
	friend std::istream& operator>>(std::istream& input, Disk& s)
	{ input >> s.r  >> s.v >> s.radius; s.radius2=s.radius*s.radius; return input; }

// prototypes 

	// geometrical 
	double rmax(){return radius;}
	double perimeter();
	double area();
	double d_volume(){return this->area();}
	
	// periodic conditions
	void in_box(void){ r.in_box();}
	void in_box(Vec2d const& L){ r.in_box(L);}
	
	// disk-disk overlap
	bool overlap(Disk const &B);
	bool overlap(Disk const &B, double const &sigma2);
	friend bool overlap(Disk const &A, Disk const &B);
	friend bool overlap(Disk const &A, Disk const &B, double const &sigma2);

	bool overlap(Disk const &B, Vec2d const &L);
	bool overlap(Disk const &B, double const &sigmaAB2, Vec2d const &L);
	friend bool overlap(Disk const &A, Disk const &B, Vec2d const &L);
	friend bool overlap(Disk const &A, Disk const &B, double const &sigma2, Vec2d const &L);

};


#endif /*DISK_H_*/

