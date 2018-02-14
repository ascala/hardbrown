#ifndef MY_PL_H_
#define MY_PL_H_


#include <cmath>
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

/*! \file wrap_plplot.h
 * \brief wrapper for PlPlot
 * */



#include <plstream.h>
#include <vectorNdim.h>
#include <Disk.h>

class wrap_plstream :  public plstream {
public:
void draw_init();

void draw_init(const Vec2d &L);

void draw_cross(Vec2d P, int col=3);

void draw_square(Vec2d P, int col=3);

void draw_x(Vec2d P, int col=3);

void draw_char(Vec2d P, char c, int col=3);

void draw_box(Vec2d L, int col=3);

void draw(Vec2d A, Vec2d B, int col=3);

void draw_disk(double Cx, double Cy, double Cr, int col=3);

void draw(Disk const &C, int col=3);

private:
};


#endif /*MY_PL_H_*/
