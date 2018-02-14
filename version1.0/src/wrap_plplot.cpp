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



#include <iostream>

#include <wrap_plplot.h>

void wrap_plstream::draw_init(){
//	sdev ("gcw"); 
	sdev ("xwin"); 
	scolbg(255,255,255);
	init ();
	PLFLT xmin(-8),xmax(8),ymin(-8),ymax(8);
	env (xmin, xmax, ymin, ymax, 0, -0);
}

void wrap_plstream::draw_init(const Vec2d &L){
	//sdev ("gcw"); 
 	sdev ("xwin"); 
    scolbg(255,255,255);
    init ();
    PLFLT side(L[0]); if(L[1]>side) side=L[1];
    env (-side/2, side/2, -side/2, side/2, 0, -0);
}


void wrap_plstream::draw_cross(Vec2d P, int col){
	PLFLT x0,y0,x1,y1;
	PLFLT eps=0.1;
	x0 = x1 = (PLFLT)P[0]; 
	y0=(PLFLT)P[1]-eps; y1=(PLFLT)P[1]+eps;
	col0 (1); join (x0,y0,x1,y1);
	y0 = y1 = (PLFLT)P[1]; 
	x0=(PLFLT)P[0]-eps; x1=(PLFLT)P[0]+eps;
	col0 (col); join (x0,y0,x1,y1);
}

void wrap_plstream::draw_square(Vec2d P, int col){
	PLFLT x0,y0,x1,y1;
	PLFLT eps=0.1; col0 (col); 
	x0=(PLFLT)P[0]-eps; x1=(PLFLT)P[0]+eps; 
	y0=(PLFLT)P[1]-eps; y1=(PLFLT)P[1]-eps;
	join (x0,y0,x1,y1);
	
	y0=(PLFLT)P[1]+eps; y1=(PLFLT)P[1]+eps;
	join (x0,y0,x1,y1);

	y0=(PLFLT)P[1]-eps; y1=(PLFLT)P[1]+eps;
	x0=(PLFLT)P[0]-eps; x1=(PLFLT)P[0]-eps; 
	join (x0,y0,x1,y1);
	
	x0=(PLFLT)P[0]+eps; x1=(PLFLT)P[0]+eps;
	join (x0,y0,x1,y1);
}

void wrap_plstream::draw_x(Vec2d P, int col){
        PLFLT x0,y0,x1,y1;
        PLFLT eps=0.1; col0 (col);
        x0=(PLFLT)P[0]-eps; x1=(PLFLT)P[0]+eps; 
        y0=(PLFLT)P[1]-eps; y1=(PLFLT)P[1]+eps;
        join (x0,y0,x1,y1);
        y0=(PLFLT)P[1]+eps; y1=(PLFLT)P[1]-eps;
        join (x0,y0,x1,y1);
}

void wrap_plstream::draw_char(Vec2d P, char c, int col){
	PLFLT x[1]={(PLFLT)P[0]}, y[1]={(PLFLT)P[1]};
	 col0 (col); poin (1, x, y, c);
}

void wrap_plstream::draw_box(Vec2d L, int col){
	PLFLT dx(L[0]/2),dy(L[1]/2); col0(col);
	join(-dx,dy,dx,dy); join(dx,dy,dx,-dy); 
	join(dx,-dy,-dx,-dy); join(-dx,-dy,-dx,dy); 
}


void wrap_plstream::draw(Vec2d A, Vec2d B, int col){
	col0(col); join ((PLFLT)A[0], (PLFLT)A[1], (PLFLT)B[0], (PLFLT)B[1]);
}

void wrap_plstream::draw_disk(double Cx, double Cy, double Cr, int col){
  PLINT n=21;  PLFLT x[21], y[21]; 
  double alpha; col0(col);
  for(int i = 0; i<(int)n-1;i++){
  	alpha=2*M_PI*i/(int)(n-1);
    x[i] = (PLFLT) (Cx+Cr*sin(alpha));
    y[i] = (PLFLT) (Cy+Cr*cos(alpha));
  }
  x[n-1] = x[0]; y[n-1] = y[0];
  col0(col); line (n, x, y);
}

void wrap_plstream::draw(Disk const &C, int col){
  const PLINT n=41;  PLFLT x[n], y[n]; 
  double alpha; col0(col);
  for(int i = 0; i<(int)n-1;i++){
  	alpha=2*M_PI*i/(int)(n-1);
    x[i] = (PLFLT) (C.r[0]+C.radius*cos(alpha));
    y[i] = (PLFLT) (C.r[1]+C.radius*sin(alpha));
  }
  x[n-1] = x[0]; y[n-1] = y[0];
  col0(col); line (n, x, y);
}

