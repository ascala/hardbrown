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



#ifndef _GOFR_H_
#define _GOFR_H_

#include <cmath>

template <class Particle, class MyVector>
class Gofr {
public:
	int Nbin;
	int Npart;
	MyVector L; // sides of the box
	double *hist;    // histogram for the radial distribution function g(r)
	double dr;    // width of a bin in the histogram 
	double norm; // normalization for the histogram 

	Gofr(int const& nbin, MyVector const &side, int const &npart){
		Nbin=nbin;
		L=side;
		dr=0.5*L.min()/Nbin;
		norm=0;
  		hist=new double[Nbin]; 
  		for(int bin=0;bin<Nbin;bin++) hist[bin]=0;
  		Npart=npart;
  	}
	~Gofr(){if(hist) delete[] hist;}
	
	void accumulate(Particle *part){
		for(int i=1;i<Npart;++i) for(int j=0;j<i;++j){
//std::cerr<<i<<" "<<j<<" : "<<std::flush;
			MyVector rij=part[i].r-part[j].r; rij.in_box(L);
			int bin=(int)(rij.norm()/dr);
//std::cerr<<rij.norm()<<" "<<dr<<" - "<<std::flush;
			if(bin<Nbin) hist[bin]+=2; // counts IJ and JI 
		}
//std::cerr<<std::endl<<std::flush;
		norm+=1;
	}
	
	void normalize3d(){ // RHO is the number density N/V  
		double rho=Npart/L.volume(), dr3=(4.0/3.0)*M_PI*dr*dr*dr;
		double bin_volume, n_idealgas; 
  	
  		for(int bin=0;bin<Nbin;bin++){
    		//double r = dr * (bin+0.5); // distance of the bin
    		int j = bin + 1;
    		bin_volume = (j*j*j - bin*bin*bin) * dr3;  // volume between I and I+1 
    		n_idealgas = bin_volume * rho;  // number of ideal-gas particles in the bin volume
    		// normalize hist so that hist[bin] = g(r)
    		hist[bin] /= (norm*n_idealgas*Npart);
		}
  	}

	void normalize2d(){ // RHO is the number density N/V  
		double rho=Npart/L.area(), dr2=4.0*M_PI*dr*dr;
		double bin_area, n_idealgas; 
  	
  		for(int bin=0;bin<Nbin;bin++){
    		// double r = dr * (bin+0.5); // distance of the bin
    		int j = bin + 1;
    		bin_area = (j*j - bin*bin) * dr2;  // volume between I and I+1 
    		n_idealgas = bin_area * rho;  // number of ideal-gas particles in the bin volume
    		// normalize hist so that hist[bin] = g(r)
    		hist[bin] /= (norm*n_idealgas*Npart);
		}
  	}

	friend std::ostream& operator<<(std::ostream& output, const Gofr& gofr){
    	for (int i = 0; i < gofr.Nbin; ++i) 
    		output << (i+0.5)*gofr.dr << " " << gofr.hist[i] << std::endl; 
    	return output; 
	}
	
};



#endif //_GOFR_H_

