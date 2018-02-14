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

/*! \file SimBox.h
 * \brief Simulation boxes for 2d and 3d
 * */


#ifndef SIMBOX_H
#define SIMBOX_H

#include <vectorNdim.h>

#include<iostream>
#include<fstream>

//!Template class for the simulation box
/*!
 * */
template <class Particle, int dim>
class SimBox{
public:
	
	int Npart;//!< Number of particles
	double tol;//!< Numerical tolerance for the overlap
	vectorNdim<double,dim> L;//!< Sides of the box
	Particle *particles;//!< Particles in the box	

	SimBox():Npart(0),tol(1e-10){}
	SimBox(int N, vectorNdim<double,dim> Side, double mytol=1e-10):
	Npart(N),L(Side),tol(mytol){particles=new Particle[Npart];}
	virtual ~SimBox(){if(Npart>0) delete[] particles; Npart=0;}
	
	void newparticles(int N){
		if(Npart>0) delete[] particles; 
		Npart=N; particles=new Particle[Npart];
	}

        // input/output
	friend std::ostream& operator<<(std::ostream& output, const SimBox& box){
		output << box.Npart << std::endl;
    	output << box.L << std::endl;
    	for (int i = 0; i < box.Npart; ++i) 
    		output << box.particles[i] << std::endl; 
    	return output; 
	}

	friend std::istream& operator>>(std::istream& input, SimBox& box){ 
        if(box.Npart>0) { delete[] box.particles; box.Npart=0;}
        input >> box.Npart;
    	input >> box.L;
    	box.particles = new Particle[box.Npart];
    	for (int i = 0; i < box.Npart; ++i) input >> box.particles[i];
        return input; 
	}

	//! Writes a configuration from file
    bool write_conf(const char *namefile){
    	std::fstream myfile;
        myfile.precision(8);
        myfile.open(namefile,std::ios::out);
        myfile<<*this;
        myfile.close();  
        return true;
    }       

	//! Read a configuration from file
    bool read_conf(const char *namefile){
        std::fstream myfile;            
        myfile.open(namefile,std::ios::in);
        myfile >> *this;
        myfile.close();
        return true;
    }       

	// other
	
	//! Calculates kinetic energy
	double kinenergy() {
  		double ekin=0;
  		for(int i=0; i<Npart; i++) ekin+=particles[i].v.norm2();
  		return ekin;
	}
	
	//! Thermalizes the velocities
	void ThermalizeVel(double Temp=1.0){ 
		for(int i=0;i<Npart;++i) particles[i].v.gauss(Temp);
	}
	
	//! Calculates the packing fraction
	double PackingFraction(){
		double Vpart=0; 
		for(int i=0;i<Npart;++i) Vpart+=particles[i].d_volume();
		return Vpart/L.d_volume();
	}
	
	void AddField(vectorNdim<double,dim> drift){
		for(int i=0;i<Npart;++i) particles[i].v+=drift;
	}
	
	//! Creates a random configuration of particles
	void CreateRND(int N, double phi, double radius=5, double tstart=0, double Temp=1.0){

		if(Npart>0) delete[] particles; 
		Npart=N; particles=new Particle[Npart];
		
  		Particle part(radius); part.t=tstart; part.ncoll=0;
  		for(int i=0; i<Npart; i++) particles[i]=part;

  		double side=pow(Npart*part.d_volume()/phi,1.0/dim); 
std::cerr<<"phi="<<Npart*part.d_volume()/pow(side,dim)<<std::endl<<std::flush;
		L=vectorNdim<double,dim>(side);
		
		long int maxattempts=10000*Npart, attempts=0;
  		for(int i=0; i<Npart; i++){
			if(!(i % 100 )&&(i!=0)) 
				std::cerr << i << " particles inserted "<< std::endl << std::flush;
  			bool any_overlap;
  			do{
  				particles[i].r.rnd(side/2-particles[i].radius-tol); 
  				any_overlap=false;
  				for(int j=0; j<i; ++j) {
					double Sigma=particles[i].radius+particles[j].radius+tol, Sigma2=Sigma*Sigma;
  					any_overlap=any_overlap || overlap(particles[i],particles[j],Sigma2,L);
				}
				attempts++;
				if(attempts>=maxattempts){
    				std::cerr<<"random insertions did not work!!!\n"<<std::flush;
    				exit(0);
    			}
    		} while(any_overlap);
  		}
		std::cerr << Npart << " particles inserted "<< std::endl << std::flush;

	}

	//! Creates a Square 2d lattice
    void CreateSQ(int Nside, double phi,double radius=5, double tstart=0, double Temp=1.0){
		if(dim!=2){
			std::cerr<<"CreateSQ needs dim=2\n"<<std::flush;
    		exit(0);
    	}
 		Particle part(radius); part.t=tstart; part.ncoll=0;
		double phimax = part.area()/(pow(2*radius,dim));
        if(phi>phimax){
        	std::cerr<<"phi too high for CreateSQ\n"<<std::flush;
    		exit(0);
    	}
    	
        Vec2d a(2*radius);
        a *= pow(phimax/phi,1.0/dim);
                
        L=Nside*a;
		if(Npart>0) delete[] particles; 
		Npart=Nside*Nside; particles=new Particle[Npart];		
 				        
        Vec2d ex(a[0],0), ey(0,a[1]);
        int count=0;
        for (int ix = 0; ix < Nside; ++ix)
        for (int iy = 0; iy < Nside; ++iy){
        	part.r=ix*ex+iy*ey+a/2;
            part.r.in_box(L);
            particles[count] = part;
            count++;
        }
        
	}
	
	//! Creates a Simple Cubic 3d lattice
    void CreateSC(int Nside, double phi,double radius=5, double tstart=0, double Temp=1.0){
		if(dim!=3){
			std::cerr<<"CreateSC needs dim=3\n"<<std::flush;
    		exit(0);
    	}

 		Particle part(radius); part.t=tstart; part.ncoll=0;
		double phimax = part.volume()/(pow(2*radius,dim));
        if(phi>phimax){
        	std::cerr<<"phi too high for CreateSC\n"<<std::flush;
    		exit(0);
    	}
        
        double sigma=2*radius;
        Vec3d a(sigma);
        a *= pow(phimax/phi,1.0/dim);
                
        L=Nside*a;
		if(Npart>0) delete[] particles; 
		Npart=Nside*Nside*Nside; particles=new Particle[Npart];
		
				        
        Vec3d ex(a[0],0,0), ey(0,a[1],0), ez(0,0,a[2]);
        int count=0;
        for (int ix = 0; ix < Nside; ++ix)
        for (int iy = 0; iy < Nside; ++iy)
        for (int iz = 0; iz < Nside; ++iz){
        	part.r=ix*ex+iy*ey+iz*ez+a/2;
            part.r.in_box(L);
            particles[count] = part;
            count++;
        }
        
	}

	//! Creates a Face Centered Cubic 3d lattice
    void CreateFCC(int Nside, double phi,double radius=5, double tstart=0, double Temp=1.0){
		if(dim!=3){
			std::cerr<<"CreateSC needs dim=3\n"<<std::flush;
    		exit(0);
    	}

 		Particle part(radius); part.t=tstart; part.ncoll=0;
        
        double sigma=2*radius, a=sqrt(2)*sigma;
        int point4cell=4;
		double phimax = point4cell*part.volume()/(a*a*a);
        if(phi>phimax){
        	std::cerr<<"phi too high for CreateFCC\n"<<std::flush;
    		exit(0);
    	}
        a *= pow(phimax/phi,1.0/dim);
        Vec3d t0(0,0,0),t1(a/2,a/2,0),t2(a/2,0,a/2),t3(0,a/2,a/2);
                
        L=Nside*a;
		if(Npart>0) delete[] particles; 
		Npart=point4cell*Nside*Nside*Nside; particles=new Particle[Npart];

				        
        int count=0;
        for (int ix = 0; ix < Nside; ++ix)
        for (int iy = 0; iy < Nside; ++iy)
        for (int iz = 0; iz < Nside; ++iz){
        	Vec3d r(ix*a,iy*a,iz*a);
        	part.r=r+t0; part.r.in_box(L); particles[count] = part; count++;
        	part.r=r+t1; part.r.in_box(L); particles[count] = part; count++;
        	part.r=r+t2; part.r.in_box(L); particles[count] = part; count++;
        	part.r=r+t3; part.r.in_box(L); particles[count] = part; count++;
        }
        
	}

	//! Creates a Body Centered Cubic 3d lattice
    void CreateBCC(int Nside, double phi,double radius=5, double tstart=0, double Temp=1.0){
		if(dim!=3){
			std::cerr<<"CreateSC needs dim=3\n"<<std::flush;
    		exit(0);
    	}

 		Particle part(radius); part.t=tstart; part.ncoll=0;
        
        double sigma=2*radius, a=2*sigma/sqrt(3);
        int point4cell=2;
		double phimax = point4cell*part.volume()/(a*a*a);
        if(phi>phimax){
        	std::cerr<<"phi too high for CreateBCC\n"<<std::flush;
    		exit(0);
    	}
        a *= pow(phimax/phi,1.0/dim);
                
        L=Nside*a;
		if(Npart>0) delete[] particles; 
		Npart=point4cell*Nside*Nside*Nside; particles=new Particle[Npart];

				        
        Vec3d t0(0,0,0),t1(a/2,a/2,a/2);
        int count=0;
        for (int ix = 0; ix < Nside; ++ix)
        for (int iy = 0; iy < Nside; ++iy)
        for (int iz = 0; iz < Nside; ++iz){
        	Vec3d r(ix*a,iy*a,iz*a);
        	part.r=r+t0; part.r.in_box(L); particles[count] = part; count++;
        	part.r=r+t1; part.r.in_box(L); particles[count] = part; count++;
        }
        
	}

	//! checks if any particle overlaps
	bool checkoverlap(){
		for(int i=1; i<Npart; i++){
  			double ri=particles[i].radius;
    		for(int j=0; j<i; j++){
      			double sigma= ri+particles[j].radius;
      			vectorNdim<double,dim> rij=particles[i].r-particles[j].r;
      			rij.in_box(L);
      			if(rij.norm2()<sigma*sigma){
std::cerr <<" OVERLAP "<<i<<" "<<j<< " " << rij.norm() - sigma << std::endl << std::flush;
				return true;
      			}
    		}
  		}

  		return false;

	}
	

};

#endif //SIMBOX_H
