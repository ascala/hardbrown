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

/*! \file EDsimul.h
 * \brief  Template class containing the event queue, 
 * collision prediction, and propagation steps for Brownian, 
 * Newtonian and Granular simulations 
 * */

#ifndef EDSIMUL_H
#define EDSIMUL_H

#include <Event.h>
#include <SimBox.h>

#include <iostream>
#include <fstream>  
#include <vector>      
#include <queue>

#include <cassert>

//! Template class to manage event-driven simulations of hyperspheres in "dim" dimensions
/*!
 * */

template <class Particle, int dim>
class EDsimul{
public: 
    
	double eps; //!< restitution coefficient
	double growthrate;//!< growth rate 
		
	double tol; //!< tolerance in overlaps	
	double timeinfinite; //!< biggest time reacheable
	
	double timenow; //!< current time	
	double endtime; //!< time at which simulation ends

	double tausync;//!< period for syncing particles
	double nexttimesync;//!< next sync time 
	
	double tautherm;//!< period for thermalizing particles
	double nexttimetherm;//!< next thermalization time 

	std::priority_queue<Event> EventQueue;//!< time ordered event queue 
	
	SimBox<Particle, dim> *box;//!< Simulation box
	int Npart;//!< Number of particles
	Particle *c;//!< Particles
	vectorNdim<double,dim> L;//!< Sides of the box
	
// constructor/destructor
	EDsimul():timenow(0),Npart(0){
		eps=1.00, growthrate=0.0, tol=1e-10, 
		timeinfinite=endtime=1.0e20, 
		tausync=tautherm=1.0e-0;
	}
		
	~EDsimul(){}

// valid collisions	
	bool CollisionIsValid(int i, Event const& e){ return (c[i].ncoll==e.nci); }
	bool CollisionIsValid(int i, int j, Event const& e){ return (c[i].ncoll==e.nci)&&(c[j].ncoll==e.ncj); }

// useful functions from SimBox	
	void ThermalizeVel(double Temp=1.0){box->ThermalizeVel(Temp);}
	double PackingFraction(){return box->PackingFraction();}
	double kinenergy(){return box->kinenergy();}
	
// Initialization from SimBox	
	void init(SimBox<Particle,dim> &simbox){
		box=&simbox;
		Npart=simbox.Npart;
		c=simbox.particles;
		L=simbox.L;
		tol=simbox.tol;
	}
	

// Check if particles will collide
	bool WillCollide(Particle &ci, Particle &cj, double &tcoll, double tmax, vectorNdim<double,dim> &L){
  		double h,p,q,w;
  		double sigma=ci.radius+cj.radius;

  		double dti=timenow-ci.t, dtj=timenow-cj.t;
  		vectorNdim<double,dim> rij=ci.r+ci.v*dti-cj.r-cj.v*dtj; rij.in_box(L);
  		vectorNdim<double,dim> vij=ci.v-cj.v;
  		double vdotr=dot(rij,vij), dt;
#if 1
  		if(vdotr>0) return false;
#else // THIS DOES NOT WORK !!!
		if(overlap(ci,cj,L)) std::cerr<<"overlap"<<std::endl<<std::flush;
  		if( (vdotr>0)&&!overlap(ci,cj,L) ) return false;
#endif
  		else{
    		double dist2 = rij.norm2();
    		double Sigma= (dist2>=sigma*sigma ? sigma : sigma-tol);
    		h=1/vij.norm2(); p=vdotr*h;
    		q=(dist2-Sigma*Sigma)*h; w=p*p-q;
    		if(w<0) return false;
    		else {
      			dt=q/(-p+sqrt(w)); 
//      			assert(dt>=0);
//      			if(dt<=0) return false;
//      			if(dt<=0) dt+=drand48()*tol;
      			tcoll=timenow+dt;
      			if(tcoll>tmax) return false;
      			else return true;
    		}
  		}
  
	}

	bool WillCollideWhileGrowing(Particle &ci, Particle &cj, double &tcoll, double tmax, vectorNdim<double,dim> &L){
  		double h,p,q,w;

  		double dti=timenow-ci.t, dtj=timenow-cj.t;

  		double sigma=ci.radius+cj.radius+growthrate*(dti+dtj);
  		vectorNdim<double,dim> rij=ci.r+ci.v*dti-cj.r-cj.v*dtj; rij.in_box(L);
  		vectorNdim<double,dim> nij=rij; nij.normalize();
  		vectorNdim<double,dim> vij=ci.v-cj.v-2*growthrate*nij;
  		double vdotr=dot(rij,vij), dt;
#if 1
  		if(vdotr>0) return false;
#else // THIS DOES NOT WORK !!!
  		if((vdotr>0)&&!overlap(ci,cj,L)) return false;
#endif
  		else{
    		double dist2 = rij.norm2();
    		double Sigma= (dist2>=sigma*sigma ? sigma : sigma-tol);
    		h=1/vij.norm2(); p=vdotr*h;
    		q=(dist2-Sigma*Sigma)*h; w=p*p-q;
    		if(w<0) return false;
    		else {
      			dt=q/(-p+sqrt(w)); 
//      			assert(dt>=0);
//      			if(dt<=0) return false;
//      			if(dt<=0) dt+=drand48()*tol;
      			tcoll=timenow+dt;
      			if(tcoll>tmax) return false;
      			else return true;
    		}
  		}
  
	}

// Elastic and inelastic collisions	
	void ElasticCollide(Particle &ci, Particle &cj, vectorNdim<double,dim> const &L){
  		vectorNdim<double,dim> nij = ci.r-cj.r; nij.in_box(L);
  		nij.normalize(); 
  		vectorNdim<double,dim> vij=ci.v-cj.v;
  		double h=dot(vij,nij);
  		ci.v -=h*nij; ci.ncoll++;
  		cj.v +=h*nij; cj.ncoll++;
	}

	void InelasticCollide(Particle &ci, Particle &cj, vectorNdim<double,dim> const &L){
  		vectorNdim<double,dim> nij = ci.r-cj.r; nij.in_box(L);
  		nij.normalize(); 
  		vectorNdim<double,dim> vij=ci.v-cj.v;
  		double h=(1+eps)*dot(vij,nij)/2;
  		ci.v -=h*nij; ci.ncoll++;
  		cj.v +=h*nij; cj.ncoll++;
	}

	void CollideWhileGrowing(Particle &ci, Particle &cj, vectorNdim<double,dim> const &L){
  		vectorNdim<double,dim> nij = ci.r-cj.r; nij.in_box(L);
  		nij.normalize(); 
  		vectorNdim<double,dim> vij=ci.v-cj.v;
  		double h=dot(vij,nij)+2*growthrate;
  		ci.v -=h*nij; ci.ncoll++;
  		cj.v +=h*nij; cj.ncoll++;
	}

// Update particle(s) up to timenow 
	void Propagate(Particle& ci, vectorNdim<double,dim> const &L){
		double dt=timenow-ci.t; ci.t=timenow;
		ci.r+=ci.v*dt; ci.r.in_box(L);
	}

	void PropagateAndGrow(Particle& ci, vectorNdim<double,dim> const &L){
		double dt=timenow-ci.t; ci.t=timenow;
		ci.r+=ci.v*dt; ci.r.in_box(L); ci.radius+=growthrate*dt; 
	}

	void PropagateAll(){for(int i=0;i<Npart;i++) Propagate(c[i],L);}
	void PropagateAndGrowAll(){for(int i=0;i<Npart;i++) PropagateAndGrow(c[i],L);}

// Check if two particles will collide and put corresponding event on the EventQueue
	void PredictCollision(int i, int j, double tmax){
			double tcoll;
			if((j!=i)&&WillCollide(c[i],c[j],tcoll,tmax,L)) 
				EventQueue.push(Event(tcoll,i,j,c[i].ncoll,c[j].ncoll));
	}

	void PredictCollisionWhileGrowing(int i, int j, double tmax){
			double tcoll;
			if((j!=i)&&WillCollideWhileGrowing(c[i],c[j],tcoll,tmax,L)) 
				EventQueue.push(Event(tcoll,i,j,c[i].ncoll,c[j].ncoll));
	}

// Predict future events for a particle and put them on the EventQueue
	void PredictEvent(int i,double tmax){
		for(int j=0; j<Npart; ++j) PredictCollision(i,j,tmax);
	}

	void PredictEventWhileGrowing(int i,double tmax){
		for(int j=0; j<Npart; ++j) PredictCollisionWhileGrowing(i,j,tmax);
	}


// Growth Process
	void InitGrowthEvents(double dt, double growth=1.0e-20, double tend=1.0e20){
		growthrate=growth;
		endtime=tend; EventQueue.push(Event(endtime,ENDTIME));
		tausync=dt; nexttimesync=timenow+tausync;
		EventQueue.push(Event(nexttimesync,GROWTHRESET));	
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollisionWhileGrowing(i,j,endtime);
	}

	void ResetGrowthEvents(){
		while (!EventQueue.empty()) EventQueue.pop();
		nexttimesync=timenow+tausync; 
		EventQueue.push(Event(nexttimesync,GROWTHRESET));
		PropagateAndGrowAll();
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollisionWhileGrowing(i,j,nexttimesync);	
	}	
	
	Event GrowthStep(){
		assert(!EventQueue.empty());
		Event currentevent=EventQueue.top(); EventQueue.pop();
		int i,j;
		timenow=currentevent.time;
		i=currentevent.i; j=currentevent.j;
		switch(currentevent.kind){
			case PARTICLECOLLISION: 
				if(CollisionIsValid(i,j,currentevent)){
					PropagateAndGrow(c[i],L); 
					PropagateAndGrow(c[j],L);
					CollideWhileGrowing(c[i],c[j],L); 
	   				PredictEventWhileGrowing(i,endtime); 
	   				PredictEventWhileGrowing(j,endtime); 
			} break; 
			case GROWTHRESET: ResetGrowthEvents(); break;
			default: break; // you should not be here 		 
		}
		return currentevent;
	}

// Newtonian Dynamics
	void InitNewtonianEvents(double dt, double tend=1.0e20){
		endtime=tend; EventQueue.push(Event(endtime,ENDTIME));
		tausync=dt; nexttimesync=timenow+tausync;
		EventQueue.push(Event(nexttimesync,NEWTONIANRESET));	
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollision(i,j,endtime);
	}

	void ResetNewtonianEvents(){
		while (!EventQueue.empty()) EventQueue.pop();
		nexttimesync=timenow+tausync; 
		EventQueue.push(Event(nexttimesync,NEWTONIANRESET));
		PropagateAll();
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollision(i,j,nexttimesync);	
	}	
	
	Event NewtonianStep(){
		assert(!EventQueue.empty());
		Event currentevent=EventQueue.top(); EventQueue.pop();
		int i,j;
		timenow=currentevent.time;
		i=currentevent.i; j=currentevent.j;
		switch(currentevent.kind){
			case PARTICLECOLLISION: 
				if(CollisionIsValid(i,j,currentevent)){
					Propagate(c[i],L); Propagate(c[j],L);
					ElasticCollide(c[i],c[j],L); 
	   				PredictEvent(i,endtime); 
	   				PredictEvent(j,endtime); 
			} break; 
			case NEWTONIANRESET: ResetNewtonianEvents(); break;
			default: break; // you should not be here 		 
		}
		return currentevent;
	}


// Brownian Dynamics
	void InitBrownianEvents(double dt,double tend=1.0e20){
		endtime=tend; EventQueue.push(Event(endtime,ENDTIME));	
		tautherm=dt; nexttimetherm=timenow+tautherm;
		EventQueue.push(Event(nexttimetherm,BROWNIANRESET));	
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollision(i,j,nexttimetherm);
	}

	void ResetBrownianEvents(){
		while (!EventQueue.empty()) EventQueue.pop();
		nexttimetherm=timenow+tautherm; 
		EventQueue.push(Event(nexttimetherm,BROWNIANRESET));
		PropagateAll();
		box->ThermalizeVel();
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollision(i,j,nexttimetherm);	
	}	
	
	Event BrownianStep(){
		assert(!EventQueue.empty());
		Event currentevent=EventQueue.top(); EventQueue.pop();
		int i,j;
		timenow=currentevent.time;
		i=currentevent.i; j=currentevent.j;
		switch(currentevent.kind){
			case PARTICLECOLLISION: 
				if(CollisionIsValid(i,j,currentevent)){
					Propagate(c[i],L); Propagate(c[j],L);
					ElasticCollide(c[i],c[j],L); 
	   				PredictEvent(i,nexttimetherm); 
	   				PredictEvent(j,nexttimetherm); 
			} break; 
			case BROWNIANRESET: ResetBrownianEvents(); break;
			default: break; // you should not be here 		 
		}
		return currentevent;
	}

// Periodically Driven Granular Dynamics
	void InitGranularEvents(double dt, double epsilon, double tend=1.0e20){
		endtime=tend; EventQueue.push(Event(endtime,ENDTIME));	
		tautherm=dt; nexttimetherm=timenow+tautherm;
		eps=epsilon;
		EventQueue.push(Event(nexttimetherm,GRANULARRESET));	
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollision(i,j,nexttimetherm);
	}

	void ResetGranularEvents(){
		while (!EventQueue.empty()) EventQueue.pop();
		nexttimetherm=timenow+tautherm; 
		EventQueue.push(Event(nexttimetherm,GRANULARRESET));
		PropagateAll();
		box->ThermalizeVel();
		for(int j=1; j<Npart; ++j)
			for(int i=0; i<j; ++i)
				PredictCollision(i,j,nexttimetherm);	
	}	
	
	Event GranularStep(){
		assert(!EventQueue.empty());
		Event currentevent=EventQueue.top(); EventQueue.pop();
		int i,j;
		timenow=currentevent.time;
		i=currentevent.i; j=currentevent.j;
		switch(currentevent.kind){
			case PARTICLECOLLISION: 
				if(CollisionIsValid(i,j,currentevent)){
					Propagate(c[i],L); Propagate(c[j],L);
					InelasticCollide(c[i],c[j],L); 
	   				PredictEvent(i,nexttimetherm); 
	   				PredictEvent(j,nexttimetherm); 
			} break; 
			case GRANULARRESET: ResetGranularEvents(); break;
			default: break; // you should not be here 		 
		}
		return currentevent;
	}



};
#endif // EDSIMUL_H
