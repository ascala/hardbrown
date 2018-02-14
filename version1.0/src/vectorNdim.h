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

/*! \file vectorNdim.h
 * \brief Basic template library for vectors in 2d and 3d.
 *  
 * Integer vectors have periodic operators (%, %=)) 
 * useful for the construction of periodic cells. 
 * */ 

#ifdef _HAVE_PETE_
#include <vectorNdimPETE.h>
#else

#ifndef VECTORNDIM_H_
#define VECTORNDIM_H_

#include <iostream>

#include<cctype>
#include<climits>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cmath>

#include "gaussdev.h"

template<class T, int dim> class vectorNdim;
//typedef vectorNdim<double,2> Vec2d ;
//typedef vectorNdim<int,2> IVec2d ;
//typedef vectorNdim<double,3> Vec3d ;
//typedef vectorNdim<int,3> IVec3d ;

//!  2d vectors 
#define Vec2d vectorNdim<double,2> 
//!  3d vectors 
#define Vec3d vectorNdim<double,3> 
//!  integer 2d vectors 
#define IVec2d vectorNdim<int,2> 
//!  integer 3d vectors 
#define IVec3d vectorNdim<int,3> 

#define FORDIM for(int i=0; i<dim;++i)

//! Template class for vectors
/*!Overloads standard operators (+,-,....) 
 * */
template <class T, int dim>
class vectorNdim {
public:
    T vec[dim];

  	T &operator[](int i)      { return vec[i]; }
  	T operator[](int i) const { return vec[i]; }

    vectorNdim(){FORDIM vec[i]=0;}
    vectorNdim(T const& v){FORDIM vec[i]=v;}
    vectorNdim(T const& vx, T const& vy){vec[0]=vx; vec[1]=vy;}
    vectorNdim(T const& vx, T const& vy, T const& vz){vec[0]=vx; vec[1]=vy; vec[2]=vz;}
    vectorNdim(vectorNdim const& vv){FORDIM vec[i]=vv.vec[i];}
    virtual ~vectorNdim() {}

    //	stream operators
    template<class U > friend std::ostream& 
    operator<<(std::ostream& pStr, const vectorNdim<U,dim>& pV)
    { FORDIM pStr <<  pV.vec[i] <<" "; return pStr; }

    template<class U> friend std::istream& 
    operator>>(std::istream& pStr, vectorNdim<U,dim>& pV)
	{ FORDIM pStr >> pV.vec[i]; return pStr; }


    // vector product/divisions (* is not the dot product !!!)
    vectorNdim operator*(vectorNdim const& vv) const {
        vectorNdim out;
        FORDIM out[i]=vec[i]*vv.vec[i];
        return out;
    }
    vectorNdim operator/(vectorNdim const& vv) const {
        vectorNdim out;
        FORDIM out[i]=vec[i]/vv.vec[i];
        return out;
    }
    vectorNdim& operator*=(vectorNdim const& vv)
    { FORDIM vec[i]*=vv.vec[i]; return *this; }

    vectorNdim& operator/=(vectorNdim const& vv)
    { FORDIM vec[i]/=vv.vec[i]; return *this; }

    vectorNdim operator%(vectorNdim<int,dim> const& vv) const { 
    	vectorNdim<int,dim> out;
    	FORDIM out[i]=vec[i]%vv.vec[i]; 
    	return out; 
    }
    
    vectorNdim& operator%=(vectorNdim<int,dim> const& vv)
    { FORDIM vec[i]%=vv.vec[i]; return *this; }

    // mapping to a lattice
    int index(vectorNdim<int,dim> const &L) { 
    	int prod=1,ret=0; 
    	FORDIM{ ret+=prod*vec[i];  prod*=L[i];}
    	return ret; 
    }
    vectorNdim from_index(int ind, vectorNdim<int,dim> const &L) { 
    	vectorNdim vv; 
    	FORDIM {vv.vec[i] = ind % L[i]; ind /= L[i];} 
        return vv;
    }


    // assignment (overloading =)
    vectorNdim &operator=(vectorNdim const& other) {
        if (this != &other) FORDIM vec[i] = other.vec[i];
        return *this;
    }

    // vector sum
    vectorNdim &operator+=(vectorNdim const& other) 
    { FORDIM vec[i]+=other.vec[i]; return *this; }

    vectorNdim const operator+(void) { return *this; }

    vectorNdim const operator+(vectorNdim const& rvalue) 
    { vectorNdim out; FORDIM out[i]=vec[i]+rvalue[i]; return out;}

    template <class U>
    friend vectorNdim<U,dim> const 
    operator+(vectorNdim<U,dim> const &l_hand, vectorNdim<U,dim> const &r_hand)
	{return vectorNdim<T,dim>(l_hand) += r_hand;}
	
    // vector difference
    vectorNdim &operator-=(vectorNdim const& other) 
    { FORDIM vec[i]-=other.vec[i]; return *this; }

    vectorNdim const operator-(void)
    { vectorNdim out; FORDIM out[i]=vec[i]; return out;}

    vectorNdim const operator-(vectorNdim const& rvalue) 
    { vectorNdim out; FORDIM out[i]=vec[i]-rvalue[i]; return out;}

    template <class U>
    friend vectorNdim<U, dim> const 
    operator-(vectorNdim<U, dim> const& l_hand, vectorNdim<U, dim> const& r_hand)
	{ return vectorNdim<T,dim>(l_hand) -= r_hand;}
	
    // dot product (NO * OVERLOADING) and norm
    T dot(const vectorNdim& rvalue)
    {T ret=0; FORDIM ret+=vec[i]*rvalue[i]; return ret;}
    friend T dot(const vectorNdim& a,const vectorNdim& b)
    {T ret=0; FORDIM ret+=a[i]*b[i]; return ret;}

	// cross product
    T cross(const vectorNdim<T,2>& rvalue)
    { return vec[0]*rvalue[1]-vec[1]*rvalue[0]; }
    template <class U> friend U cross(const vectorNdim<U,2>& a,const vectorNdim<U,2>& b);
    vectorNdim<T,3> cross(const vectorNdim<T,3>& rvalue){
    	vectorNdim<T,3> out; 
    	out[0]=vec[1]*rvalue[2]-vec[2]*rvalue[1]; 
    	out[1]=vec[2]*rvalue[0]-vec[0]*rvalue[2]; 
    	out[2]=vec[0]*rvalue[1]-vec[1]*rvalue[0]; 
    	return out;
    }
    template <class U> friend vectorNdim<U,3> cross(const vectorNdim<U,3>& a,const vectorNdim<U,3>& b);

    // normalization
    T norm(void) 
    { T ret=0; FORDIM ret+=vec[i]*vec[i]; return sqrt(ret); }
    T norm2(void) 
    { T ret=0; FORDIM ret+=vec[i]*vec[i]; return ret; }
    T normalize(void) 
    { T ret=0; FORDIM ret+=vec[i]*vec[i]; ret=sqrt(ret); FORDIM vec[i]/=ret; return ret; }

    friend T norm(vectorNdim const &v) 
    { T ret=0; FORDIM ret+=v[i]*v[i]; return sqrt(ret); }
    friend T norm2(vectorNdim const &v) 
    { T ret=0; FORDIM ret+=v[i]*v[i]; return ret; }
    friend T normalize(vectorNdim &v) 
    { T ret=0; FORDIM ret+=v[i]*v[i]; ret=sqrt(ret); FORDIM v[i]/=ret; return ret; }

    // multiplication by a scalar
    template <class U> friend vectorNdim<T,dim> const 
    operator*(U const &l_hand, vectorNdim<T,dim> const& r_hand)
    { vectorNdim<T,dim> out; FORDIM out[i]=l_hand*r_hand.vec[i]; return out;}
    
    template <class U> friend vectorNdim<T,dim> const 
    operator*(vectorNdim<T,dim> const& l_hand, U const &r_hand)
    { vectorNdim<T,dim> out; FORDIM out[i]=l_hand[i]*r_hand; return out;}
    
    vectorNdim const operator*=(T const &scalar)
    { FORDIM vec[i]*=scalar; return *this; }
	
    // division by a scalar
    template <class U> friend vectorNdim<U,dim> const 
    operator/(vectorNdim<U,dim> const& l_hand, U const &r_hand)
    { vectorNdim<T,dim> out; FORDIM out[i]=l_hand[i]/r_hand; return out;}
    
    vectorNdim const operator/=(T const &scalar)
    { FORDIM vec[i]/=scalar; return *this; }
    
    // equality & inequality
    bool const operator==(vectorNdim const &rvalue) {
        FORDIM if (vec[i]!=rvalue[i]) return false;
        return true;
    }
    template <class U> friend bool const 
    operator==(vectorNdim<U,dim> const& l_hand, vectorNdim<U,dim> const& r_hand) 
    { return vectorNdim<T,dim>(l_hand)==r_hand;}
    
    bool const operator!=(vectorNdim const &rvalue) {
        FORDIM if (vec[i]!=rvalue[i]) return true;
        return false;
    }
    
    template <class U> friend bool const 
    operator!=(vectorNdim<U,dim> const& l_hand, vectorNdim <U,dim>const& r_hand) 
    { return vectorNdim<T,dim>(l_hand)!=r_hand;}
    
    // "area"
    T const area(){if(dim==2) return vec[0]*vec[1]; else return -1;}
    T const volume(){if(dim==3) return vec[0]*vec[1]*vec[2]; else return -1;}
    T const d_volume(){double dvol=vec[0]; for(int i=1; i<dim;++i) dvol*=vec[i]; return dvol;}
    
    // min,max
    T const min(){double m=vec[0]; for(int i=1; i<dim;++i) if(vec[i]<m) m=vec[i]; return m;}
    T const max(){double m=vec[0]; for(int i=1; i<dim;++i) if(vec[i]>m) m=vec[i]; return m;}


    // 2d rotations
 /*   void rotate(T const angle) {
        T sa=sin(angle), ca=cos(angle);
        T xr(ca*vec[0] + sa*vec[1]), yr(-sa*vec[0] + ca*vec[1]);
        vec[0]=xr;
        vec[1]=yr;
    }
*/

    // random vectors
    void rnd(void) { FORDIM vec[i]=drand(); }
    void rnd(T const s) { FORDIM vec[i]=s*(2*drand48()-1); }
    void rnd(vectorNdim const &s) { FORDIM vec[i]=s.vec[i]*(2*drand48()-1); }
    void rnd(T const sx, T const sy) { 
    	vec[0]=sx*(2*drand48()-1);
        vec[1]=sy*(2*drand48()-1);
    }
    void rnd(T const sx, T const sy, T const sz) { 
    	vec[0]=sx*(2*drand48()-1);
        vec[1]=sy*(2*drand48()-1);
        vec[2]=sz*(2*drand48()-1);
    }

    void gauss(void) { FORDIM vec[i]=gaussdev(); }
    void gauss(T const W) { FORDIM vec[i]=W*gaussdev(); }
    void gauss(vectorNdim const &W) { FORDIM vec[0]=W.vec[0]*gaussdev();}
    void gauss(T const sx, T const sy) { 
    	vec[0]=sx*(2*gaussdev()-1);
        vec[1]=sy*(2*gaussdev()-1);
    }
    void gauss(T const sx, T const sy, T const sz) { 
    	vec[0]=sx*(2*gaussdev()-1);
        vec[1]=sy*(2*gaussdev()-1);
        vec[2]=sz*(2*gaussdev()-1);
    }

    // random displacements
    void rnd_displ(void) { FORDIM vec[i]+=drand48(); }
    void rnd_displ(T const s) { FORDIM vec[i]+=s*(2*drand48()-1); }
    void rnd_displ(vectorNdim const &s) { FORDIM vec[i]=s.vec[i]*(2*drand48()-1); }
    void rnd_displ(T const sx, T const sy) {
        vec[0]+=sx*(2*drand48()-1);
        vec[1]+=sy*(2*drand48()-1);
    }
    void rnd_displ(T const sx, T const sy, T const sz) {
        vec[0]+=sx*(2*drand48()-1);
        vec[1]+=sy*(2*drand48()-1);
        vec[2]+=sy*(2*drand48()-1);
    }

    void gauss_displ(void) { FORDIM vec[i]+=gaussdev(); }
    void gauss_displ(T const W) { FORDIM vec[i]+=W*gaussdev(); }
    void gauss_displ(vectorNdim const &W) { FORDIM vec[0]+=W.vec[0]*gaussdev();}
    void gauss_displ(T const sx, T const sy) { 
    	vec[0]+=sx*(2*gaussdev()-1);
        vec[1]+=sy*(2*gaussdev()-1);
    }
    void gauss_displ(T const sx, T const sy, T const sz) { 
    	vec[0]+=sx*(2*gaussdev()-1);
        vec[1]+=sy*(2*gaussdev()-1);
        vec[2]+=sz*(2*gaussdev()-1);
    }

    // periodic boundaries
    void in_box(void) { FORDIM vec[i]-=rint(vec[i]); }
    void in_box(T const L) { FORDIM vec[i]-=L*rint(vec[i]/L); }
    void in_box(T const Lx, T const Ly) 
    { vec[0]-=Lx*rint(vec[0]/Lx); 
      vec[1]-=Ly*rint(vec[1]/Ly); }
    void in_box(T const Lx, T const Ly, T const Lz) 
    { vec[0]-=Lx*rint(vec[0]/Lx); 
      vec[1]-=Ly*rint(vec[1]/Ly); 
      vec[2]-=Lz*rint(vec[2]/Lz); }
    void in_box(vectorNdim const L) 
    { FORDIM vec[i]-=L.vec[i]*rint(vec[i]/L.vec[i]); }
    
};

template <class T> T cross(const vectorNdim<T,2>& a,const vectorNdim<T,2>& b) 
{ return a.vec[0]*b.vec[1]-a.vec[1]*b.vec[0]; }

template <class T>  vectorNdim<T,3> cross(const vectorNdim<T,3>& a,const vectorNdim<T,3>& b){
    vectorNdim<T,3> out; 
    out[0]=a[1]*b[2]-a[2]*b[1]; 
    out[1]=a[2]*b[0]-a[0]*b[2]; 
    out[2]=a[0]*b[1]-a[1]*b[0]; 
    return out;
}

#undef FORDIM

#endif /*VECTORNDIM_H_*/


#endif // _HAVE_PETE_
