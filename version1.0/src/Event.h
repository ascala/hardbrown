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

/*! \file Event.h
 * \brief Events for Event-Driven simulations
 * */

#ifndef EVENT_H_
#define EVENT_H_

#include <iostream>

//! possible events
enum EventKind { 
	PARTICLECOLLISION, 
	BROWNIANRESET,
	NEWTONIANRESET,
	GROWTHRESET,
	GRANULARRESET,
	ENDTIME,
	NOACTION
};

//! names of possible events
const char *KindName[]={ 
	"PARTICLECOLLISION", 
	"BROWNIANRESET",
	"NEWTONIANRESET",
	"GROWTHRESET",
	"GRANULARRESET",
	"ENDTIME",
	"NOACTION"
};

//! Class for the scheduled events
/*!Events are scheduled to occurr at time "time"; the can involve 
 * one particle, two particles or be general collective events
 * 
 * If the event involves a single particle, "i" is the index of 
 * the particle and "nci" is the number of collision suffered by 
 * the i-th particle at the scheduling time
 * 
 * If the event involves two particles, their indices are "i" and
 * "j" and "nci", "ncj" are the number of collision suffered by the
 * involved particles at the scheduling time
 * */
class Event{
public:		
		
        double time; //!< time for which the event is scheduled
        
        
        int i; //!< index of first particle involved in an event
        int j; //!< index of second particle involved in an event
        
        long int nci; //!< # of collisions of first particle
        long int ncj; //!< # of collisions of second particle

        EventKind kind; //!< Event kind
        
		// constructor/destructor
        Event(){};
        Event(Event const& e)
                :time(e.time),i(e.i),j(e.j),nci(e.nci),ncj(e.ncj),kind(e.kind){};
        Event(double tt, int ii, int jj, int ci, int cj, EventKind kk=PARTICLECOLLISION)
                :time(tt),i(ii),j(jj),nci(ci),ncj(cj),kind(kk){};
        Event(double tt, EventKind kk=NOACTION):time(tt),kind(kk){};
        virtual ~Event(){};
        
        // "<" operator for time ordering
        bool operator<(Event const &right) const { return time > right.time;}
    
    	// formatted output
    	void print(std::ostream &os) const { 
    		os<<KindName[kind]<<" "<<time;
    		if(kind==PARTICLECOLLISION)
    			os<<" "<<i<<" "<<j<<" "<<nci<<" "<<ncj;
    	}

};

// overloading of <<
std::ostream &operator<<(std::ostream &os, const Event &e)
        { e.print(os); return os;}


#endif // EVENT_H_
