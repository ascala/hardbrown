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

/*! \file gaussdev.h
 * \brief gaussian random number generator
 * 
 * uses Box-Muller method
 * */

#ifndef GAUSSDEV_H
#define GAUSSDEV_H

double drand(void);
double gaussdev(void);

double gsum(void);
double gsum48(void);
double grand(void);
double grand48(void);

#endif //GAUSSDEV_H
