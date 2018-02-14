/*! \mainpage HardBrown
 *
 * \section intro_sec Introduction
 * The simulation of hard particles has been longly regarded as a tool for checking idealized theories. Nowadays, Colloidal Science, Granular Matter and Nano Fluidics deal with objects that can be described as hard particles subject to Brownian motion.
 * 
 * While there is a vast choice of simulation tools and libraries for atoms, molecules and in general "soft" particles, this is not the case for "hard" particles. In particular, the standard integration methods for Newtonian and Brownian dynamics are based on the continuity of the inter-particle potential and fail for the steep interactions among hard particles.
 * 
 * Therefore, while simulation methods for systems of many particle have begone a standard in the curriculum of numerical scientists, it is still difficult for them to get an "hand-on" knowledge of hard particles simulations. HardBrown aims to fill such gaps.
 * 
 * For hard particles, it is necessary to resort on Event-Driven simulation methods where the possible collisions among the particles are predicted and then processed through the use of a time-ordered Priority Queue. HardBrown contains programs that can be used as templates for:
 * 
 *  - Event-Driven Newtonian Dynamics
 *  - Event-Driven Granular Dynamics
 *  - Event-Driven Brownian Dynamics
 * 
 * Both for educational and debugging purposes, it is important to access visually to the ongoing simulation. HardBrown has interfaces and template programs for 2d and 3d real-time animated simulations via PlPlot and Inventor libraries.
 * 
 * HardBrown has been developed on an Ubuntu system, and should therefore be easy to install on Ubuntu and Debian based architecture.
 * 
 * The language chosen is C++.
 * 
 * Project can be downloaded at <A HREF="http://download.gna.org/hardshap/"> gna.org </A>
 *
 * \section compile_sec Compilation
 *  
 *  HardBrown is not meant as a library/utility to be installed, but as 
 *  an educational and research tool. It is therefore meant to compile 
 *  and work in a local directory. Everything is written in C++
 *  
 *  The skeletons of Brownian simulation of 2d disks and 3d spheres should 
 *  compile under the requirement that Standard Template Libraries are 
 *  present on your machine.
 *  
 *  A standard Makefile allows the compilation of the simulation 
 *  examples BROWNIAN2d.x and BROWNIAN3d.x:
 *   - make allSIMUL
 *  
 *  For Ubuntu systems where PlPlot libraries are installed, it is 
 * possible to compile real-time 2d animated simulations of Brownian 
 * and Granular dynamics:
 *   - make allPLPLOT
 *   
 *  For Ubuntu systems where Inventor libraries are installed, it is 
 *  possible to compile real-time 3d animated simulations of Newtonian, 
 *  Brownian, Granular dynamics and random growth of bouncing spheres:
 *   - make alINVENTOR
 *  	
 * If both PlPlot and Inventor are present, it is possible to 
 * compile all real-time animated simulations:
 *   - make allGRAPHIC
 *   
 *  
 *  \subsection compile_ubuntu Ubuntu/Debian:
 *  Check that the following libraries are installed: 
 *  - libstlportX.X-dev
 *  (Program have been developed with X.X=5.1, but should not matter)
 *  
 *  For 2d graphic, check for:
 *   - libplplot-dev
 *   - libplplot-c++
 *   - plplot9-driver-xwin
 *  For 3d graphic, check for:
 *   - libinventor
 *   - inventor-dev
 * .
 *  It should be possible to use FreeGlut as an alternative to Inventor
 *  On some Ubuntu systems, it could be also necessary to install
 *   - freeglut3-dev 
 *  
 *  \subsection compile_other other distributions:
 *  For the real-time animations, it is necessary to modify the 
 *  PLPLOT and INVENTOR flags in the Makefile according to the 
 *  hosting system in order to link PlPlot and Inventor libraries.
 * 
 *  \section dirs Files:
 *  - headers
 *  	- gaussdev.h: Random numbers generators
 *  	- vectorNdim.h: Basic template library for vectors in 2d and 3d. 
 *  		      Integer vectors have periodic operators (%, %=)) 
 *  		      useful for the construction of periodic cells. 
 *  	- Disk.h: 2d disks
 *  	- Sphere.h: 3d spheres 
 *  	- SimBox.h: Template simulation boxes for 2d and 3d
 *  	- Event.h: Events for Event Driven simulations
 *  	- EDsimul.h: Core of the project. Template class containing the 
 *                     event queue, collision prediction, and propagation 
 *                     steps for Brownian, Newtonian and Granular simulations  
 *  	- wrap_plplot.h: wrapper for the 2d graphic libraries 
 *  	- InventorShow.h: wrapper for the 3d graphic libraries
 *  	- Gofr.h: calculation of radial distribution function in 2d and 3d
 *  - sources
 *  	- gaussdev.cpp: Random numbers generators
 *  	- Disk.cpp: 2d disks
 *  	- Sphere.cpp: 3d spheres 
 *  	- wrap_plplot.cpp: wrapper for the 2d graphic libraries 
 *  - examples
 *   -# general	
 *  	- BROWNIAN2d.cpp: skeleton of a 2d event driven simulation of
 *                          hard Brownian disks
 *  	- BROWNIAN3d.cpp: skeleton of a 3d event driven simulation of 
 *                          hard Brownian spheres
 *   -# 2d/PlPlot										
 *  	- show2dconf.cpp: shows 2d configuration using PlPlot
 *  	- GRANULARplplot.cpp: Animated event-driven simulation of 
 *                              2d shaken granular disks
 *  	- BROWNIANplplot.cpp: Animated event-driven simulation of 
 *                              hard Brownian disks
 *   -# 3d/OpenInventor						
 *  	- show3dconf.cpp: shows 3d configuration using Open Inventor
 *  	- BROWNIANinventor.cpp: Animated event-driven simulation of 
 *                                hard Brownian spheres 
 *  	- NEWTONIANinventor.cpp: Animated event-driven simulation of 
 *                                 hard spheres following Newtonian dynamics 
 *  	- GROWTHinventor.cpp: Animation of the growth of a low density
 *                              configuration of spheres up to a final
 *                              packing fraction 
 *    .
 *  
 */

