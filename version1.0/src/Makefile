#CXXFLAGS = -O6
CXXFLAGS = -g -Wall
CXX = g++ $(CXXFLAGS)

CXXINCLUDE = -I./

# Define here the flags to compile with the PlPlot libraries 
PLPLOTINC = `pkg-config --cflags plplotd-c++`
PLPLOT = $(PLPLOTINC) `pkg-config --libs plplotd-c++`
#PLPLOT = -I/usr/include/plplot  -lplplotcxxd -lplplotd -lltdl -ldl -lm -lcsirocsa -lcsironn -lqhull -lfreetype  

# Define here the flags to compile with the Inventor libraries 
INVENTOR = -lInventor -lInventorXt

OBJS = gaussdev.cpp Disk.cpp Sphere.cpp 

default:
	@echo
	@echo "*** HARD particles BROWNian and granular simulation ***"
	@echo
	@echo "		A general skeleton for 2d/3d simulations is given"
	@echo "		in BROWNIAN2d.cpp and BROWNIAN3d.cpp "
	@echo
	@echo "		If PlPlot libraries are available on your system, "
	@echo "		it is possible to animate 2d simulation of Brownian"
	@echo "		and Granular particles (see BROWNIANplplot.cpp ,"
	@echo "		GRANULARplplot.cpp)"
	@echo
	@echo "		If Inventor libraries are available on your system, "
	@echo "		it is possible to animate 3d simulation of Brownian"
	@echo "		and Newtonian particles (see BROWNIANinventor.cpp ,"
	@echo "		NEWTONIANinventor.cpp); also, it is possible to"
	@echo "		simulate growing particles (GROWTHinventor.cpp)"
	@echo 
	@echo "				options:"
	@echo "make all"
	@echo "make allSIMUL"
	@echo "make allGRAPHIC"
	@echo "make allPLPLOT"
	@echo "make allINVENTOR"
	@echo

all: allSIMUL allGRAPHIC

allSIMUL: BROWNIAN2d BROWNIAN3d
allGRAPHIC: allPLPLOT allINVENTOR  
allINVENTOR: GROWTHinventor BROWNIANinventor NEWTONIANinventor show3dconf
allPLPLOT: BROWNIANplplot GRANULARplplot show2dconf
	
show2dconf: wrap_plplot
	$(CXX) $@.cpp $(OBJS) wrap_plplot.o -o $@.x $(CXXINCLUDE) $(PLPLOT) 

show3dconf: 
	$(CXX) $@.cpp $(OBJS) -o $@.x $(CXXINCLUDE) $(INVENTOR)

NEWTONIANinventor: 
	$(CXX) $@.cpp $(OBJS) -o $@.x $(CXXINCLUDE) $(INVENTOR)

GROWTHinventor: 
	$(CXX) $@.cpp $(OBJS) -o $@.x $(CXXINCLUDE) $(INVENTOR)

BROWNIANinventor: 
	$(CXX) $@.cpp $(OBJS) -o $@.x $(CXXINCLUDE) $(INVENTOR)

BROWNIANplplot: wrap_plplot
	$(CXX) $@.cpp $(OBJS) wrap_plplot.o -o $@.x $(CXXINCLUDE) $(PLPLOT) 

GRANULARplplot: wrap_plplot
	$(CXX) $@.cpp $(OBJS) wrap_plplot.o -o $@.x $(CXXINCLUDE) $(PLPLOT) 


BROWNIAN3d: 
	$(CXX) $@.cpp $(OBJS) -o $@.x $(CXXINCLUDE)

BROWNIAN2d: wrap_plplot
	$(CXX) $@.cpp $(OBJS) -o $@.x $(CXXINCLUDE) $(PLPLOT) 

wrap_plplot:
	$(CXX) -c $@.cpp $(CXXINCLUDE) $(PLPLOTINC) 
	
clean:
	rm -f *.o *~ core *.out *.x *.cnf
