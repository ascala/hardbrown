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

/*! \file InventorShow.h
 * \brief wrapper for Open Inventor
 * */


#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>

#include <Sphere.h>
#include <SimBox.h>

#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoGroup.h>

#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoDrawStyle.h>

#include <Inventor/sensors/SoTimerSensor.h>

SoMaterial *colorTransparent(double r, double g, double b){ 
	SoMaterial *m= new SoMaterial;
	m->ambientColor.setValue(r, g, b);
	m->diffuseColor.setValue(r, g, b);
	m->specularColor.setValue(0.5, 0.5, 0.5);
	m->transparency = 0.95;
	m->shininess = 0.5;
	return m;
}
SoMaterial *colorPlastic(double r, double g, double b){ 
	SoMaterial *m= new SoMaterial;
	m->ambientColor.setValue(r, g, b);
	m->diffuseColor.setValue(r, g, b);
	m->specularColor.setValue(0.5, 0.5, 0.5);
	m->shininess = 0.5;
	return m;
}
SoMaterial *whitePlastic(){ return colorPlastic(1,1,1); }
SoMaterial *redPlastic(){ return colorPlastic(1,0,0); }
SoMaterial *greenPlastic(){ return colorPlastic(0,1,0); }
SoMaterial *bluePlastic(){ return colorPlastic(0,0,1); }
SoMaterial *cyanPlastic(){ return colorPlastic(0,1,1); }
SoMaterial *magentaPlastic(){ return colorPlastic(1,0,1); }
SoMaterial *yellowPlastic(){ return colorPlastic(1,1,0); }
SoMaterial *bronze(){ 
	SoMaterial *m= new SoMaterial;
	m->ambientColor.setValue(.33, .22, .27);
	m->diffuseColor.setValue(.78, .57, .11);
	m->specularColor.setValue(.99, .94, .81);
	m->shininess = .28;
	return m;
}
SoMaterial *silver(){ 
	SoMaterial *m= new SoMaterial;
	m->ambientColor.setValue(.2, .2, .2);
	m->diffuseColor.setValue(.6, .6, .6);
	m->specularColor.setValue(.5, .5, .5);
	m->shininess = .5;
	return m;
}


SoSeparator *MySoSimbox(){
  SoSeparator *root =  new SoSeparator; root->ref();   
  
  root->addChild(colorTransparent(1,1,0));
  
  SoDrawStyle *ds = new SoDrawStyle;
  ds->lineWidth=4;
  ds->style=SoDrawStyle::LINES;
  root->addChild(ds);
  SoTransform *myXform = new SoTransform; 
  myXform->scaleFactor.setValue(box.L[0]/2,box.L[1]/2,box.L[2]/2);
  root->addChild(myXform); 
  root->addChild(new SoCube); 
  
  root->unrefNoDelete(); return root;
}

const int Nmax=10000;
SoTransform *position[Nmax]; 


void AddSoParticles(SoSeparator *root){
//SoTransform *position; 
  for(int i=0;i<box.Npart;i++){
	position[i] = new SoTransform;
	double radius=box.particles[i].radius;
	position[i]->scaleFactor.setValue(radius, radius, radius);
	position[i]->translation.setValue(
		box.particles[i].r[0],
		box.particles[i].r[1],
		box.particles[i].r[2]);

  	SoSeparator *particle = new SoSeparator;
	particle->addChild(position[i]); 
//	particle->addChild(colorTransparent(0,1,0));
  	particle->addChild(bronze());
	particle->addChild(new SoSphere); 

   	root->addChild(particle);

   }

}

void InventorShow()
{

//	   SoDB::init();
	
	
   Widget myWindow = SoXt::init("3d spheres");
   if (myWindow == NULL) exit(1);

   SoSeparator *root = new SoSeparator;
   root->ref();
   
   SoMaterial *myMaterial = new SoMaterial;
   myMaterial->diffuseColor.setValue(1.0, 0.0, 0.0);
   root->addChild(myMaterial);

 
// simulation box
  root->addChild(MySoSimbox()); 

//particles 
  AddSoParticles(root); 
	
	// Render Area
	SoXtRenderArea *myRenderArea = new SoXtRenderArea(myWindow);
    myRenderArea->setSceneGraph(root);
//    myRenderArea->setTitle("2d");
//    myRenderArea->show();
 
   // Set up viewer:
   SoXtExaminerViewer *myViewer = 
            new SoXtExaminerViewer(myWindow);
   myViewer->setSceneGraph(root);
   myViewer->setTitle("Examine 3d Spheres");
   myViewer->show();
	
	
   
	SoXt::show(myWindow);
   	SoXt::mainLoop();


}
	
void simulstep();
void animate(void *data, SoSensor *) {
  simulstep();	
  for(int i=0;i<box.Npart;i++){ 
  	double radius=box.particles[i].radius;
	position[i]->scaleFactor.setValue(radius, radius, radius);
  	position[i]->translation.setValue(
		box.particles[i].r[0],
		box.particles[i].r[1],
		box.particles[i].r[2]);
  }
}


void InventorAnimate(float timer){

//	   SoDB::init();
	
	
   Widget myWindow = SoXt::init("EDBD3d");
   if (myWindow == NULL) exit(1);

   SoSeparator *root = new SoSeparator;
   root->ref();
   
   SoMaterial *myMaterial = new SoMaterial;
   myMaterial->diffuseColor.setValue(1.0, 0.0, 0.0);
   root->addChild(myMaterial);

 
// simulation box
  root->addChild(MySoSimbox()); 

//particles 
  AddSoParticles(root); 
	
	// Render Area
	SoXtRenderArea *myRenderArea = new SoXtRenderArea(myWindow);
    myRenderArea->setSceneGraph(root);
//    myRenderArea->setTitle("2d");
//    myRenderArea->show();

 

SoTimerSensor *sensor = new SoTimerSensor(animate, root);
sensor->setInterval(SbTime(timer)); // 1.0/timer frames per second
sensor->schedule();
 
   // Set up viewer:
   SoXtExaminerViewer *myViewer = 
            new SoXtExaminerViewer(myWindow);
   myViewer->setSceneGraph(root);
   myViewer->setTitle("Examine 3d EDBD");
   myViewer->show();
	
	
   
	SoXt::show(myWindow);
   	SoXt::mainLoop();


}
