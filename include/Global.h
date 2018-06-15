#ifndef GLOBAL_H
#define GLOBAL_H

/*Floating point precisions*/
float se; 		//Surface Energy 		
float mobility; // Mobility of the interface
float dx; 		// spacing between grid points
float total_time; // Total time for the simulation
float dt; //timestep
float beta; //beta for grain force

/*Integer points*/
int gridsize; 	// Gridsize in points
int eta;		// Interface thickness in grid points


int grainforce;

#define PI 3.14159
#endif
