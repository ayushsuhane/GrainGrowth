#ifndef GLOBAL_H
#define GLOBAL_H

/*Floating point precisions*/
float se; 		//Surface Energy 		
float mobility; // Mobility of the interface
float dx; 		// spacing between grid points
float total_time; // Total time for the simulation
float dt; //timestep
float beta; //beta for grain force
float alpha; // Dependence of omega on composition
float comp;
float diffusivity;
float temperature;

float epsilon; // Free energy of the form =  `sum_{a,b} (epsilon^2/2)(\nabla \phi_a)(\nabla \phi_b) 
float omega;	//  + omega(\phi_a)(\phi_b)(1 - \alpha*c)
float mob_phi; // parametrized mobility

/*Integer points*/
int gridx; 	// Gridpoints in x 
int gridy;  // Gridpoints in y
int gridsize;
int eta;		// Interface thickness in grid points

char *bc;


int grainforce;

#define PI 3.14159
#endif
