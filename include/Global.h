#ifndef GLOBAL_H
#define GLOBAL_H

/*Floating point precisions*/
float real_se; 		//Surface Energy 		
float real_mobility; // Mobility of the interface
float real_dx; 		// spacing between grid points
float real_totaltime; // Total time for the simulation
float real_dt; //timestep
float real_diffusivity; //Diffusivity
float real_beta; //beta from free energy density

//Non dimensional parameters
/*******************************/
float nd_se;
float nd_mobility; // Mobility of the interface
float nd_dx; 		// spacing between grid points
float nd_totaltime; // Total time for the simulation
float nd_dt; //timestep
float nd_diffusivity;
float nd_beta;
/*****************************/

float beta; //beta for grain force
int beta_timestep;
int phi_timestep;
float alpha; // Dependence of omega on composition
float comp;
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
