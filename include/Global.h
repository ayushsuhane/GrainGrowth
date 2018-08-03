#ifndef GLOBAL_H
#define GLOBAL_H


/*doubleing point precisions*/
double real_se; 		//Surface Energy 		
double real_mobility; // Mobility of the interface
double real_dx; 		// spacing between grid points
double real_totaltime; // Total time for the simulation
double real_dt; //timestep
double real_diffusivity; //Diffusivity
double real_beta; //beta from free energy density

//Non dimensional parameters
/*******************************/
double c_time, c_length, c_diffusivity, c_mobility, c_energy;
double nd_se;
double nd_mobility; // Mobility of the interface
double nd_dx; 		// spacing between grid points
double nd_totaltime; // Total time for the simulation
double nd_dt; //timestep
double nd_diffusivity;
double nd_beta;
/*****************************/
double dt;
double beta; //beta for grain force
int beta_timestep;
int phi_timestep;
double alpha; // Dependence of omega on composition
double initcomp;
double temperature;

double epsilon; // Free energy of the form =  `sum_{a,b} (epsilon^2/2)(\nabla \phi_a)(\nabla \phi_b) 
double omega;	//  + omega(\phi_a)(\phi_b)(1 - \alpha*c)
double mob_phi; // parametrized mobility

/*Integer points*/
int gridx; 	// Gridpoints in x 
int gridy;  // Gridpoints in y
int gridsize;
int eta;		// Interface thickness in grid points

char *bc;
char inputfilename[100], controlinputfilename[100];
int noutput;

char init[100], geometry[100],  phisolver[100], csolver[100], initialfile[100];

int grainforce;
int flag_readfromfile, flag_writetofile;
double time_constant;

#define PI 3.14159
#endif
