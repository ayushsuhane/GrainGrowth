#include "Functions.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "Global.h"
#include "nodefunctions.h"
#include<math.h>

int Initialize(char *init, char *geom)
{
	printf("Initialization : %s\n", init);
	printf("Initializing geometry : %s\n", geom);
	/*
	Default initialization for two grains in a 1D subsystem
	*/
	
	//Read the input parameters from file
	Readinputs();
	
	//Print Dimensionalization consants for scaling back to real values
	NonDimensionalize();

	// Evaluate the parameters used in the phase field simulations
	Parametrize(init);

	//Print the settings on the screen
	Printsettings();
	
	//Allocate the relevant memory 
	Allocate();
	
	//Setup the geometry
	Setup(geom);

	//randomTest();
	//Check_initialstate();
	return 1;

}

void Readinputs(char *init)
{
	/* 
	Extern Variable declarations
	*/
	extern float real_mobility;
	extern float real_se;
	extern float real_dx;
	extern float real_totaltime;
	extern float real_dt;
	extern float real_diffusivity;
	extern float real_beta; 	// Artificial force coefficient
	extern int beta_timestep;		//timestep after which beta gets activated
	extern int phi_timestep; 		//timestep to form the interfaces
	
	extern float alpha; 		// Parameter representing interaction between solute atoms and Grain Boundary
	extern float comp;			// Initial Composition
	extern float temperature;


	extern int gridx, gridy;
	extern int gridsize;
	extern int eta;				// Number of grid points for the interface
	extern int grainforce;		// Grain in which additional artificial force needs to be activated

	extern char *bc;			// String for boundary condition

	bc = (char *)malloc(15*sizeof(char));
	
	grainforce = 0;
	real_beta = 0.0;
	alpha = 0.0; //If not present, then no segregation potential
	real_diffusivity = 0; //If not specified
	comp = 0.0; 
	temperature = 1273.0; //Initialize at 273 K
	real_dt = 0.01; //Just in case
	phi_timestep = 0;
	/******/

	char FILENAME[] = "../input/parameters.txt";
	FILE *fp;
	fp = fopen(FILENAME, "r");
	
	//Error check
	if(fp == NULL) 
	{
		printf("File to read parameters not  found");
		return;
	}
	printf("%s : Reading Parameters\n", FILENAME);

	char tmpstr[15], tmpval[15];
	char tempbuffer[100];
	while(!feof(fp))
	{
		if(fgets(tempbuffer, 100, fp))
		{
			sscanf(tempbuffer, "%15s : %15[^;]", tmpstr, tmpval);
			tmpstr[strcspn(tmpstr, "\r\n")] = 0;
			tmpval[strcspn(tmpval, "\r\n")] = 0;

			if(!strcmp(tmpstr, "Mobility")) real_mobility = atof(tmpval);
			if(!strcmp(tmpstr, "SurfaceEnergy")) real_se = atof(tmpval);
			if(!strcmp(tmpstr, "dx")) real_dx = atof(tmpval);
			if(!strcmp(tmpstr, "GridX")) gridx = atoi(tmpval);
			if(!strcmp(tmpstr, "GridY")) gridy = atoi(tmpval);
			if(!strcmp(tmpstr, "Eta")) eta = atoi(tmpval);
			if(!strcmp(tmpstr, "TotalTime")) real_totaltime = atof(tmpval);
			if(!strcmp(tmpstr, "dt")) real_dt = atof(tmpval);
			if(!strcmp(tmpstr, "GrainForce")) grainforce = atoi(tmpval);
			if(!strcmp(tmpstr, "Beta")) real_beta = atof(tmpval);
			if(!strcmp(tmpstr, "Betatimestep")) beta_timestep = atoi(tmpval);
			if(!strcmp(tmpstr, "Phitimestep")) phi_timestep = atoi(tmpval);
			if(!strcmp(tmpstr, "Alpha")) alpha = atof(tmpval);
			if(!strcmp(tmpstr, "Diffusivity")) real_diffusivity = atof(tmpval);
			if(!strcmp(tmpstr, "Composition")) comp = atof(tmpval);
			if(!strcmp(tmpstr, "Temperature")) temperature = atof(tmpval);
			if(!strcmp(tmpstr, "BC")) strcpy(bc,tmpval);
		}
	}
	fclose(fp);
	gridsize = gridx*gridy;
	return;
}

void NonDimensionalize()
{
	/*
	Contains the non-dimensional parameters. Check the input with these parameters first
	*/
	extern float temperature;
	float gas_constant = 8.314; // J/K
	float molar_volume = 7.0e-6; // 1/m^3
	/*Non dimensionalization Constants*/
	float c_length, c_energy, c_mobility, c_time;
	float c_diffusivity;
	float c_se;
	/********************************/
	c_diffusivity = 1e-14; // m^2/s
	c_length = 1e-9;     // m
	c_energy = (gas_constant*temperature/molar_volume) ; // J/m^2
	c_se = c_energy*c_length;
	c_time = (c_length*c_length)/c_diffusivity;
	c_mobility = c_length/(c_time*c_energy);      //  m^4/(Js)
	
	printf("Non Dimensional Constants: \n");
	printf("For more information : check Dimensionalize function\n");
	/*
	These constants are defined to scale the quantities back to the real values 
	such that 
	A = a* X a
	where the values here are a*, the values in parameters.txt are a and actual values can be obtained 
	by multiplying these values depending on their respective scaling behaviour.
	*/
	printf("Length Scale : %0.2e\n", c_length);
	printf("Time Scale : %0.2e\n", c_time);
	printf("Energy Scale : %0.2e\n", c_energy);
	printf("Surface Energy Scale Constant : %0.2e\n", c_se);
	printf("Mobility scaling constant : %0.2e\n", c_mobility);
	printf("Diffusion scaling constant : %0.2e\n", c_diffusivity);

	extern float nd_mobility, real_mobility;
	extern float nd_se, real_se;
	extern float nd_dx, real_dx;
	extern float nd_dt, real_dt;
	extern float nd_totaltime, real_totaltime;
	extern float nd_diffusivity, real_diffusivity;
	extern float nd_beta, real_beta;

	nd_mobility =  real_mobility/c_mobility;
	nd_se = real_se/c_se;
	nd_dx = real_dx/c_length;
	//nd_dt = real_dt/c_time;
	nd_totaltime = real_totaltime/c_time;
	nd_diffusivity = real_diffusivity/c_diffusivity;
	nd_dt = nd_dx*nd_dx/(100.0*nd_diffusivity);
	nd_beta = real_beta/c_energy;

	return;
}

void Parametrize(char *init)
{
	// Define the parameters like epsilon and omega, mobility_phi here, which 
	// can be used in solver functions
	extern float nd_se;
	extern float nd_mobility;
	extern int eta; // Total width of interface
	extern float nd_dx;

	extern float epsilon;
	extern float omega;
	extern float mob_phi;

	//Free energy of the form =  `sum_{a,b} (epsilon^2/2)(\nabla \phi_a)(\nabla \phi_b) + omega(\phi_a)(\phi_b)(1 - \alpha*c)
	if(!strcmp(init, "default"))
	{
		omega = 2.0*nd_se/(eta*nd_dx);
		epsilon = (4.0/PI)*sqrt(eta*nd_dx*nd_se);
		mob_phi = (PI*PI/16.0)*(nd_mobility/(eta*nd_dx));

		printf("Free energy parameters\n");
		printf("epsilon : %0.2e \t omega : %0.2e \t mob_phi : %0.2e\n",epsilon, omega, mob_phi);
	}
	else if(!strcmp(init, "equaldiff"))
	{
		printf("%s - Not implemented yet. Change the init parameter keyword for initialization\n", init);
		exit(2);
	}
	else if(!strcmp(init, "equalcomp"))
	{
		printf("%s -  Not implemented yet. Change the init parameter keyword for initialization\n", init);
		exit(2);
	}
	
}

void Printsettings()
{
	extern float nd_mobility, nd_se, nd_dx, nd_totaltime, nd_dt, beta;
	extern float real_mobility, real_se, real_dx, real_totaltime, real_dt;
	extern int gridx, gridy, gridsize, eta, grainforce;
	extern char *bc;

	printf("\n");
	printf("Settings for simulation:\n");
	printf("Parameter ---- Dimensional ---- Non Dimensional Values\n");
	printf("Mobility ---- %0.2e (m^4/Js)---- %0.2e\n", real_mobility, nd_mobility);
	printf("Surface Energy ---- %0.2e (J/m^2)---- %0.2e\n", real_se, nd_se);
	printf("dx ---- %0.2e (m)---- %0.2e\n", real_dx, nd_dx);
	printf("dt ---- %0.2e (s)---- %0.2e\n", real_dt, nd_dt);
	printf("Total Time ---- %0.2e (s)---- %0.2e\n", real_totaltime, nd_totaltime);
	printf("Diffusivity ---- %0.2e (m^2/s)---- %0.2e\n", real_diffusivity, nd_diffusivity);
	printf("Total Domain length ---- %0.2e (m)---- %0.2e\n", gridx*real_dx, gridx*nd_dx);

	printf("Other Settings\n");
	printf("Rectangular Grid : %d X %d points\n", gridx, gridy);
	printf("Boundary Conditions : %s\n", bc);
	return;
}


void Allocate()
{
	
	/*Extern Variable declarations*/
	extern node *grid;

	extern int gridsize;
	grid = node_alloc(gridsize);
	return;
}

void Setup(char *geom)
{
	/*
	Setup the geometry for 1D grain growth with two phase field parameters
	*/
	extern node *grid;
	extern int gridsize;
	extern int gridx, gridy;
	extern float comp;

	int index;
	if(!strcmp(geom, "linear"))
	{
		for(int i=0;i<gridx;i++)
		{
			for(int j = 0; j < gridy; j++)
			{
				index = j*gridx + i;
				grid[index].activegrain = (int *)malloc(1 * sizeof(int));
				grid[index].phi = (float *)malloc(1 * sizeof(float));
				grid[index].phase = (char *)malloc(1 * sizeof(char));
				grid[index].nactive = 1;
				if(i<gridx/2)
				{		
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = comp;
				}
				else					
				{
					grid[index].activegrain[0] = 2;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = comp;
				}
			}	
		}
	}
	else if(!strcmp(geom, "circle"))
	{
		printf("Setting the geometry to circle with radius 0.3*size of grid points. Change the radius from the initialize.c, if required\n");
		float radius = 0.3*gridx;
		if(gridy < 2*radius) 
		{
			printf("Grid Size not enough for the required geometry\n");
			exit(1);
		}
		for(int i=0;i<gridx;i++)
		{
			for(int j = 0; j < gridy; j++)
			{
				index = j*gridx + i;
				grid[index].activegrain = (int *)malloc(1 * sizeof(int));
				grid[index].phi = (float *)malloc(1 * sizeof(float));
				grid[index].phase = (char *)malloc(1 * sizeof(char));
				grid[index].nactive = 1;
				if ((i - gridx/2)*(i - gridx/2) + (j - gridy/2)*(j - gridy/2) < radius*radius) 
				{		
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = comp;
				}
				else					
				{
					grid[index].activegrain[0] = 2;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = comp;
				}
			}	
		}	
	}
	else 
	{
		printf("Geometry parameter not supported\n");
		exit(2);
	}
	
	return;
}



/*
Test Functions 
*/
void Check_initialstate()
{
	extern int gridsize;
	extern node *grid;

	int i;
	for(i=0;i<gridsize;i++)
	{
		printf("%d ", *grid[i].activegrain);
	}
	return;
}

void randomTest()
{
	node *grid1;
	extern int gridsize;

	grid1 = node_alloc(gridsize);
	printf("SUCCESS in Allocation\n");

	for(int i = 0; i< gridsize; i++)
	{
		node_dealloc(grid1[i]);	
	}
	
	printf("Successfully deallocated\n");
}

void FreeMemory()
{
	extern char *bc;
	free(bc);

	extern node *grid;
	extern int gridsize;
	for(int i = 0 ; i < gridsize ; i++) node_dealloc(grid[i]);
	free(grid);

}