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
	Dimensionalize();

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
	extern float mobility;
	extern float se;
	extern float dx;
	extern float total_time;
	extern float dt;
	extern float beta;  		// Artificial force coefficient
	extern float alpha; 		// Parameter representing interaction between solute atoms and Grain Boundary
	extern float diffusivity;
	extern float comp;			// Initial Composition
	extern float temperature;


	extern int gridx, gridy;
	extern int gridsize;
	extern int eta;				// Number of grid points for the interface
	extern int grainforce;		// Grain in which additional artificial force needs to be activated

	extern char *bc;			// String for boundary condition

	bc = (char *)malloc(15*sizeof(char));
	
	grainforce = 0;
	beta = 0.0;
	alpha = 0.0; //If not present, then no segregation potential
	diffusivity = 0; //If not specified
	comp = 0.0; 
	temperature = 273.0; //Initialize at 273 K
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

			if(!strcmp(tmpstr, "Mobility")) mobility = atof(tmpval);
			if(!strcmp(tmpstr, "SurfaceEnergy")) se = atof(tmpval);
			if(!strcmp(tmpstr, "dx")) dx = atof(tmpval);
			if(!strcmp(tmpstr, "GridX")) gridx = atoi(tmpval);
			if(!strcmp(tmpstr, "GridY")) gridy = atoi(tmpval);
			if(!strcmp(tmpstr, "Eta")) eta = atoi(tmpval);
			if(!strcmp(tmpstr, "TotalTime")) total_time = atof(tmpval);
			if(!strcmp(tmpstr, "dt")) dt = atof(tmpval);
			if(!strcmp(tmpstr, "GrainForce")) grainforce = atoi(tmpval);
			if(!strcmp(tmpstr, "Beta")) beta = atof(tmpval);
			if(!strcmp(tmpstr, "Alpha")) alpha = atof(tmpval);
			if(!strcmp(tmpstr, "Diffusivity")) diffusivity = atof(tmpval);
			if(!strcmp(tmpstr, "Composition")) comp = atof(tmpval);
			if(!strcmp(tmpstr, "Temperature")) temperature = atof(tmpval);
			if(!strcmp(tmpstr, "BC")) strcpy(bc,tmpval);
		}
	}
	fclose(fp);
	gridsize = gridx*gridy;
	return;
}

void Dimensionalize()
{
	/*
	Contains the non-dimensional parameters. Check the input with these parameters first
	*/
	extern float temperature;
	float gas_constant = 8.314; // J/K
	float molar_volume = 1.06e-5; // 1/m^3
	float length, energy, mob, time;
	float diffusion;
	float se_scaling;
	diffusion = 1e-15; // m^2/s
	length = 1e-9;     // m
	energy = (gas_constant*temperature/molar_volume) ; // J/m^2
	se_scaling = energy*length;
	time = (length*length)/diffusion;
	mob = time*length*length/se_scaling;      //  m^4/(Js)
	
	printf("Non Dimensional Constants: \n");
	printf("For more information : check Dimensionalize function\n");
	/*
	These constants are defined to scale the quantities back to the real values 
	such that 
	A = a* X a
	where the values here are a*, the values in parameters.txt are a and actual values can be obtained 
	by multiplying these values depending on their respective scaling behaviour.
	*/
	printf("Length Scale : %e\n", length);
	printf("Time Scale : %e\n", time);
	printf("Energy Scale : %e\n", energy);
	printf("Surface Energy Scale Constant : %e\n", se_scaling);
	printf("Mobility scaling constant : %e\n", mob);
	printf("Diffusion scaling constant : %e\n", diffusion);

	return;
}

void Parametrize(char *init)
{
	// Define the parameters like epsilon and omega, mu here, which 
	// can be used in solver functions
	extern float se;
	extern float mobility;
	extern int eta; // Total width of interface
	extern float dx;

	extern float epsilon;
	extern float omega;
	extern float mob_phi;

	//Free energy of the form =  `sum_{a,b} (epsilon^2/2)(\nabla \phi_a)(\nabla \phi_b) + omega(\phi_a)(\phi_b)(1 - \alpha*c)
	if(!strcmp(init, "default"))
	{
		omega = 2.0*se/(eta*dx/2.0);
		epsilon = (4.0/PI)*sqrt(eta*dx*se/2.0);
		mob_phi = (PI*PI/16.0)*(mobility/(eta*dx/2.0));

		printf("Free energy parameters\n");
		printf("epsilon : %f \t omega : %f \t mob_phi : %f\n",epsilon, omega, mob_phi);
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
	extern float mobility, se, dx, total_time, dt, beta;
	extern int gridx, gridy, gridsize, eta, grainforce;
	extern char *bc;

	printf("\n");
	printf("Settings for simulation:\n");
	printf("Rectangular Grid : %d X %d points\n", gridx, gridy);
	printf("Grid Spacing : %f\n", dx);
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