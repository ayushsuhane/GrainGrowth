#include "Functions.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "Global.h"
#include "nodefunctions.h"
#include<math.h>

int Readcontrols(char *filename)
{
	extern char init[100], geometry[100], phisolver[100], csolver[100], initialfile[100];
	extern int noutput;
	extern int flag_readfromfile, flag_writetofile;
	extern double time_constant;
	flag_readfromfile = 0; // 0 : Not read 1: Reading
	flag_writetofile = 1;
	FILE *fp;
	fp = fopen(filename, "r");
	if(fp == NULL) 
	{
		printf("File to read control parameters not  found");
		return 0;
	}
	printf("%s : Reading Control Parameters\n", filename);

	char tmpstr[20], tmpval[20];
	char tempbuffer[100];
	while(!feof(fp))
	{
		if(fgets(tempbuffer, 100, fp))
		{
			sscanf(tempbuffer, "%20s : %20[^;]", tmpstr, tmpval);
			tmpstr[strcspn(tmpstr, "\r\n")] = 0;
			tmpval[strcspn(tmpval, "\r\n")] = 0;

			if(!strcmp(tmpstr, "init")) strcpy(init,tmpval);
			if(!strcmp(tmpstr, "geometry")) strcpy(geometry,tmpval);
			if(!strcmp(tmpstr, "phisolver")) strcpy(phisolver,tmpval);
			if(!strcmp(tmpstr, "csolver")) strcpy(csolver,tmpval);
			/*To start from a restart file*/
			if(!strcmp(tmpstr, "readfromfile")) 
			{
				strcpy(initialfile,tmpval);
				flag_readfromfile = 1;
			}

			if(!strcmp(tmpstr, "flagreadfromfile")) flag_readfromfile = atoi(tmpval);
			if(!strcmp(tmpstr, "flagwritetofile")) flag_writetofile = atoi(tmpval);
			if(!strcmp(tmpstr, "noutput")) noutput = atoi(tmpval);
			if(!strcmp(tmpstr, "time_constant")) time_constant = atof(tmpval);
		}
	}
	printf("%d Visual Output will be generated\n",noutput);
	printf("Read Control Parameters\n");
	fclose(fp);
	return 1;
}

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
	Setup_single(geom);
	//Setup(geom);

	//randomTest();
	//Check_initialstate();
	return 1;

}



void Readinputs(char *init)
{
	/* 
	Extern Variable declarations
	*/
	extern double real_mobility;
	extern double real_se;
	extern double real_dx;
	extern double real_totaltime;
	extern double real_dt;
	extern double real_diffusivity;
	extern double real_beta; 	// Artificial force coefficient
	extern int beta_timestep;		//timestep after which beta gets activated
	extern int phi_timestep; 		//timestep to form the interfaces
	
	extern double alpha; 		// Parameter representing interaction between solute atoms and Grain Boundary
	extern double initcomp;			// Initial Composition
	extern double temperature;


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
	initcomp = 0.0; 
	temperature = 1273.0; //Initialize at 273 K
	real_dt = 0.01; //Just in case
	phi_timestep = 0;
	/******/
    extern char inputfilename[100];
	//char FILENAME[] = "../input/parameters.txt";
	FILE *fp;
	fp = fopen(inputfilename, "r");
	
	//Error check
	if(fp == NULL) 
	{
		printf("File to read parameters not  found");
		return;
	}
	printf("%s : Reading Parameters\n", inputfilename);

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
			if(!strcmp(tmpstr, "Composition")) initcomp = atof(tmpval);
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
	extern double temperature;
	double gas_constant = 8.314; // J/K
	double molar_volume = 7.1e-6; // 1/m^3
	/*Non dimensionalization Constants*/
	extern double c_length, c_time, c_diffusivity, c_energy, c_mobility;
	/********************************/
	extern double real_dx, real_diffusivity; 
	extern int eta;
	extern int gridx;
	//c_diffusivity = 1e-14; // m^2/s
	c_diffusivity = real_diffusivity;
	c_length = real_dx;
	//c_length = 1e-10;     // m
	c_energy = (gas_constant*temperature/molar_volume) ; // J/m^3
	//c_se = c_energy*c_length;
	c_time = (c_length*c_length)/c_diffusivity;
	c_mobility = 1/(c_time*c_energy);      //  m^3/(Js)
	
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
	//printf("Surface Energy Scale Constant : %0.2e\n", c_se);
	printf("Mobility Phi scaling constant : %0.2e\n", c_mobility);
	printf("Diffusion scaling constant : %0.2e\n", c_diffusivity);

	extern double nd_mobility, real_mobility;
	//extern double nd_se, real_se;
	extern double nd_dx, real_dx;
	extern double nd_dt, real_dt;
	extern double nd_totaltime, real_totaltime;
	extern double nd_diffusivity, real_diffusivity;
	extern double nd_beta, real_beta;

	extern double time_constant;

	printf("Time constant = %e\n", time_constant);
	//nd_mobility =  real_mobility/c_mobility;
	//nd_se = real_se/c_se;
	nd_dx = real_dx/c_length;
	//nd_dt = real_dt/c_time;
	nd_totaltime = real_totaltime/c_time;
	nd_diffusivity = real_diffusivity/c_diffusivity;
	nd_dt = nd_dx*nd_dx/(time_constant*nd_diffusivity);
	nd_beta = real_beta;

	return;
}

void Parametrize(char *init)
{
	// Define the parameters like epsilon and omega, mobility_phi here, which 
	// can be used in solver functions
	//extern double nd_se;
	//extern double nd_mobility;
	extern int eta; // half width of interface
	extern double real_dx;
	extern double real_se;
	extern double c_energy, c_length, c_mobility;
	//extern double nd_dx;

	extern double epsilon;
	extern double omega;
	extern double mob_phi;

	//Free energy of the form =  `sum_{a,b} (epsilon^2/2)(\nabla \phi_a)(\nabla \phi_b) + omega(\phi_a)(\phi_b)(1 - \alpha*c)
	if(!strcmp(init, "default"))
	{
		omega = 2.0*real_se/(eta*real_dx);
		epsilon = (4.0/PI)*sqrt(eta*real_dx*real_se);
		mob_phi = (PI*PI/16.0)*(real_mobility/(eta*real_dx));

		printf("Free energy parameters\n");
		printf("Real - epsilon : %0.2e \t omega : %0.2e \t mob_phi : %0.2e\n",epsilon, omega, mob_phi);
		omega = omega/c_energy;
		epsilon = epsilon/(c_length*sqrt(c_energy));
		mob_phi = mob_phi/c_mobility;
		printf("Scaled: -epsilon : %0.2e \t omega : %0.2e \t mob_phi : %0.2e\n",epsilon, omega, mob_phi);
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
	extern double nd_mobility, nd_se, nd_dx, nd_totaltime, nd_dt, beta;
	extern double real_mobility, real_se, real_dx, real_totaltime, real_dt;
	extern int gridx, gridy, gridsize, eta, grainforce;
	extern char *bc;

	printf("\n");
	printf("Settings for simulation:\n");
	printf("Parameter ---- Dimensional ---- Non Dimensional Values\n");
	//printf("Mobility ---- %0.2e (m^4/Js)---- %0.2e\n", real_mobility, nd_mobility);
//	printf("Surface Energy ---- %0.2e (J/m^2)---- %0.2e\n", real_se, nd_se);
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
	extern double initcomp;
	extern int flag_readfromfile;
	extern char initialfile[100];

	printf("File to readfrom : %s\n", initialfile);
	int index;
	if(!strcmp(geom, "linear") && !flag_readfromfile)
	{
		for(int i=0;i<gridx;i++)
		{
			for(int j = 0; j < gridy; j++)
			{
				index = j*gridx + i;
				//grid[index].activegrain = (int *)malloc(1 * sizeof(int));
				//grid[index].phi = (double *)malloc(1 * sizeof(double));
				//grid[index].phase = (char *)malloc(1 * sizeof(char ));
				grid[index].nactive = 1;
				if(i<gridx/2)
				{		
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
				else					
				{
					grid[index].activegrain[0] = 2;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
			}	
		}
	}
	else if(!strcmp(geom, "circle") && !flag_readfromfile)
	{
		printf("Setting the geometry to circle with radius 0.3*size of grid points. Change the radius from the initialize.c, if required\n");
		double radius = 0.3*gridx;
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
				//grid[index].activegrain = (int *)malloc(1 * sizeof(int));
				//grid[index].phi = (double *)malloc(1 * sizeof(double));
				//grid[index].phase = (char *)malloc(1 * sizeof(char));
				grid[index].nactive = 1;
				if ((i - gridx/2)*(i - gridx/2) + (j - gridy/2)*(j - gridy/2) < radius*radius) 
				{		
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
				else					
				{
					grid[index].activegrain[0] = 2;
					grid[index].phi[0] = 1;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
			}	
		}	
	}
	else if(flag_readfromfile)
	{
		/*Read file here*/
		FILE *initsetup = fopen(initialfile, "r");
		if (initsetup == NULL)
		{
			printf("FILE NOT FOUND. Resetting back to running from the scratch");
			flag_readfromfile = 0;
			Setup(geom);
			return;
		}
		extern int gridsize, gridx, gridy;
		for (int i=0; i<gridx ; i++)
		{
			for (int j=0; j<gridy; j++)
			{
				index = j*gridx + i;
				fscanf(initsetup, "%d\n", &index);
				fscanf(initsetup, "%d\n", &grid[index].nactive);
				/**Allocate memory for filling data**/
				//grid[index].activegrain = (int *)malloc(grid[index].nactive * sizeof(int));
				//grid[index].phi = (double *)malloc(grid[index].nactive * sizeof(double));
				//grid[index].phase = (char *)malloc(grid[index].nactive * sizeof(char));
				/*fill data*/
				for(int k=0; k< grid[index].nactive; k++)
					fscanf(initsetup, "%d ", &grid[index].activegrain[k]);
				fscanf(initsetup, "\n");
				for(int k=0; k< grid[index].nactive; k++)
					fscanf(initsetup, "%le ", &grid[index].phi[k]);
				fscanf(initsetup, "\n");
				for(int k=0; k< grid[index].nactive; k++)
					fscanf(initsetup, "%c ", &grid[index].phase[k]);
				fscanf(initsetup, "\n");
				fscanf(initsetup, "%le\n", &grid[index].comp);
				fscanf(initsetup, "\n");
			}
			fscanf(initsetup,"\n");
		}
		fclose(initsetup);
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
		printf("%d %lf\n", *grid[i].activegrain, *grid[i].phi);
	}
	return;
}

void randomTest()
{
	node *grid1;
	extern int gridsize;

	grid1 = node_alloc(gridsize);
	printf("SUCCESS in Allocation\n");
	/*
	for(int i = 0; i< gridsize; i++)
	{
		node_dealloc(grid1[i]);	
	}
	*/
	free(grid1);
	
	printf("Successfully deallocated\n");
}

void FreeMemory()
{
	extern char *bc;
	free(bc);

	extern node *grid;
	extern int gridsize;
	//for(int i = 0 ; i < gridsize ; i++) node_dealloc(grid[i]);
	free(grid);

}



/********************************************************/
/*****Single Setup *************************************/
/**********Only for one phi*****************************/
void Setup_single(char *geom)
{
	/*
	Setup the geometry for 1D grain growth with two phase field parameters
	*/
	extern node *grid;
	extern int gridsize;
	extern int gridx, gridy;
	extern double initcomp;
	extern int flag_readfromfile;
	extern char initialfile[100];

	printf("File to readfrom : %s\n", initialfile);
	int index;
	if(flag_readfromfile == 0) printf("Not reading any previous file");
	if(!strcmp(geom, "linear") && !flag_readfromfile)
	{
		for(int i=0;i<gridx;i++)
		{
			for(int j = 0; j < gridy; j++)
			{
				index = j*gridx + i;
				//grid[index].activegrain = (int *)malloc(1 * sizeof(int));
				//grid[index].phi = (double *)malloc(1 * sizeof(double));
				//grid[index].phase = (char *)malloc(1 * sizeof(char ));
				grid[index].nactive = 1;
				if(i<gridx/2)
				{		
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 1.0;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
				else					
				{
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 0.0;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
			}	
		}
	}
	else if(!strcmp(geom, "circle") && !flag_readfromfile)
	{
		printf("Setting the geometry to circle with radius 0.3*size of grid points. Change the radius from the initialize.c, if required\n");
		double radius = 0.3*gridx;
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
				//grid[index].activegrain = (int *)malloc(1 * sizeof(int));
				//grid[index].phi = (double *)malloc(1 * sizeof(double));
				//grid[index].phase = (char *)malloc(1 * sizeof(char));
				grid[index].nactive = 1;
				if ((i - gridx/2)*(i - gridx/2) + (j - gridy/2)*(j - gridy/2) < radius*radius) 
				{		
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 1.0;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
				else					
				{
					grid[index].activegrain[0] = 1;
					grid[index].phi[0] = 0.0;
					grid[index].phase[0] = 'f';
					grid[index].comp = initcomp;
				}
			}	
		}	
	}
	else if(flag_readfromfile)
	{
		/*Read file here*/
		FILE *initsetup = fopen(initialfile, "r");
		if (initsetup == NULL)
		{
			printf("FILE NOT FOUND. Resetting back to running from the scratch");
			flag_readfromfile = 0;
			Setup(geom);
			return;
		}
		extern int gridsize, gridx, gridy;
		for (int i=0; i<gridx ; i++)
		{
			for (int j=0; j<gridy; j++)
			{
				index = j*gridx + i;
				fscanf(initsetup, "%d\n", &index);
				fscanf(initsetup, "%d\n", &grid[index].nactive);
				/**Allocate memory for filling data**/
				//grid[index].activegrain = (int *)malloc(grid[index].nactive * sizeof(int));
				//grid[index].phi = (double *)malloc(grid[index].nactive * sizeof(double));
				//grid[index].phase = (char *)malloc(grid[index].nactive * sizeof(char));
				/*fill data*/
				for(int k=0; k< grid[index].nactive; k++)
					fscanf(initsetup, "%d ", &grid[index].activegrain[k]);
				fscanf(initsetup, "\n");
				for(int k=0; k< grid[index].nactive; k++)
					fscanf(initsetup, "%le ", &grid[index].phi[k]);
				fscanf(initsetup, "\n");
				for(int k=0; k< grid[index].nactive; k++)
					fscanf(initsetup, "%c ", &grid[index].phase[k]);
				fscanf(initsetup, "\n");
				fscanf(initsetup, "%le\n", &grid[index].comp);
				fscanf(initsetup, "\n");
			}
			fscanf(initsetup,"\n");
		}
		fclose(initsetup);
	}
	else 
	{
		printf("Geometry parameter not supported\n");
		exit(2);
	}
	
	return;
}