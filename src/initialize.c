#include "Functions.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "Global.h"
#include "nodefunctions.h"

int Initialize(char *init)
{
	printf("Initialization : %s\n", init);
	/*
	Default initialization for two grains in a 1D subsystem
	*/
	/*
	Read the values of mobility, surface energy, number of grid points, gridspacing
	*/
	Readinputs();
	/*
	Print Dimensionalization consants for scaling back to real values
	*/
	Dimensionalize();
	/*
	Allocate the relevant memory 
	*/
	Allocate();
	/*
	Setup the geometry
	*/
	Setup();


	//randomTest();
	//Check_initialstate();
	return 1;

}

void Readinputs()
{
	/* 
	Extern Variable declarations
	*/
	extern float mobility;
	extern float se;
	extern float dx;
	extern float total_time;
	extern float dt;
	extern float beta;

	extern int gridsize;
	extern int eta;
	extern int grainforce;
	grainforce = 0;
	beta = 0.0;

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
	printf("%s : Reading Parameters", FILENAME);

	char tmpstr[15], tmpval[15];
	char tempbuffer[100];
	while(!feof(fp))
	{
		if(fgets(tempbuffer, 100, fp))
		{
			sscanf(tempbuffer, "%15s : %15[^;]", tmpstr, tmpval);
			if(!strcmp(tmpstr, "Mobility")) mobility = atof(tmpval);
			if(!strcmp(tmpstr, "SurfaceEnergy")) se = atof(tmpval);
			if(!strcmp(tmpstr, "dx")) dx = atof(tmpval);
			if(!strcmp(tmpstr, "GridSize")) gridsize = atoi(tmpval);
			if(!strcmp(tmpstr, "Eta")) eta = atoi(tmpval);
			if(!strcmp(tmpstr, "TotalTime")) total_time = atof(tmpval);
			if(!strcmp(tmpstr, "dt")) dt = atof(tmpval);
			if(!strcmp(tmpstr, "GrainForce")) grainforce = atoi(tmpval);
			if(!strcmp(tmpstr, "Beta")) beta = atof(tmpval);
		}
	}
	fclose(fp);
	return;
}

void Dimensionalize()
{
	/*
	Contains the non-dimensional parameters. Check the input with these parameters first
	*/
	float length, energy, mob, time;
	length = 1e-9;
	mob = 1e-14;
	time = (length*length)/mob;
	energy = 1.0/length;

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
	printf("Mobility scaling constant : %e\n", mob);
	return;
}

void Allocate()
{
	
	/*Extern Variable declarations*/
	extern node *grid;

	extern int gridsize;


	/*************************/
	//grid = (node *)malloc(gridsize * sizeof(node)); //Only one dimensional assignment
	

	//Initialize
	int i = 0;
	// for(i=0;i<gridsize;i++)
	// {
		
	// 	This is an array of pointers each pointing towards a node type variable
	// 	can be accessed through 
	// 	grid[i].activegrain, grid[i].phase, grid[i].phi, grid[i].nactive
		
	// 	grid[i] = node_initialize();
	// 	new_grid[i] = node_initialize();
	// }
	grid = node_alloc(gridsize);
	

	return;
}

void Setup()
{
	/*
	Setup the geometry for 1D grain growth with two phase field parameters
	*/
	extern node *grid;
	extern int gridsize;

	int i;
	for(i=0;i<gridsize;i++)
	{
		grid[i].activegrain = (int *)malloc(1 * sizeof(int));
		grid[i].phi = (float *)malloc(1 * sizeof(float));
		grid[i].phase = (char *)malloc(1 * sizeof(char));
		grid[i].nactive = 1;
		if(i<gridsize/2)
		{
			grid[i].activegrain[0] = 1;
			grid[i].phi[0] = 1;
			grid[i].phase[0] = 'f';
		}
		else
		{
			grid[i].activegrain[0] = 2;
			grid[i].phi[0] = 1;
			grid[i].phase[0] = 'f';
		}
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