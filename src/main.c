#include<stdio.h>
#include "Global.h"
#include "Functions.h"
#include "nodefunctions.h"
#include "analysisfunctions.h"
#include <stdlib.h>

void main()
{
	/*
	This code deals with normal grain growth in 1D with an intention to replicate Cahn's model in steady state
	*/
	// Initialization string
	/*    'default', 'equalcomp', 'equaldiff'    */
	char init[] = "default";
	
	//Geometry intitialization flag 
	/*   'linear', 'circle'   */
	char geometry[] = "linear"; 
	
	// Analysis string to track the point with phi = 0.5 for grain 
	// interface between grain 1 and 2 
	/*    'NULL', 'track 1 2 0.5'    */
	char analysisinit[] = "track : 1 2 0.5"; 
	
	// MPF for multiphase field formulation
	/*    'MPF'    */
	char phisolver[] = "MPF"; 
	
	// Solver string for composition solver
	// 'equalcomp' for equal composition
	// 'equaldiff' for equal diffusion potential 
	// 'None' for not solving composition
	/*    'None', 'equalcomp', 'equaldiff'    */
	char csolver[] = "equalcomp"; 
	
	// To turn on the keywords associalted with Initialize
	int flag; 
	
	// To check if analysis keyword was activated
	int analysisflag; 
	
	// Initialize the simulation 
	/* 
	flag = 0 => Initialize was unsuccessful
		 = 1 => Successfully Initialized
	*/
	flag = Initialize(init, geometry);

	// Initialize the analysis 
	/* 
	flag = 0 => Initialize was unsuccessful Analysis
		 = 1 => Successfully Initialized Analysis
	*/
	analysisflag = Analysisinit(analysisinit);
	

	if(flag != 1)
	{
		printf("NOT INITIALIZED PROPERLY, Please check");
		exit(0);
	}

	int t = 0;
	int timesteps = 0;
	extern float total_time,dt;	

	timesteps = (total_time/dt);
	printf("Time = %f, Total Timesteps: %d\n", total_time, timesteps );
	printf("\n");

	while(t < timesteps + 1)
	{
		if (t%1000 == 0) printf("Present Time = %f\n", t*dt);
		
		// Main solver function //
		/*    
		phisolver: string with values to solve phase field parameter
		csolver  : string with values to solve the composition field
		*/
		Solver(phisolver, csolver);

		// All the analysis included in the Analysic function are
		// evaluated. Currently include only tracking of the grain boundary
		// at a certain value of phi
		if(t%1000 == 0) Analysis(t*dt);

		// Returns desired outputs in files
		if(t == timesteps - 1) Output(t);

		// Different function for output
		if(t % ((int)(timesteps/2)) == 0) output_visual(t, 1);
		t += 1;
	}

	// End the analysis function
	Analysisend();
	// Frees all the allocated memory
	FreeMemory();
}