#include<stdio.h>
#include "Global.h"
#include "Functions.h"
#include "nodefunctions.h"
#include "analysisfunctions.h"
#include<stdlib.h>
#include<string.h>

int main()
{
	/*
	This code deals with normal grain growth in 1D with an intention to replicate Cahn's model in steady state
	*/
	// Initialization string
	/*    'default', 'equalcomp', 'equaldiff'    */
	extern char inputfilename[100], controlinputfilename[100];
	strcpy(inputfilename,"../input/parameters.txt");
	strcpy(controlinputfilename,"../input/control.txt");

	//char init[] = "default";
	
	//Geometry intitialization flag 
	/*   'linear', 'circle'   */
	//char geometry[] = "linear"; 
	
	// Analysis string to track the point with phi = 0.5 for grain 
	// interface between grain 1 and 2 
	/*    'NULL', 'track 1 2 0.5'    */
	char analysisinit[] = "track : 0.5"; 
	
	// MPF for multiphase field formulation
	/*    'MPF'    */
	//char phisolver[] = "MPF"; 
	
	// Solver string for composition solver
	// 'equalcomp' for equal composition
	// 'equaldiff' for equal diffusion potential 
	// 'None' for not solving composition
	/*    'None', 'equalcomp', 'equaldiff'    */
	//char csolver[] = "equalcomp"; 
	
	// To turn on the keywords associalted with Initialize
	int flag;
	int flagcontrol = Readcontrols(controlinputfilename);
	if(flagcontrol != 1) 
		{
			printf("Couldnot read the controls file");
			exit(2);
		}
	extern int noutput;
	extern char init[100], geometry[100], phisolver[100], csolver[100]; 
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
	analysisflag = Analysisinit_single(analysisinit);
	

	if(flag != 1 && analysisflag != 1)
	{
		printf("NOT INITIALIZED PROPERLY, Please check");
		exit(0);
	}

	int t = 0;
	long unsigned int timesteps = 0;
	double time = 0.0;
	extern double nd_totaltime,nd_dt;	
	extern double real_totaltime;
	extern double nd_dt;
	//extern double dt, nd_dt, nd_dx, nd_diffusivity;
	//extern int phi_timestep, beta_timestep;

	//int count = 0;
	timesteps = (int)(nd_totaltime/nd_dt);
	printf("Real Values : Total Time = %0.2e, Timestep : %0.2e, Number of Timesteps: %lu\n", real_totaltime, nd_dt*real_totaltime/nd_totaltime,  timesteps);
	printf("\n");
	printf("Non Dimensional Values : Time = %e, Timestep : %0.2e, Number of Timesteps : %lu\n", nd_totaltime, nd_dt, timesteps);


	//while(time < nd_totaltime)
	while(t < timesteps + 1)
	{
		/**For numerical reasons**/
		/*
		if(count<(phi_timestep+beta_timestep)) 
		{
			dt = nd_dx*nd_dx/(100.0*nd_diffusivity);
			count += 1;
		}
		else dt = nd_dt;

		if(count == phi_timestep+beta_timestep) 
		{
			t = 0;
			count += 1;
		}
		*/
		/*************************/
		if (t%10000 == 0) printf("Present Time = %0.2e, timestep = %d\n", t*nd_dt*real_totaltime/nd_totaltime, t);
		
		// Main solver function //
		/*    
		phisolver: string with values to solve phase field parameter
		csolver  : string with values to solve the composition field
		*/
		Solver_single(phisolver, csolver, t);


		// All the analysis included in the Analysic function are
		// evaluated. Currently include only tracking of the grain boundary
		// at a certain value of phi
		if(t%10000 == 0) 
		{
			Analysis_single(t*nd_dt*real_totaltime/nd_totaltime);
			printf("Analysis completed at timestep %d\n", t);
		}
		// Returns desired outputs in files
		if(t == timesteps) 
		{
			Output(t);
			//printf("Outside Output");
		}

		// Different function for output
		if(t % ((int)((timesteps+1)/noutput)) == 0) 
		{
			output_visual(t, 1);
			//printf("Outside Visual Output\n");
		}
		/*
		if(t > 4000000)
		{
			t += 5;
			nd_dt *= 5;
		}
		*/
		//else
		//{
			t += 1;
		//}
		
	}

	// End the analysis function
	Analysisend();
	// Frees all the allocated memory
	FreeMemory();
	return 0;
}