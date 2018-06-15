#include<stdio.h>
#include "Global.h"
#include "Functions.h"
#include "nodefunctions.h"

void main()
{
	/*
	This code deals with normal grain growth in 1D with an intention to replicate Cahn's model in steady state
	*/
	char init[] = "default";
	int flag; 

	flag = Initialize(init);
	/*
	This function handles the parameter setting, 
	memory allocation and setup of the domain
	*/ 
	if(flag != 1)
	{
		printf("NOT INITIALIZED PROPERLY, Please check");
	}

	int t = 0;
	int timesteps = 0;
	extern float total_time,dt;	

	timesteps = (total_time/dt);
	printf("Time = %f, Total Timesteps: %d\n", total_time, timesteps);
	printf("\n");

	while(t < timesteps)
	{
		printf("Present Time = %f\n", t*dt);
		Solver();
		if(t == timesteps - 1) Output(t);
		t += 1;
	}
}