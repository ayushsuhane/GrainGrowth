#include<stdio.h>
#include<stdlib.h>
#include"Global.h"
#include"nodefunctions.h"
#include"Functions.h"



void Solver()
{
	/*External Variables*/
	extern int gridsize;
	extern node *grid;
	extern float dt;
	extern float beta;

	int i, active_grain;
	int intf_flag = 1;

	float lambda = 0;
	float *derivative;
	float laplace = 0;
	float sumterms=0;
	float *term ;
	node list;

	node *new_grid;
	new_grid = node_alloc(gridsize);
	/*
	Drichlet Boundary conditions
	End points fixed
	*/
	for(i=1;i<gridsize-1;i++)
	{
		intf_flag = node_checkgradient(grid[i-1], grid[i], grid[i+1]); //Checks whether the gradient in terms of phi is present

		//if(intf_flag) printf("Grid ID : %d, left grain : %d, right grain : %d, mid grain : %d\n", i, grid[i-1].activegrain[0], grid[i+1].activegrain[0], grid[i].activegrain[0]);

		if(!intf_flag)
		{
			new_grid[i] = node_update(grid[i], i);
			continue; // skips the point if this point is not a apart of the interface
		}

		sumterms = 0;
		lambda = 0; // Langrange Parameter

		/* If at the interface, build a list to calculate laplacian*/
		

		list = node_buildlist(grid[i-1], grid[i], grid[i+1]);
		//node_print(list, i);

		//for(int j = 0; j< list.nactive; j++) printf("Phi : %f, GrainID: %d, gridid : %d \n", list.phi[j], list.activegrain[j], i);
		/*Allocate an array for storing the laplacian of all the active phase fields*/
		term = (float *)malloc(list.nactive*sizeof(float));
		for(int temp = 0; temp < list.nactive; temp++) term[temp] = 0.0;
		/*
		Loop over all the active phase fields
		*/
		for(int j = 0; j < list.nactive; j++)
		{
			laplace = node_laplacian1D(list.activegrain[j], grid[i-1], grid[i], grid[i+1], dx);
			//printf("Laplacian : %f, GrainNum : %d, nodenum: %d\n", laplace, list.activegrain[j], i);
			term[j] = node_interface(laplace, list.phi[j]) + node_addforce(list.activegrain[j], list.phi[j], beta);
			//printf("InterfaceTerm: %f, ExtraForce : %f, GrainID: %d\n", node_interface(laplace, list.phi[j]), node_addforce(list.activegrain[j], list.phi[j]), list.activegrain[j] );
			sumterms += term[j]; 
		}
		//substract the lambda
		lambda = sumterms/list.nactive;
		for(int j = 0; j < list.nactive; j++) term[j] -= lambda;
		for(int j = 0; j< list.nactive; j++) 
		{
			/*term represents phi in PF eqution here*/
			term[j] = dt*term[j] + list.phi[j];
			//printf("Phi : %f, GrainID: %d, gridid : %d \n", list.phi[j], list.activegrain[j], i);
			/*change the value of list as old list.phi[j] is not needed anymore*/
			list.phi[j] = term[j];
		}

		list = node_simplex(list);
		/*An updated list for the node*/
		new_grid[i] = node_update(list, i); //Updates the list into new grid
		node_print(new_grid[i], i);

		free(term);
		node_dealloc(list);
	}
	//node_print(grid[497], 497);
	for(int i = 1; i< gridsize-1; i++) 
	{
		grid[i] = node_transfer(new_grid[i], grid[i]);
		node_dealloc(new_grid[i]);
	}
	return;
}

