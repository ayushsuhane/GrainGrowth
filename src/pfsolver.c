#include<stdio.h>	
#include<stdlib.h>
#include"Global.h"
#include"nodefunctions.h"
#include"Functions.h"
#include<math.h>


void Solver(char *phisolver, char *csolver, int t)
{
	/*External Variables*/
	extern int gridx, gridy;
	extern int gridsize;
	extern node *grid;
	extern double nd_dx;
	extern char *bc;
	extern int phi_timestep;
	extern int beta_timestep;
	extern double nd_dt;

	node list;

	int up, left, centre, right, down, index;
	node *new_grid;
	new_grid = node_alloc(gridsize);

	/*Just a way to handle single nodes in y direction for 1D study*/
	int county_start, county_end;
	if(gridy == 1)
	{
		county_start = 0;
		county_end = 1;
	}
	else
	{
		county_start = 1;
		county_end = gridy - 1;
	}
	//printf("Gridy points start: %d stop :%d\n",county_start, county_end);
	//For a smaller grid, boundary conditions are dealt later
	for(int gidx = 1;gidx < gridx-1;gidx++)
	{

		for(int gidy = county_start;gidy < county_end;gidy++)
		{
			//Assignment of stencil indices
			index = gidy*gridx + gidx;
			//printf("Index : %d\n", index);
			if(gridy < 2 )
			{
				//Assign centre, up, down, left, right index
				centre = index;
				up = index;
				down = index;
				left = gidy*gridx + gidx - 1;
				right = gidy*gridx + gidx + 1;
			}
			else
			{
				centre = index;
				up = (gidy+1)*gridx + gidx;
				down = (gidy-1)*gridx + gidx;
				left = gidy*gridx + gidx - 1;
				right = gidy*gridx + gidx + 1;
			}
			//Check the gradient. If no gradient exist, don't solve
			new_grid[index] = node_phisolver(phisolver, index, grid[up], grid[left], grid[centre], grid[right], grid[down], new_grid[index], nd_dx, t);
			if(t > phi_timestep) new_grid[index].comp = node_csolver(csolver, index, grid[up], grid[left], grid[centre], grid[right], grid[down], new_grid[index], nd_dx, t);
			new_grid[index] = node_update(new_grid[index], index);
			/**Check for potential*/
			/*********************/
			if((index == 101 || index == 100|| index == 109 || index == 108 || index == 107 || index == 106 || index == 104 || index == 105 || index == 103) && (t%10000 == 0))
			{
				double sumphi=0, sumgrad=0;
				extern double nd_dx;
				for(int i=0; i<grid[index].nactive; i++)
				{
					sumgrad += (1-grid[index].phi[i])*(node_gradient2D(grid[index].activegrain[i], grid[right], grid[left], nd_dx));
					for(int j=0; j< grid[index].nactive;j++)
					{
						if(i!=j) sumphi += grid[index].phi[i]*grid[index].phi[j];
					}
				}
				
				double gradcx = (grid[right].comp - grid[left].comp)/(2.0*nd_dx);

				extern double nd_diffusivity;  
				extern double omega, alpha;
				double potential = log(grid[index].comp/(1-grid[index].comp)) - omega*alpha*sumphi;
				double H = gradcx - grid[index].comp*(1-grid[index].comp)*omega*alpha*sumgrad;
				printf("Index : %d Potential: %e H:%e\n", index, potential, H);		
			}
			
			/*********************/

		}
	}
	//For corners, copy any of new_grid node, for rest copy the consecutive new_grid value
	if(gridy > 2) node_handleBC(gridx, gridy, new_grid, bc); //drichlet copies the boundary element with its consecutive entry
	else
	{
		//just handle x axis`
		new_grid[0] = node_update(new_grid[gridx - 2], 1);
		new_grid[gridx - 1] = node_update(new_grid[1], gridx - 2);
	}
	//Transfer the new_grid to grid : node by node
	for(int gidx = 0;gidx < gridx ;gidx++)
	{
		for(int gidy = 0;gidy < gridy ;gidy++)
		{
			index = gidy*gridx + gidx;
			grid[index] = node_transfer(new_grid[index], grid[index]);
			//node_dealloc(new_grid[index]);
		}
	}
	free(new_grid);
	return;
}	
		
		




/******************************************************************************/
/*************Single***********************************************************/
void Solver_single(char *phisolver, char *csolver, int t)
{
	/*External Variables*/
	extern int gridx, gridy;
	extern int gridsize;
	extern node *grid;
	extern double nd_dx;
	extern char *bc;
	extern int phi_timestep;
	extern int beta_timestep;
	extern double nd_dt;

	node list;

	int up, left, centre, right, down, index;
	int nextleft, nextright;
	node *new_grid;
	new_grid = node_alloc(gridsize);

	/*Just a way to handle single nodes in y direction for 1D study*/
	int county_start, county_end;
	if(gridy == 1)
	{
		county_start = 0;
		county_end = 1;
	}
	else
	{
		county_start = 1;
		county_end = gridy - 1;
	}
	//printf("Gridy points start: %d stop :%d\n",county_start, county_end);
	//For a smaller grid, boundary conditions are dealt later
	for(int gidx = 2;gidx < gridx-2;gidx++)
	{

		for(int gidy = county_start;gidy < county_end;gidy++)
		{
			//Assignment of stencil indices
			index = gidy*gridx + gidx;
			//printf("Index : %d\n", index);
			if(gridy < 2 )
			{
				//Assign centre, up, down, left, right index
				centre = index;
				up = index;
				down = index;
				left = gidy*gridx + gidx - 1;
				right = gidy*gridx + gidx + 1;
				nextleft = gidy*gridx + gidx - 2;
				nextright = gidy*gridx + gidx + 2;

			}
			else
			{
				centre = index;
				up = (gidy+1)*gridx + gidx;
				down = (gidy-1)*gridx + gidx;
				left = gidy*gridx + gidx - 1;
				right = gidy*gridx + gidx + 1;
				nextleft = gidy*gridx + gidx - 2;
				nextright = gidy*gridx + gidx + 2;

			}
			//Check the gradient. If no gradient exist, don't solve
			new_grid[index] = node_phisolver_single(phisolver, index, grid[up], grid[nextleft], grid[left], grid[centre], grid[right], grid[nextright], grid[down], new_grid[index], nd_dx, t);
			if(t > phi_timestep) new_grid[index].comp = node_csolver_singlecons(csolver, index, grid[up], grid[nextleft], grid[left], grid[centre], grid[right], grid[nextright], grid[down], nd_dx, t);
			else new_grid[index].comp = grid[centre].comp;
			/**Check for potential*/
			/*********************/
			/*********************/

		}
	}
	//For corners, copy any of new_grid node, for rest copy the consecutive new_grid value
	if(gridy > 2) node_handleBC(gridx, gridy, new_grid, bc); //drichlet copies the boundary element with its consecutive entry
	else
	{
		//just handle x axis
		new_grid[0] = node_update_single(new_grid[4]);
		new_grid[1] = node_update_single(new_grid[4]);
		new_grid[2] = node_update_single(new_grid[4]);
		new_grid[3] = node_update_single(new_grid[4]);
		
		//new_grid[2] = node_update_single(new_grid[3]);
		new_grid[gridx - 1] = node_update_single(new_grid[gridx - 5]);
		new_grid[gridx - 2] = node_update_single(new_grid[gridx - 5]);
		new_grid[gridx - 3] = node_update_single(new_grid[gridx - 5]);
		new_grid[gridx - 4] = node_update_single(new_grid[gridx - 5]);
		
		//new_grid[gridx - 3] = node_update_single(new_grid[gridx - 4]);
		//new_grid[0].comp = initcomp;
		//new_grid[1].comp = initcomp;
		//new_grid[gridx - 1].comp = initcomp;
		//new_grid[gridx - 2].comp = initcomp;
		//new_grid[0].comp = new_grid[gridx - 2].comp;
		//new_grid[gridx - 1].comp = new_grid[1].comp;
	}
	
	/**Impose Mass Balance**/
	/*
	double sumcomp = 0.0;
	extern double initcomp;
	for (int i=0; i<gridx ; i++)
		for(int j=0; j<gridy; j++)
			sumcomp += new_grid[j*gridx + i].comp;

	for (int i=0; i<gridx ; i++)
		for(int j=0; j<gridy; j++)
			new_grid[j*gridx + i].comp += (initcomp - sumcomp/(gridx*gridy));
	*/
	/***********************/
	//Transfer the new_grid to grid : node by node
	for(int gidx = 0;gidx < gridx ;gidx++)
	{
		for(int gidy = 0;gidy < gridy ;gidy++)
		{
			index = gidy*gridx + gidx;
			grid[index] = node_update_single(new_grid[index]);
			//node_dealloc(new_grid[index]);
		}

	}
	free(new_grid);
	return;
}	