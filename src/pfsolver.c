#include<stdio.h>	
#include<stdlib.h>
#include"Global.h"
#include"nodefunctions.h"
#include"Functions.h"



void Solver(char *phisolver, char *csolver)
{
	/*External Variables*/
	extern int gridx, gridy;
	extern int gridsize;
	extern node *grid;
	extern float dx;
	extern char *bc;

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
			new_grid[index] = node_phisolver(phisolver, index, grid[up], grid[left], grid[centre], grid[right], grid[down], new_grid[index], dx);
			new_grid[index] = node_csolver(csolver, index, grid[up], grid[left], grid[centre], grid[right], grid[down], new_grid[index], dx);
		}
	}
	//For corners, copy any of new_grid node, for rest copy the consecutive new_grid value
	if(gridy > 2) node_handleBC(gridx, gridy, new_grid, bc); //drichlet copies the boundary element with its consecutive entry
	else
	{
		//just handle x axis
		new_grid[0] = node_update(new_grid[1], 1);
		new_grid[gridx - 1] = node_update(new_grid[gridx - 2], gridx - 2);
	}
	//Transfer the new_grid to grid : node by node
	for(int gidx = 0;gidx < gridx ;gidx++)
	{
		for(int gidy = 0;gidy < gridy ;gidy++)
		{
			index = gidy*gridx + gidx;
			grid[index] = node_transfer(new_grid[index], grid[index]);
			node_dealloc(new_grid[index]);
		}
	}
	free(new_grid);
	return;
}	
		
		
