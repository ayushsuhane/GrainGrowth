#include<stdio.h>
#include"nodefunctions.h"
#include<stdlib.h>

node node_initialize()
{
	// Creating a node pointer and assining values
	// returning the dereferenced value to assing the initialized pointer 
	// address. The return type assigns the dereferenced pointer such that the
	// handling variable points to the node initialized.
	node *a;
	if((a = malloc(sizeof *a)) != NULL)
	{
		a->nactive = 0;
		a->activegrain = NULL;
		a->phi = NULL;
		a->phase = NULL;
	}
	return *a;

}

node * node_alloc(int size)
{	
	//Does a same job to node_initialize except for a given size 
	// and returns the array of pointers

	node *a;
	/*
	Memory Handler
	*/
	if((a = (node *)malloc(size * sizeof(node))) != NULL)
	{
		for(int i = 0; i < size; i++)
		{	
			a[i].nactive = 0;
			a[i].activegrain = NULL;
			a[i].phi = NULL;
			a[i].phase = NULL;
		}
	}
	else printf("Different Types");
	
	return a;
}

void node_dealloc(node  free_node)
{
	/*
	Takes in a pointer to node array 
	to deallocate the memory
	*/
	int i = 0;
	/* First free the memory allocated by pointers of node -  active_grain, phi, phase.
	Then clear the memory occupied by node. Generally preferred to go in a reverse order, otherwise 
	causes heap memory leak
	*/
		free(free_node.phi);
		free(free_node.phase);
		free(free_node.activegrain);
		//Check Whether to free the free_node as well??????
	//free(free_node);

	return;
}

int node_checkgradient(node left, node mid, node right)
{
	/*Check the gradient by checking the number of active grains*/
	/*
	returns 1 : if active grain of any of the node is more than 2
			  : if active grain is 1 and activegrain[0] is different
	*/

	if(left.nactive > 1 || mid.nactive > 1 || right.nactive > 1 ) return 1;
	else if(left.activegrain[0] != mid.activegrain[0] || left.activegrain[0] != right.activegrain[0]) return 1;
	else return 0; 
}


node  node_buildlist(node left, node mid, node right)
{
	/*
	Takes in neighbouring nodes, with sorted lists, and append the list to be checked in the result list
	*/
	// Allocate memory for result
	node *result;

	float SIZE_f = sizeof(float);
	int SIZE_i = sizeof(int);
	char SIZE_c = sizeof(char);

	int count = 0, change_flag = 0; //length of list
	result = node_alloc(1);

	//Start by copying the entries in mid node
	// Allocate size to copy the list
	result[0].activegrain = (int *)malloc(5 * SIZE_i); //5 is the maximum limit of neihbourhood grains
	result[0].phi = (float *)malloc(5 * SIZE_f);
	result[0].phase = (char *)malloc(5 * SIZE_c);
	//Assignment
	for(int i=0;i < mid.nactive; i++)
	{
		result[0].activegrain[i] = mid.activegrain[i];
		result[0].phi[i] = mid.phi[i];
		result[0].phase[i] = mid.phase[i];
		count += 1;
	}
	result[0].nactive = mid.nactive;
	//Now move to the points to check any unknown active grain
	/*Go to left*/
	/****************************Repeatative Code : can be built into a different function*************************/
	for (int i=0; i < left.nactive; i++)
	{	
		change_flag = 1;
		for(int j = 0; j< result[0].nactive; j++) 
		{
			if(result[0].activegrain[j] == left.activegrain[i])
			{
				change_flag = 0;
				break;
			} 	
		}		 //Break the current loop and go to next indices in left
		if(change_flag == 1)
		{
			result[0].activegrain[count] = left.activegrain[i];
			result[0].phi[count] = 0.0;
			result[0].phase[count] = left.phase[i];
			result[0].nactive += 1;
			count += 1;
		}
			
	}

	/* Right */
	for (int i=0; i< right.nactive; i++)
	{	
		change_flag = 1;
		for(int j = 0; j< result[0].nactive; j++) 
		{
			if(result[0].activegrain[j] == right.activegrain[i] ) 
			{
				change_flag = 0;
				break; //Break the current loop and go to next indices in right
			}
		}
		if(change_flag == 1)
		{
			result[0].activegrain[count] = right.activegrain[i];
			result[0].phi[count] = 0.0;
			result[0].phase[count] = right.phase[i];
			result[0].nactive += 1;
			count += 1;
		}	
	}
	/****************************************************************************************************/
	result[0] = node_membersort(result[0]);

	return result[0];
}

node node_membersort(node n)
{
	/*Insertion Sort for small arrays : More efficient*/
	int i = 0, j = 0, key = 0;
	int temp_grain;
	char temp_phase;
	float temp_phi;
	for( i = 1; i < n.nactive; i++)
	{
		key = n.activegrain[i];
		j = i - 1;
		while(j >= 0 && n.activegrain[j] > key)
		{
			temp_grain = n.activegrain[j+1];
			n.activegrain[j+1] = n.activegrain[j];
			n.activegrain[j] = temp_grain;

			temp_phi = n.phi[j+1];
			n.phi[j+1] = n.phi[j];
			n.phi[j] = temp_phi;

			temp_phase = n.phase[j+1];
			n.phase[j+1] = n.phase[j];
			n.phase[j] = temp_phase;
			j = j-1;
		}
	}
	return n;
}

float node_laplacian1D(int grainnum, node left, node mid, node right, float dx)
{
	float val_left, val_mid, val_right;

	val_left = node_phival(left, grainnum);
	val_mid = node_phival(mid, grainnum);
	val_right = node_phival(right, grainnum);

	return ((val_left + val_right - 2.0*val_mid)/(dx*dx));
}

float node_phival(node n, int grainnum)
{
	for(int i = 0; i < n.nactive; i++)
	{
		if(n.activegrain[i] == grainnum) return n.phi[i];
		else if(n.activegrain[i] > grainnum) return 0.0;
	}
	return 0.0;
}


void node_print(node n, int gid)
{
	printf("Grid id: %d\n", gid );
	printf("Number of Grains: %d\n", n.nactive);
	printf("Active Grains: ");
	for(int i = 0; i < n.nactive; i++) printf("%d ", n.activegrain[i]);
	printf("\n");
	printf("Phi: ");
	for(int i = 0; i < n.nactive; i++) printf("%f ", n.phi[i]);
	printf("\n");
	printf("Phases: ");
	for(int i = 0; i < n.nactive; i++) printf("%c ", n.phase[i]);
	printf("\n \n");
	return;
}



void node_outputphi_tofile(char *FILENAME, int grainnum)
{
	extern node *grid;
	extern int gridsize;
	float data[gridsize];
	int flag = 0;
	for(int i=0; i<gridsize ; i++)
	{
		flag = 0;
		for(int j=0; j< grid[i].nactive; j++)
		{
			if (grid[i].activegrain[j] == grainnum)  
				{
					data[i] = grid[i].phi[j];
					flag = 1;
					break;
				}
		}
		if(flag == 0) data[i] = 0.0;
	}

	FILE *fp = fopen(FILENAME, "w");
	for(int i=0; i<gridsize ; i++) fprintf(fp, "%d %f\n", i, data[i]);
	fclose(fp);
}




/*Not actually dependent on node, but solver functions*/
float node_interface(float laplace, float phi)
{
	extern float mobility;
	extern float se;
	extern float dx;
	extern int eta;
	float PI = 3.14159;

	return ((4.0*mobility*se/(eta*dx)) * ( ( eta*dx/PI)*(eta*dx/PI)*laplace  +  (phi) ) );
}

float node_addforce(int grainnum, float phi, float beta)
{
	extern int grainforce;
	if(grainforce == 0) return 0.0;
	else if (grainforce == grainnum ) return (beta*6*phi*(1.0 - phi));
	else return 0;
}

node node_update(node list, int gid)
{
	//Count how many grains have non zero phi
	node copyto;
	copyto.nactive = 0;
	float SIZE_f = sizeof(float);
	int SIZE_i = sizeof(int);
	char SIZE_c = sizeof(char);
	copyto.activegrain = (int *)malloc(5*SIZE_i);
	copyto.phi = (float *)malloc(5*SIZE_f);
	copyto.phase = (char *)malloc(5*SIZE_c);

	for(int i = 0; i< list.nactive; i++) 
	{
		if(list.phi[i] > 0.0 && list.phi[i] <= 1.0 && list.activegrain[i] != 0) 
		{
			copyto.activegrain[copyto.nactive] = list.activegrain[i];

			copyto.phi[copyto.nactive] = list.phi[i];

			copyto.phase[copyto.nactive] = list.phase[i];
			
			copyto.nactive += 1;			
		}
	}
	// if(gid == 498)
	// {	
	// 	printf("After Update\n");
	// 	printf("Nactive: %d\n", list.nactive);
	// 	node_print(list, gid);
	// 	printf("copyto\n");
	// 	node_print(copyto, gid);
	// } 
		 
	return copyto;
}

node node_transfer(node new_list, node list)
{
	/*If unchanged, return the same*/
	if(new_list.nactive == 1 && list.nactive == 1) return list;
	/*Otherwise return the new list*/
	else 
	{	
		float SIZE_f = sizeof(float);
		int SIZE_i = sizeof(int);
		char SIZE_c = sizeof(char);
		int size;
		node_dealloc(list);
		node *list;
		list = node_alloc(1);
		/*Allocate Memory*/
		size = new_list.nactive;
		list[0].activegrain = (int *)malloc(size*SIZE_i);
		list[0].phi = (float *)malloc(size*SIZE_f);
		list[0].phase = (char *)malloc(size*SIZE_c);
		list[0].nactive = new_list.nactive;

		for(int i=0; i<new_list.nactive; i++)
		{
			list[0].activegrain[i] = new_list.activegrain[i];
			list[0].phi[i] = new_list.phi[i];
			list[0].phase[i] = new_list.phase[i];
		}
	return list[0];	
	}

}

node node_simplex(node list)
{
	float sum= 0;
	for(int i = 0; i< list.nactive; i++)
	{
		if(list.phi[i] > 1.0) list.phi[i] = 1.0;
		if(list.phi[i] < 0.0) list.phi[i] = 0.0;
		sum += list.phi[i];
	}
	for(int i = 0;i < list.nactive; i++)
	{
		list.phi[i] /= sum;
	}
	return list;
}

