#include<stdio.h>
#include"nodefunctions.h"
#include<stdlib.h>
#include<string.h>

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
		a->comp = 0;
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
			a[i].comp = 0;
		}
	}
	else printf("Memory not available");
	
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


int node_checkgradient2D(node up, node left, node centre, node right, node down)
{
	/*Check the gradient by checking the number of active grains*/
	/*
	returns 1 : if active grain of any of the node is more than 2
			  : if active grain is 1 and activegrain[0] is different
	*/

	if(up.nactive > 1 || left.nactive > 1 || centre.nactive > 1 || right.nactive > 1 || down.nactive > 1) return 1;
	else if(left.activegrain[0] != centre.activegrain[0] || right.activegrain[0] != centre.activegrain[0] || 
		up.activegrain[0] != centre.activegrain[0] || centre.activegrain[0] != centre.activegrain[0]) return 1;
	else return 0; 
}

int node_checkgradient2D_comp(float c_up, float c_left, float c_centre,  float c_right, float c_down)
{
	/*Check the gradient of composition within the domain*/
	/* 
	returns 1 : if composition has a finite gradient and not equal to 0
			0 : if the gradient is non-zero
	*/
	if(c_up != c_centre || c_left != c_centre || c_right != c_centre || c_down != c_centre) return 1;
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

node node_buildlist2D(node up, node left, node centre, node right, node down)
{
	node result; //Create a node
	result.nactive = 0;
	//result = node_alloc(1); //Allocate memory to the pointer of size 1, list is only one node
	result = _node_checklist(centre, result); //Always pass the centre first, since we need to copy the phi of mid node
	result = _node_checklist(up, result);
	result = _node_checklist(left, result);
	result = _node_checklist(right, result);
	result = _node_checklist(down, result);
	result = node_membersort(result); // Sort the array for fast laplace calculations

	return result;
	 //length of list
	
}

node _node_checklist(node n, node result)
{
	int count = result.nactive;
	int change_flag = 0;
	float SIZE_f = sizeof(float);
	int SIZE_i = sizeof(int);
	char SIZE_c = sizeof(char);
	if(result.nactive == 0)
	{
		result.activegrain = (int *)malloc(5 * SIZE_i); //5 is the maximum limit of neihbourhood grains
		result.phi = (float *)malloc(5 * SIZE_f);
		result.phase = (char *)malloc(5 * SIZE_c);
		//Blindly copy the node
		for(int i=0;i < n.nactive; i++)
		{
			result.activegrain[count] = n.activegrain[i];
			result.phi[count] = n.phi[i];
			result.phase[count] = n.phase[i];
			count += 1;
		}
		result.nactive = count;
		result.comp = n.comp;
	}
	else
	{
		for (int i=0; i< n.nactive; i++)
		{	
			change_flag = 1;
			for(int j = 0; j< result.nactive; j++) 
			{
				if(result.activegrain[j] == n.activegrain[i] ) 
				{
					change_flag = 0;
					break; //Break the current loop and go to next indices in right
				}
			}
			if(change_flag == 1)
			{
				result.activegrain[count] = n.activegrain[i];
				result.phi[count] = 0.0;
				result.phase[count] = n.phase[i];
				result.nactive += 1;
				count += 1;
			}	
		}
	}
	return result;
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

float node_gradient2D(int grainnum, node second, node first, float dx)
{
	float val_second, val_first;
	val_second = node_phival(second, grainnum);
	val_first = node_phival(first, grainnum);

	return (val_second - val_first)/(2.0*dx);
}

float node_laplacian1D(int grainnum, node left, node mid, node right, float dx)
{
	float val_left, val_mid, val_right;

	val_left = node_phival(left, grainnum);
	val_mid = node_phival(mid, grainnum);
	val_right = node_phival(right, grainnum);

	return ((val_left + val_right - 2.0*val_mid)/(dx*dx));
}

float node_laplacian2D(int grainnum, node up, node left, node centre, node right, node down, float dx)
{
	float val_up, val_left, val_centre, val_right, val_down;

	val_up = node_phival(up, grainnum);
	val_left = node_phival(left, grainnum);
	val_centre = node_phival(centre, grainnum);
	val_right = node_phival(right, grainnum);
	val_down = node_phival(down , grainnum);

	return ((val_left + val_right - 2*val_centre)/(dx*dx) + (val_up + val_down - 2*val_centre)/(dx*dx));
}

float node_laplacian2D_comp(float c_up, float c_left, float c_centre, float c_right, float c_down, float dx)
{
	return ((c_up + c_down - 2.0*c_centre)/(dx*dx) + (c_right + c_left - 2.0*c_centre)/(dx*dx));
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

void node_handleBC(int gx, int gy, node *new, char *bc)
{
	int copyfrom, copyto;
	int index_copyfrom, index_copyto;
	char impbc[] = "drichlet";
	if(!strcmp(bc,"drichlet"))
	{
		//First move in x direction
		for(int i = 1; i < gx-1; i++)
		{
			copyfrom = 1;
			copyto = 0;
			index_copyfrom = copyfrom*gx + i;
			index_copyto = copyto*gx + i;
			new[index_copyto] = node_update(new[index_copyfrom], index_copyto); 
			copyfrom = gy - 2;
			copyto = gy - 1;
			index_copyfrom = copyfrom*gx + i;
			index_copyto = copyto*gx + i;
			new[index_copyto] = node_update(new[index_copyfrom], index_copyto);
		}
		// Handle the y direction
		for(int i = 0; i < gy; i++)
		{
			copyfrom = 1;
			copyto = 0;
			index_copyfrom = i*gx + copyfrom;
			index_copyto = i*gx + copyto;
			new[index_copyto] = node_update(new[index_copyfrom], index_copyto); 
			copyfrom = gx - 2;
			copyto = gx - 1;
			index_copyfrom = i*gx + copyfrom;
			index_copyto = i*gx + copyto;
			new[index_copyto] = node_update(new[index_copyfrom], index_copyto);
		}
	}
	else 
	{
		printf("CHange the boundary conditon, this has not been implemented yet\n");
		return;
	}
}


void node_print(node n, int gid)
{
	printf("Grid id: %d\n", gid );
	printf("Number of Grains: %d\n", n.nactive);
	printf("Composition: %f\n", n.comp);
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

void node_outputcomp_tofile(char *FILENAME)
{
	extern node *grid;
	extern int gridsize;
	float c[gridsize];
	int flag = 0;
	for(int i=0; i<gridsize ; i++) c[i] = grid[i].comp;

	FILE *fp = fopen(FILENAME, "w");
	for(int i=0; i<gridsize ; i++) fprintf(fp, "%d %f\n", i, c[i]);
	fclose(fp);
}

void node_visual(char *FILENAME, int grainnum)
{	
	//printf("Inside Visual Node\n");
	extern node *grid;
	extern int gridx, gridy;
	extern int gridsize;

	float data[gridsize];
	int flag = 0;
	int index;

	for(int i=0; i<gridsize; i++) data[i] = 0;
	FILE *fp = fopen(FILENAME, "w");
	for(int i = 0; i < gridx; i++)
	{
		for(int j = 0; j < gridy; j++)
		{
			flag = 0;
			index = j*gridx + i;
			for(int k = 0; k < grid[index].nactive; k++)
			{
				if(grid[index].activegrain[k] == grainnum) 
				{
					data[index] = grid[index].phi[k];
					flag = 1;
					break;
				}
			}
			if(flag == 0) data[index] = 0.0;
		}
	}

	//output
	for(int i=0; i<gridx; i++)
	{
		for(int j = 0; j<gridy; j++)
		{
			index = j*gridx + i;
			fprintf(fp, "%d %d %f\n", i, j, data[index]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}




/*Not actually dependent on node, but solver functions*/
float node_interface(float laplace, float phi)
{
	extern float nd_mobility;
	extern float nd_se;
	extern float nd_dx;
	extern int eta;
	float PI = 3.14159;

	return ((4.0*nd_mobility*nd_se/(eta*nd_dx)) * ( ( eta*nd_dx/PI)*(eta*nd_dx/PI)*laplace  +  (phi) ) );
}

float node_addforce(int grainnum, float phi, float beta, int t)
{
	extern int grainforce;
	extern int beta_timestep;
	if(grainforce == 0 || t < beta_timestep) return 0.0;
	else if (grainforce == grainnum ) return (-beta*6*phi*(1.0 - phi));
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
	int max_n = list.nactive;
	copyto.activegrain = (int *)malloc(max_n*SIZE_i);
	copyto.phi = (float *)malloc(max_n*SIZE_f);
	copyto.phase = (char *)malloc(max_n*SIZE_c);

	for(int i = 0; i< max_n; i++) 
	{
		if(list.phi[i] > 0.0 && list.phi[i] <= 1.0 && list.activegrain[i] != 0) 
		{
			copyto.activegrain[copyto.nactive] = list.activegrain[i];

			copyto.phi[copyto.nactive] = list.phi[i];

			copyto.phase[copyto.nactive] = list.phase[i];
			
			copyto.nactive += 1;			
		}
	}
	copyto.comp = list.comp;
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
	if(new_list.nactive == 1 && list.nactive == 1 && list.comp == new_list.comp) return list;
	/*Otherwise return the new list*/
	else 
	{	
		float SIZE_f = sizeof(float);
		int SIZE_i = sizeof(int);
		char SIZE_c = sizeof(char);
		int size;
		//node_dealloc(list);
		//node *list;
		//list = node_alloc(1);
		/*Allocate Memory*/
		size = new_list.nactive;
		//list.activegrain = (int *)malloc(size*SIZE_i);
		//list.phi = (float *)malloc(size*SIZE_f);
		//list.phase = (char *)malloc(size*SIZE_c);
		list.nactive = new_list.nactive;
		list.comp = new_list.comp;

		for(int i=0; i<new_list.nactive; i++)
		{
			list.activegrain[i] = new_list.activegrain[i];
			list.phi[i] = new_list.phi[i];
			list.phase[i] = new_list.phase[i];
		}
	return list;	
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

node node_phisolver(char *solver, int index, node up, node left, node centre, node right, node down, node update, float dx, int t)
{
	float *laplacian;
	float *term;
	float PI = 3.14159;

	int intf_flag = 1;

	node list;
	extern float omega;
	extern float epsilon, mob_phi;
	extern float nd_dt;
	extern float alpha;
	extern int grainforce;
	extern float nd_beta;

	if(!strcmp(solver,"MPF"))
	{
		intf_flag = node_checkgradient2D(up, left, centre, right, down);
		if(!intf_flag)
		{
			update = node_update(centre, index);
			return update;
		}
		list = node_buildlist2D(up, left, centre, right, down);
		//printf("index : %d\n", index);
		//Allocate memory and initialize laplacian and overall derivative
		laplacian = (float *)malloc(list.nactive*sizeof(float));
		for(int temp = 0; temp < list.nactive; temp++) laplacian[temp] = 0.0;

		term = (float *)malloc(list.nactive*sizeof(float));
		for(int temp = 0; temp < list.nactive; temp++) term[temp] = 0.0;

		//Build Laplacians
		for(int j = 0; j < list.nactive; j++)
			laplacian[j] = node_laplacian2D(list.activegrain[j], up, left, centre, right, down, dx);
		//for(int j = 0; j < list.nactive; j++) printf("Laplacian : %d Grain %d Phi:%f laplace : %f\n", j, list.activegrain[j], list.phi[j], laplacian[j]);
		for(int j = 0; j< list.nactive; j++)
		{

			for(int k = 0; k < list.nactive; k++)
			{
				term[j] += epsilon*epsilon/2.0*(laplacian[j] - laplacian[k]) 
				+ omega*(1 - alpha*list.comp)*(list.phi[j] - list.phi[k]) 
				+ node_addforce(list.activegrain[j], list.phi[j], nd_beta, t) - node_addforce(list.activegrain[k], list.phi[k], nd_beta, t); //Add the chemical energy term here
			}
		}
		//for(int j = 0; j < list.nactive; j++) printf("term = %d Grain %d Phi:%f laplace : %f, term : %f\n", j, list.activegrain[j], list.phi[j], laplacian[j], term[j]);
		for (int j = 0; j < list.nactive ; j++)
		{
			term[j] *= 2.0*mob_phi/list.nactive;
			term[j] = nd_dt*term[j] + list.phi[j];
			list.phi[j] = term[j];
		}
		list = node_simplex(list); //Limit to zero and one
		update = node_update(list, index);
		free(term);
		//node_print(list, index);
		free(laplacian);
		node_dealloc(list);
		return update;
	}
}

node node_csolver(char *solver, int index, node up, node left, node centre, node right, node down, node update, float dx, int t)
{
	
	extern float nd_diffusivity, omega;
	extern float nd_dt;

	if(!strcmp(solver, "equalcomp"))
	{
		int c_flag, intf_flag;
		float laplacian = 0.0;
		float term = 0.0;
		float c = centre.comp;
		float K = 1.0; //Should be taken out into input/initialize for scaling the thickness
		c_flag = node_checkgradient2D_comp(up.comp, left.comp, centre.comp, right.comp, down.comp);
		intf_flag = node_checkgradient2D(up, left, centre, right, down);
		if(!c_flag && !intf_flag)
		{
			//Composition is already updated in phi_solver, so no need to change here
			return update;
		}

		
		laplacian = node_laplacian2D_comp(up.comp, left.comp, centre.comp, right.comp, down.comp, dx);
		//K is scaling factor used by Kim
		// Assuming diffusivity is constant everywhere
		term = nd_diffusivity*laplacian - (nd_diffusivity*K)*(node_termsumphi_comp(up, left, centre, right, down, dx));
		if(index == 500 && t%10000 == 0) printf("Index : %d, term phi = %f, term lap = %f, diffusivity = %f\n", index, node_termsumphi_comp(up, left, centre, right, down, dx), laplacian, nd_diffusivity ); 
		update.comp = centre.comp + nd_dt*(term);
		if(index == 500 && t%10000 == 0) printf("Extra = %e, Composition old= %e, new= %e, timestep = %f new_comp = %e\n", nd_dt*(term), centre.comp, update.comp, nd_dt, centre.comp + nd_dt*(term) );
		//printf("index : %d, comp : %f\n", index, update.comp);
		return update;
	}
	else if(!strcmp(solver, "equaldiff"))
	{
		printf("%s - select an appropriate csolver keyword\n", solver);
		exit(5);
	}
	else if(!strcmp(solver, "None"))
	{
		return update; //Dont update the index and return as is
	}
	else 
	{
		printf("Choose an appropriate condition for composition solver");
		exit(5);
	}
}

float node_termsumphi_comp(node up, node left, node centre, node right, node down, float dx)
{
	extern float alpha;
	extern float omega;
	float c = centre.comp;
	/*Don't calculate anything if not required*/
	if (alpha <= 0 || omega <= 0) return 0.0;

	float term_laplace, term_grad;
	

	//Using central difference scheme for gradient calculation
	float grad_cx = (right.comp - left.comp)/(2.0*dx);
	float grad_cy = (up.comp - down.comp)/(2.0*dx);

	//COMPUTING LAPLACIAN AGAIN - can be merged with solver_phi to not to recompute again
	/**************************************************************************/
	node list;
	list = node_buildlist2D(up, left, centre, right, down);
	float laplacian[list.nactive];
	float gradx[list.nactive], grady[list.nactive];
	for(int j=0; j<list.nactive; j++)
	{
		laplacian[j] = node_laplacian2D(list.activegrain[j], up, left, centre, right, down, dx);
		gradx[j] = node_gradient2D(list.activegrain[j], right, left, dx);
		grady[j] = node_gradient2D(list.activegrain[j], up, down, dx);
	}
	/***************************************************************************/
	float gradterm = 0, lapterm = 0, sqgrad = 0;
	for (int j=0; j<list.nactive; j++)
	{
		gradterm += omega*alpha*(1.0 - 2.0*c)*(grad_cx*gradx[j] + grad_cy*grady[j])*(1.0 - list.phi[j]);
		lapterm  += omega*alpha*(c)*(1.0 - c)*(1.0 - list.phi[j])*laplacian[j];
		sqgrad   += omega*alpha*(c)*(1.0 - c)*(gradx[j]*gradx[j] + grady[j]*grady[j]);
	}
	return (gradterm + lapterm - sqgrad);



//	float prodphi[5]; 
	/***********************/
	//        up[0]
	//left[1] centre[2] right[3]
	//        down[4]
	/***********************/
	/*
	prodphi[0] = node_prodphi_comp(up);
	prodphi[1] = node_prodphi_comp(left);
	prodphi[2] = node_prodphi_comp(centre);
	prodphi[3] = node_prodphi_comp(right);
	prodphi[4] = node_prodphi_comp(down);

	float grad_prodphix = (prodphi[3] - prodphi[1])/(2.0*dx);
	float grad_prodphiy = (prodphi[0] - prodphi[4])/(2.0*dx);

	float laplace_prodphi = (prodphi[0] + prodphi[4] - 2.0*prodphi[2])/(dx*dx) + (prodphi[1] + prodphi[3] - 2.0*prodphi[2])/(dx*dx);

	term_laplace = c*(1-c)*omega*alpha*laplace_prodphi;
	term_grad = (1 - 2*c)*omega*alpha*(grad_cx*grad_prodphix + grad_cy*grad_prodphiy);

	return (term_grad + term_laplace);
	*/
}

float node_prodphi_comp(node n)
{
	float sum = 0;
	for (int i=0; i<n.nactive; i++)
	{
		for(int j=i; j<n.nactive; j++)
		{
			sum += n.phi[i]*n.phi[j];
		}
	}
	return sum;
}