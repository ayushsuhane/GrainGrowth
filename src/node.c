#include<stdio.h>
#include"nodefunctions.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>

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
		a->activegrain[0] = 0;
		a->phi[0] = 0;
		a->phase[0] = 0;
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
			a[i].activegrain[0] = 0;
			a[i].phi[0] = 0;
			a[i].phase[0] = 0;
			a[i].comp = 0;
		}
	}
	else printf("Memory not available");
	
	return a;
}

//
//void node_dealloc(node  free_node)
//{
	/*
	Takes in a pointer to node array 
	to deallocate the memory
	*/
	/* First free the memory allocated by pointers of node -  active_grain, phi, phase.
	Then clear the memory occupied by node. Generally preferred to go in a reverse order, otherwise 
	causes heap memory leak
	*/
///		free(free_node.phase);
///		free(free_node.phi);
///		free(free_node.activegrain);
		//Check Whether to free the free_node as well??????
	//free(free_node);

//	return;
//}
//

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
	else if((left.activegrain[0] != centre.activegrain[0]) || (right.activegrain[0] != centre.activegrain[0]) || 
		(up.activegrain[0] != centre.activegrain[0]) || (down.activegrain[0] != centre.activegrain[0])) return 1;
	else return 0; 
}

int node_checkgradient2D_comp(double c_up, double c_left, double c_centre,  double c_right, double c_down)
{
	/*Check the gradient of composition within the domain*/
	/* 
	returns 1 : if composition has a finite gradient and not equal to 0
			0 : if the gradient is non-zero
	*/
	double c_tol = 1e-15;
	if( (fabs(c_left-c_centre) > c_tol) || (fabs(c_right-c_centre) > c_tol) || (fabs(c_down-c_centre) > c_tol) || (fabs(c_up-c_centre) > c_tol)) return 1;
	else return 0;
}

/*Not in USE*/
node  node_buildlist(node left, node mid, node right)
{
	/*
	Takes in neighbouring nodes, with sorted lists, and append the list to be checked in the result list
	*/
	// Allocate memory for result
	node *result;

	//double SIZE_f = sizeof(double);
	//int SIZE_i = sizeof(int);
	//char SIZE_c = sizeof(char);

	int count = 0, change_flag = 0; //length of list
	result = node_alloc(1);

	//Start by copying the entries in mid node
	// Allocate size to copy the list
	//result[0].activegrain = (int *)malloc(5 * SIZE_i); //5 is the maximum limit of neihbourhood grains
	//result[0].phi = (double *)malloc(5 * SIZE_f);
	//result[0].phase = (char *)malloc(5 * SIZE_c);
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

//Have to free the pointer whenever used
node _node_checklist(node n, node result)
{
	int count = result.nactive;
	int change_flag = 0;
	//double SIZE_f = sizeof(double);
	//int SIZE_i = sizeof(int);
	//char SIZE_c = sizeof(char);
	if(result.nactive == 0)
	{
		//node new;
		//new.activegrain = (int *)malloc(5 * SIZE_i); //5 is the maximum limit of neihbourhood grains
		//new.phi = (double *)malloc(5 * SIZE_f);
		//new.phase = (char *)malloc(5 * SIZE_c);
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
		//result = node_update(new, 0);
		//node_dealloc(new);
	}
	else
	{
		for (int i=0; i< n.nactive; i++)
		{	
			change_flag = 1;
			for(int j = 0; j< result.nactive; j++) 
			{
				if(n.activegrain[i] == result.activegrain[j]) 
				{
					change_flag = 0;
					break; //Break the current loop and go to next indices
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
	double temp_phi;
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

double node_gradient2D(int grainnum, node second, node first, double dx)
{
	double val_second, val_first;
	val_second = node_phival(second, grainnum);
	val_first = node_phival(first, grainnum);

	return ((val_second - val_first)/(2.0*dx));
	//return ((val_second - val_first)/(dx));
}

double node_laplacian1D(int grainnum, node left, node mid, node right, double dx)
{
	double val_left, val_mid, val_right;

	val_left = node_phival(left, grainnum);
	val_mid = node_phival(mid, grainnum);
	val_right = node_phival(right, grainnum);

	return ((val_left + val_right - 2.0*val_mid)/(dx*dx));
}

double node_laplacian2D(int grainnum, node up, node left, node centre, node right, node down, double dx)
{
	double val_up, val_left, val_centre, val_right, val_down;

	val_up = node_phival(up, grainnum);
	val_left = node_phival(left, grainnum);
	val_centre = node_phival(centre, grainnum);
	val_right = node_phival(right, grainnum);
	val_down = node_phival(down , grainnum);

	return ((val_left + val_right - 2.0*val_centre)/(dx*dx) + (val_up + val_down - 2.0*val_centre)/(dx*dx));
}

double node_laplacian2D_comp(double c_up, double c_left, double c_centre, double c_right, double c_down, double dx)
{
	return ((c_up + c_down - 2.0*c_centre)/(dx*dx) + (c_right + c_left - 2.0*c_centre)/(dx*dx));
}

double node_phival(node n, int grainnum)
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
		for(int i = 0; i < gy ; i++)
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
	printf("Composition: %e\n", n.comp);
	printf("Active Grains: ");
	for(int i = 0; i < n.nactive; i++) printf("%d ", n.activegrain[i]);
	printf("\n");
	printf("Phi: ");
	for(int i = 0; i < n.nactive; i++) printf("%e ", n.phi[i]);
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
	extern int gridx, gridy;
	double data[gridsize];
	int flag = 0, index;
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
	for(int i=0; i<gridx; i++)
	{
		for(int j=0; j<gridy; j++)
		{
			index = j*gridx + i;
			fprintf(fp, "%d %d %e\n", i, j, data[index]);
		}
		if(gridy>1) fprintf(fp, "\n");
	}
	fclose(fp);
}

void node_outputcomp_tofile(char *FILENAME)
{
	extern node *grid;
	extern int gridsize;
	extern int gridx, gridy;
	double c[gridsize];
	int index;
	for(int i=0; i<gridsize ; i++) c[i] = grid[i].comp;

	FILE *fp = fopen(FILENAME, "w");
	for(int i=0; i<gridx; i++)
	{
		for(int j=0; j<gridy; j++)
		{
			index = j*gridx + i;
			fprintf(fp, "%d %d %e\n", i, j, c[index]);
		}
		if(gridy>1) fprintf(fp, "\n");
	}
	//for(int i=0; i<gridsize ; i++) fprintf(fp, "%d %e\n", i, c[i]);
	fclose(fp);
}

void node_visual(char *FILENAME, int grainnum)
{	
	//printf("Inside Visual Node\n");
	extern node *grid;
	extern int gridx, gridy;
	extern int gridsize;

	double data[gridsize];
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
			fprintf(fp, "%d %d %e\n", i, j, data[index]);
		}
		if(gridy > 1)fprintf(fp,"\n");
	}
	fclose(fp);
}




/*Not actually dependent on node, but solver functions*/
double node_interface(double laplace, double phi)
{
	extern double nd_mobility;
	extern double nd_se;
	extern double nd_dx;
	extern int eta;
	double PI = 3.14159;

	return ((4.0*nd_mobility*nd_se/(eta*nd_dx)) * ( ( eta*nd_dx/PI)*(eta*nd_dx/PI)*laplace  +  (phi) ) );
}

double node_addforce(int grainnum, double phi, double beta, int t)
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
	//double SIZE_f = sizeof(double);
	//int SIZE_i = sizeof(int);
	//char SIZE_c = sizeof(char);
	int max_n = list.nactive;
	double tol = 1e-15;
	//copyto.activegrain = (int *)malloc(max_n*SIZE_i);
	//copyto.phi = (double *)malloc(max_n*SIZE_f);
	//copyto.phase = (char *)malloc(max_n*SIZE_c);

	for(int i = 0; i< max_n; i++) 
	{
		if((list.phi[i] > tol) && ((list.phi[i]) <= 1.0) && list.activegrain[i] != 0) 
		{
			copyto.activegrain[copyto.nactive] = list.activegrain[i];

			copyto.phi[copyto.nactive] = list.phi[i];

			copyto.phase[copyto.nactive] = list.phase[i];
			
			copyto.nactive += 1;			
		}
	}
	copyto.comp = list.comp;
	//node new;
	//node_transfer(new, copyto);
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
	double tol = 1e-15;
	if((new_list.nactive == 1) && (list.nactive == 1) && (new_list.activegrain[0] == list.activegrain[0]) && (fabs(list.comp - new_list.comp) < tol) ) return list;
	/*Otherwise return the new list*/	
	else 
	{	
		int size;
		//node_dealloc(list);
		//node *list;
		//list = node_alloc(1);
		/*Allocate Memory*/
		size = new_list.nactive;
		//list.activegrain = (int *)malloc(size*SIZE_i);
		//list.phi = (double *)malloc(size*SIZE_f);
		//list.phase = (char *)malloc(size*SIZE_c);
		list.nactive = new_list.nactive;
		list.comp = new_list.comp;

		for(int i=0; i<size; i++)
		{
			list.activegrain[i] = new_list.activegrain[i];
			list.phi[i] = new_list.phi[i];
			list.phase[i] = new_list.phase[i];
		}
		return list;	
	}
	return list;
}

node node_simplex(node list)
{
	double sum= 0;
	for(int i = 0; i< list.nactive; i++)
	{
		if(list.phi[i] > 1.0) list.phi[i] = 1.0;
		if(list.phi[i] < 0.0) list.phi[i] = 0.0;
		sum += list.phi[i];
	}
	for(int i = 0;i < list.nactive; i++)
	{
		list.phi[i] /= sum;
		//list.phi[i] = list.phi[i] + (1 - sum)/list.nactive;
	}
	return list;
}

node node_phisolver(char *phisolver, int index, node up, node left, node centre, node right, node down, node update, double dx, int t)
{
	double *laplacian;
	double *term;

	int intf_flag = 1;

	
	extern double omega;
	extern double epsilon, mob_phi;
	extern double nd_dt;
	extern double alpha;
	extern int grainforce;
	extern double nd_beta;
	//double nd_dt = dt;
	if(!strcmp(phisolver,"MPF"))
	{
		node list;
		intf_flag = node_checkgradient2D(up, left, centre, right, down);
		if(!intf_flag)
		{
			update = node_update(centre, index);
			return update;
		}
		list = node_buildlist2D(up, left, centre, right, down);
		//printf("index : %d\n", index);
		//Allocate memory and initialize laplacian and overall derivative
		laplacian = (double *)malloc(list.nactive*sizeof(double));
		for(int temp = 0; temp < list.nactive; temp++) laplacian[temp] = 0.0;

		term = (double *)malloc(list.nactive*sizeof(double));
		for(int temp = 0; temp < list.nactive; temp++) term[temp] = 0.0;

		//Build Laplacians
		for(int j = 0; j < list.nactive; j++)
			laplacian[j] = node_laplacian2D(list.activegrain[j], up, left, centre, right, down, dx);
		//for(int j = 0; j < list.nactive; j++) printf("Laplacian : %d Grain %d Phi:%e laplace : %e\n", j, list.activegrain[j], list.phi[j], laplacian[j]);
		for(int j = 0; j<list.nactive; j++)
		{
			for(int k = 0; k < list.nactive; k++)
			{
				term[j] += (epsilon*epsilon/2.0)*(laplacian[j] - laplacian[k]) 
				+ omega*(1.0 - alpha*centre.comp)*(list.phi[j] - list.phi[k]) 
				+ node_addforce(list.activegrain[j], list.phi[j], nd_beta, t) - node_addforce(list.activegrain[k], list.phi[k], nd_beta, t); //Add the chemical energy term here
			}
		}
		//for(int j = 0; j < list.nactive; j++) printf("term = %d Grain %d Phi:%e laplace : %e, term : %e\n", j, list.activegrain[j], list.phi[j], laplacian[j], term[j]);
		for (int j = 0; j < list.nactive ; j++)
		{
			term[j] *= 2.0*mob_phi/list.nactive;
			term[j] = nd_dt*term[j] + list.phi[j];
			list.phi[j] = term[j];
		}
		list = node_simplex(list); //Limit to zero and one

		///update = node_update(list, index);
		free(term);
		//node_print(list, index);
		free(laplacian);
		///node_dealloc(list);
		
		return list;
	}
	return update;
}

double node_csolver(char *solver, int index, node up, node left, node centre, node right, node down, node update, double dx, int t)
{
	
	extern double nd_diffusivity, omega;
	extern double nd_dt;
	extern int beta_timestep;
	//double nd_dt = dt;

	if(!strcmp(solver, "equalcomp"))
	{
		int c_flag, intf_flag = 0;
		double laplacian = 0.0;
		double term = 0.0;
		double term_sumphi_comp = 0.0;
		double c;

		double K = 1.0; //Should be taken out into input/initialize for scaling the thickness
		c_flag = node_checkgradient2D_comp(up.comp, left.comp, centre.comp, right.comp, down.comp);
		if(update.nactive > 1) intf_flag = 1;
		//if(!c_flag && !intf_flag)
		/*
		if (t < beta_timestep)
		{
			if(intf_flag) c_flag = 1;
			else c_flag = 0;
		}
		*/
		if(!intf_flag && !c_flag)
		{
			//Composition is already updated in phi_solver, so no need to change here
			return centre.comp;
		}
		//Caculate second term only if gradient of phi exist
		//if((index == 501 || index == 500|| index == 509 || index == 508 || index == 507 || index == 506 || index == 504 || index == 505 || index == 503) && (t%10000 == 0 || t > 200000)  ) printf("Index : %d\n", index);
		//if(intf_flag) term_sumphi_comp = node_termsumphi_comp(up, left, centre, right, down, update, dx, t, index);
		term_sumphi_comp = node_termsumphi_comp(up, left, centre, right, down, update, dx, t, index);
		//if(t%100000 == 0) node_print(centre, index);
		/*
		if(t%10000 == 0) 
			{
				node_print(centre, index);
				printf("phi:%d c:%d", intf_flag, c_flag);
			}
		*/
		laplacian = node_laplacian2D_comp(up.comp, left.comp, centre.comp, right.comp, down.comp, dx);
		//K is scaling factor used by Kim
		// Assuming diffusivity is constant everywhere
		term = nd_diffusivity*laplacian - (nd_diffusivity*K)*(term_sumphi_comp);
		if((index == 101 || index == 100|| index == 109 || index == 108 || index == 107 || index == 106 || index == 104 || index == 105 || index == 103) && (t%10000 == 0 || t > 400000000)  ) printf("phi_1 = %e, grain = %d, term lap = %e, diffusivity = %e\n", centre.phi[0], centre.activegrain[0], laplacian, nd_diffusivity ); 
		c = centre.comp + nd_dt*(term);
		if((index == 101 || index == 100|| index == 109 || index == 108 || index == 107 || index == 106 || index == 104 || index == 105 || index == 103) && (t%10000 == 0 || t > 400000000)) printf("Extra = %e, r_comp = %e, l_comp=%e, Composition old= %e, new= %e, timestep = %e new_comp = %e flux = %e\n", nd_dt*(term), right.comp, left.comp, centre.comp, update.comp, nd_dt, centre.comp + nd_dt*(term), node_calcflux_x(up, left, centre, right, down) );
		//printf("index : %d, comp : %e\n", index, update.comp);
		//printf("\n");
		return c;
	}
	else if(!strcmp(solver, "equaldiff"))
	{
		printf("%s - select an appropriate csolver keyword\n", solver);
		exit(5);
	}
	else if(!strcmp(solver, "None"))
	{
		return centre.comp; //Dont update the index and return as is
	}
	else 
	{
		printf("Choose an appropriate condition for composition solver");
		exit(5);
	}
}

double node_termsumphi_comp(node up, node left, node centre, node right, node down, node list, double dx, int t, int index)
{
	extern double alpha;
	extern double omega;
	double c = centre.comp;
	double tol = 1e-15;
	/*Don't calculate anything if not required*/
	if (alpha <= tol || omega <= tol) return 0.0;
	

	//Using central difference scheme for gradient calculation
	//double grad_cx = (right.comp - centre.comp)/(dx);
	//double grad_cy = (up.comp - centre.comp)/(dx);

	double grad_cx = (right.comp - left.comp)/(2.0*dx);
	double grad_cy = (up.comp - down.comp)/(2.0*dx);

	//COMPUTING LAPLACIAN AGAIN - can be merged with solver_phi to not to recompute again
	/**************************************************************************/
	//list = node_buildlist2D(up, left, centre, right, down);
	double laplacian[list.nactive];
	double gradx[list.nactive], grady[list.nactive];
	double sumphi=0;
	for(int j=0; j<list.nactive; j++)
	{
		laplacian[j] = node_laplacian2D(list.activegrain[j], up, left, centre, right, down, dx);
		//gradx[j] = node_gradient2D(list.activegrain[j], right, centre, dx);
		//grady[j] = node_gradient2D(list.activegrain[j], up, centre, dx);
		gradx[j] = node_gradient2D(list.activegrain[j], right, left, dx);
		grady[j] = node_gradient2D(list.activegrain[j], up, down, dx);
	}
	/***************************************************************************/
	//list = node_update(list, index);
	double gradterm = 0.0, lapterm = 0.0, sqgrad = 0.0;
	for (int j=0; j<list.nactive; j++)
	{
		gradterm += omega*alpha*(1.0 - 2.0*c)*(grad_cx*gradx[j] + grad_cy*grady[j])*(1.0 - list.phi[j]);
		lapterm  += omega*alpha*(c)*(1.0 - c)*(1.0 - list.phi[j])*laplacian[j];
		sqgrad   += omega*alpha*(c)*(1.0 - c)*(gradx[j]*gradx[j] + grady[j]*grady[j]);
		if((index == 101 || index == 100|| index == 109 || index == 108 || index == 107 || index == 106 || index == 104 || index == 105 || index == 103) && (t%10000 == 0 || t > 400000000)) 
			printf("index = %d, grain = %d, phi = %e, grad :%e, gradc :%e, lap=%e, t1 = %e, t2=%e, t3=%e cumul=%e\n", index, list.activegrain[j], list.phi[j], gradx[j], grad_cx,laplacian[j], omega*alpha*(1.0 - 2.0*c)*(grad_cx*gradx[j] + grad_cy*grady[j])*(1.0 - list.phi[j]), omega*alpha*(c)*(1.0 - c)*(1.0 - list.phi[j])*laplacian[j], omega*alpha*(c)*(1.0 - c)*(gradx[j]*gradx[j] + grady[j]*grady[j]), gradterm+lapterm-sqgrad);
	}
	if((index == 101 || index == 100|| index == 109 || index == 108 || index == 107 || index == 106 || index == 104 || index == 105 || index == 103) && (t%10000 == 0 || t > 400000000)) printf("term=%e\n", gradterm + lapterm - sqgrad);
	
	//node_dealloc(list);
	//return (gradterm + lapterm - sqgrad);
	return ( gradterm + lapterm - sqgrad);



//	double prodphi[5]; 
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

	double grad_prodphix = (prodphi[3] - prodphi[1])/(2.0*dx);
	double grad_prodphiy = (prodphi[0] - prodphi[4])/(2.0*dx);

	double laplace_prodphi = (prodphi[0] + prodphi[4] - 2.0*prodphi[2])/(dx*dx) + (prodphi[1] + prodphi[3] - 2.0*prodphi[2])/(dx*dx);

	term_laplace = c*(1-c)*omega*alpha*laplace_prodphi;
	term_grad = (1 - 2*c)*omega*alpha*(grad_cx*grad_prodphix + grad_cy*grad_prodphiy);

	return (term_grad + term_laplace);
	*/
}

double node_prodphi_comp(node n)
{
	double sum = 0;
	for (int i=0; i<n.nactive; i++)
	{
		for(int j=i; j<n.nactive; j++)
		{
			sum += n.phi[i]*n.phi[j];
		}
	}
	return sum;
}

double node_calcflux_x(node up, node left, node centre, node right, node down)
{
	double gradcx, gradcy, gradphix;
	extern double nd_dx, omega, alpha;
	gradcx = (right.comp - left.comp)/(2.0*nd_dx);
	node list;
	list = node_buildlist2D(up, left, centre, right, down);
	double sum = 0;
	for(int i=0; i<list.nactive; i++)
	{
		gradphix = (node_phival(right, list.activegrain[i]) - node_phival(left, list.activegrain[i]))/(2.0*nd_dx);
		sum = (1 - list.phi[i])*(gradphix);
	}
	return (gradcx - omega*alpha*centre.comp*(1.0 - centre.comp)*sum);
}




/************************Single***********************************/
int node_checkgradient2D_single(node up, node nextleft, node left, node centre, node right, node nextright, node down)
{
	/*Check the gradient by checking the number of active grains*/
	/*
	returns 1 : if active grain of any of the node is more than 2
			  : if active grain is 1 and activegrain[0] is different
	*/

	double tol = 1e-12;
	if((fabs(left.phi[0] - centre.phi[0]) > tol) ||
		(fabs(right.phi[0] - centre.phi[0]) > tol) ||
		(fabs(up.phi[0] - centre.phi[0]) > tol) ||
		(fabs(down.phi[0] - centre.phi[0]) > tol) ||
		(fabs(nextleft.phi[0] - centre.phi[0]) > tol)||
		(fabs(nextright.phi[0] - centre.phi[0]) > tol)) return 1;
	else return 0; 
}

int node_checkgradient2D_comp_single(double c_up, double c_nextleft, double c_left, double c_centre, double c_right, double c_nextright, double c_down)
{
	/*Check the gradient by checking the number of active grains*/
	/*
	returns 1 : if active grain of any of the node is more than 2
			  : if active grain is 1 and activegrain[0] is different
	*/

	double tol = 1e-12;
	if((fabs(c_left - c_centre) > tol) ||
		(fabs(c_right- c_centre) > tol) ||
		(fabs(c_up - c_centre) > tol) ||
		(fabs(c_down - c_centre) > tol) ||
		(fabs(c_nextleft - c_centre) > tol) ||
		(fabs(c_nextright - c_centre) > tol)) return 1;
	else return 0; 
}


node node_update_single(node n)
{
	return n;
}
node node_phisolver_single(char *phisolver, int index, node up, node nextleft, node left, node centre, node right, node nextright, node down, node update, double dx, int t)
{
	
	int intf_flag = 1;

	
	extern double omega;
	extern double epsilon, mob_phi;
	extern double nd_dt;
	extern double alpha;
	extern int grainforce;
	extern double nd_beta;
	extern int beta_timestep;

	if(!strcmp(phisolver,"MPF"))
	{
		node list;
		list.nactive = centre.nactive;
		list.comp = centre.comp;
		for(int i=0; i<list.nactive; i++) 
		{
			list.activegrain[i] = centre.activegrain[i];
			list.phase[i] = centre.phase[i];
		}
		intf_flag = node_checkgradient2D_single(up, nextleft, left, centre, right, nextright, down);
		if(!intf_flag) return centre;
		double laplacian;
		double k = 0.0;
		if(t > beta_timestep) k = 1.0;
		laplacian = ((right.phi[0] + left.phi[0] - 2.0*centre.phi[0])/(dx*dx) + (up.phi[0] + down.phi[0] - 2.0*centre.phi[0])/(dx*dx));
		//laplacian = (-nextright.phi[0] + 16.0*right.phi[0] - 30.0*centre.phi[0] + 16.0*left.phi[0] - nextleft.phi[0])/(12.0*dx*dx);
		for (int i=0; i< list.nactive; i++)
		{ 
			list.phi[i] = centre.phi[i] + nd_dt*mob_phi*(epsilon*epsilon*laplacian - omega*(1.0 - alpha*centre.comp)*(1.0 - 2.0*centre.phi[i]) - k*nd_beta*6.0*centre.phi[0]*(1.0-centre.phi[0]));
			if(list.phi[i] < 0.0) list.phi[i] = 0.0;
			if(list.phi[i] > 1.0) list.phi[i] = 1.0;
		}
		return list;
	}
}

double node_csolver_single(char *solver, int index, node up, node nextleft, node left, node centre, node right, node nextright, node down, double dx, int t)
{
	if(!strcmp(solver, "equalcomp"))
	{
		double c;
		double gradphi, gradc;
		double laplacian_phi, laplacian_c;
		extern double nd_dt;
		extern double nd_diffusivity;
		extern double omega, alpha;
		extern int phi_timestep;
		double tol = 1e-12;
		if(alpha < tol) return centre.comp;
		int intf_flag;
		int c_flag;
		intf_flag = node_checkgradient2D_single(up, nextleft, left, centre, right, nextright, down);
		c_flag = node_checkgradient2D_comp_single(up.comp, nextleft.comp, left.comp, centre.comp, right.comp, nextright.comp, down.comp);
		if(!intf_flag && !c_flag) return centre.comp;

		//gradphi = (right.phi[0] - left.phi[0])/(2.0*dx);
		//gradc = (right.comp - left.comp)/(2.0*dx);
		//gradphi = (right.phi[0] - centre.phi[0])/(dx);
		//gradc = (right.comp - centre.comp)/(dx);
		gradphi = (-nextright.phi[0] + 8.0*right.phi[0] -  8.0*left.phi[0] + nextleft.phi[0])/(12.0*dx);
		gradc = (-nextright.comp + 8.0*right.comp -  8.0*left.comp + nextleft.comp)/(12.0*dx);
		
		
		double k = 1.0;
		laplacian_phi = (-nextright.phi[0] + 16.0*right.phi[0] - 30.0*centre.phi[0] + 16.0*left.phi[0] - nextleft.phi[0])/(12.0*dx*dx);
		laplacian_c = (-nextright.comp + 16.0*right.comp - 30.0*centre.comp + 16.0*left.comp - nextleft.comp)/(12.0*dx*dx);
		//laplacian_phi = ((right.phi[0] + left.phi[0] - 2.0*centre.phi[0])/(dx*dx) + (up.phi[0] + down.phi[0] - 2.0*centre.phi[0])/(dx*dx));
		//laplacian_c = ((right.comp + left.comp - 2.0*centre.comp)/(dx*dx) + (up.comp + down.comp - 2.0*centre.comp)/(dx*dx));
		double t1 = (gradc*gradphi*(1 - 2.0*centre.comp)*(1 - 2.0*centre.phi[0]));
		double t2 = centre.comp*(1 - centre.comp)*( (1.0 - 2.0*centre.phi[0])*laplacian_phi);
		double t3 = (-1)*centre.comp*(1 - centre.comp)*(2.0*gradphi*gradphi);
		c = centre.comp + nd_dt*(nd_diffusivity*laplacian_c  - k*nd_diffusivity*omega*alpha*(t1 + t2 + t3));
		if (((index == 94) || (index==2) || (index ==3)|| (index == 495) || (index == 496) || (index == 497) || (index == 498) || (index == 499) || (index == 500) || (index == 491) || (index == 492) || (index == 493) || 
			(index == 494) || (index == 485) || (index == 486) || (index == 487) || (index == 488) || (index == 489) || (index == 490) || (index == 4) || (index == 5) || (index == 6))
			 && (t%10000 == 1 ))
			printf("%d Phi : %lf, composition : %lf, updatecomp: %lf, gradphi : %lf, gradc : %lf, laplacephi : %lf, laplace_c : %e, lcomp:%e, rcomp:%e, term laplace : %e, other: %e, total : %e, t1 : %lf, t2 : %lf, t3 : %lf Flux = %e Potential = %e\n",
				index, centre.phi[0], centre.comp, c, gradphi, gradc, laplacian_phi, laplacian_c, left.comp, right.comp, nd_diffusivity*laplacian_c, nd_diffusivity*omega*alpha*(t1 + t2 + t3) ,
			    (c - centre.comp)/nd_dt,
				(t1), (t2), (t3), (nd_diffusivity*gradc - nd_diffusivity*omega*alpha*centre.comp*(1.0 - centre.comp)*(1.0 - 2.0*centre.phi[0])*gradphi), log(centre.comp/(1.0-centre.comp)) -omega*alpha*centre.phi[0]*(1.0 - centre.phi[0])  );
		
		return c;
	}
	return 0.0;
}

/******************************************************************/

double node_calcmu(node n)
{
	extern double omega, alpha;
	return(log(n.comp/(1 - n.comp)) - omega*alpha*n.phi[0]*(1 - n.phi[0]));
}

double node_csolver_singlecons(char *solver, int index, node up, node nextleft, node left, node centre, node right, node nextright, node down, double dx, int t)
{
	if(!strcmp(solver, "equalcomp"))
	{
		double c;
		double gradphi, gradc;
		double laplacian_phi, laplacian_c;
		double mu[3];
		extern double nd_dt;
		extern double nd_diffusivity;
		extern double omega, alpha;
		extern int phi_timestep;
		double tol = 1e-12;
		if(alpha < tol) return centre.comp;
		int intf_flag;
		int c_flag;

		/*###########*/
		/*
		/*mu[0]  mu[1]  mu[2]
		/*
		/*###########*/
		mu[0] = node_calcmu(left);
		mu[1] = node_calcmu(centre);
		mu[2] = node_calcmu(right);

		intf_flag = node_checkgradient2D_single(up, nextleft, left, centre, right, nextright, down);
		c_flag = node_checkgradient2D_comp_single(up.comp, nextleft.comp, left.comp, centre.comp, right.comp, nextright.comp, down.comp);
		if(!intf_flag && !c_flag) return centre.comp;

		double M[3];
		/** M[0]  M[1]   M[2] */
		M[0] = nd_diffusivity*left.comp*(1 - left.comp);
		M[1] = nd_diffusivity*centre.comp*(1 - centre.comp);
		M[2] = nd_diffusivity*right.comp*(1 - right.comp);



		//gradphi = (right.phi[0] - left.phi[0])/(2.0*dx);
		//gradc = (right.comp - left.comp)/(2.0*dx);
		//gradphi = (right.phi[0] - centre.phi[0])/(dx);
		//gradc = (right.comp - centre.comp)/(dx);

		c = centre.comp + nd_dt*( ((M[2]+M[1])/2.0)*(mu[2] - mu[1]) - ((M[0]+M[1])/2.0)*(mu[1] - mu[0]))/(dx*dx);
		
		
		//c = centre.comp + nd_dt*(nd_diffusivity*laplacian_c  - k*nd_diffusivity*omega*alpha*(t1 + t2 + t3));
		//if (((index == 94) || (index==2) || (index ==3)|| (index == 495) || (index == 496) || (index == 497) || (index == 498) || (index == 499) || (index == 500) || (index == 491) || (index == 492) || (index == 493) || 
		//	(index == 494) || (index == 485) || (index == 486) || (index == 487) || (index == 488) || (index == 489) || (index == 490) || (index == 4) || (index == 5) || (index == 6))
		//	 && (t%10000 == 1 ))
		//	printf("%d Phi : %lf, composition : %lf, updatecomp: %lf, gradphi : %lf, gradc : %lf, laplacephi : %lf, laplace_c : %e, lcomp:%e, rcomp:%e, term laplace : %e, other: %e, total : %e, t1 : %lf, t2 : %lf, t3 : %lf Flux = %e Potential = %e\n",
		//		index, centre.phi[0], centre.comp, c, gradphi, gradc, laplacian_phi, laplacian_c, left.comp, right.comp, nd_diffusivity*laplacian_c, nd_diffusivity*omega*alpha*(t1 + t2 + t3) ,
		//	    (c - centre.comp)/nd_dt,
		//		(t1), (t2), (t3), (nd_diffusivity*gradc - nd_diffusivity*omega*alpha*centre.comp*(1.0 - centre.comp)*(1.0 - 2.0*centre.phi[0])*gradphi), log(centre.comp/(1.0-centre.comp)) -omega*alpha*centre.phi[0]*(1.0 - centre.phi[0])  );
		
		return c;
	}
	return 0.0;
}