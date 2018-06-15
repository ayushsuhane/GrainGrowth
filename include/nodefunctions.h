
#ifndef NODEFUNCTIONS_H
#define NODEFUNCTIONS_H
//Node Type
/***Structure**/
/*
Every grid point is a struct type which contains information about the active phases 
at that node. 
*/
typedef struct 
{
	int nactive;
	int *activegrain;
	float *phi;
	char *phase;
} node;



node *grid; //A one dimensional pointer to hold the values of phase field parameters at each node

node node_initialize();
node * node_alloc(int);
void node_dealloc(node);

int node_checkgradient(node , node, node);
node  node_buildlist(node, node, node);
node node_membersort(node );
float node_laplacian1D(int, node, node, node, float);
float node_phival(node, int);
void node_print(node, int); 
node node_transfer(node, node);
node node_simplex(node);
void node_outputphi_tofile(char *, int);

//Not actually node functions but phase field functions
float node_interface(float, float);
float node_addforce(int, float, float);
node node_update(node ,int);
#endif