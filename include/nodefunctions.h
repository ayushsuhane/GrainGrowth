
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
	int activegrain[2];
	double phi[2];
	char phase[2];
	double comp;
} node;



node *grid; //A one dimensional pointer to hold the values of phase field parameters at each node

node node_initialize();
node * node_alloc(int);
void node_dealloc(node);

int node_checkgradient(node , node, node);
int node_checkgradient2D(node, node, node, node, node);
int node_checkgradient2D_comp(double, double, double, double, double);

int node_checkgradient2D_single(node, node, node, node, node, node, node);
int node_checkgradient2D_comp_single(double, double, double, double, double, double, double);

node  node_buildlist(node, node, node);
node node_buildlist2D(node, node, node, node, node);
node _node_checklist(node, node);
node node_membersort(node );

double node_gradient2D(int, node, node, double);
double node_laplacian1D(int, node, node, node, double);
double node_laplacian2D(int, node, node, node, node, node, double);
double node_laplacian2D_comp(double, double, double, double, double, double);
double node_phival(node, int);
double node_termsumphi_comp(node, node, node, node, node, node, double, int, int);
double node_prodphi_comp(node );
void node_handleBC(int, int, node *, char * );

void node_print(node, int); 
node node_transfer(node, node);
node node_simplex(node);
void node_outputphi_tofile(char *, int);
void node_outputcomp_tofile(char *);

void node_visual(char *, int);

//Not actually node functions but phase field functions
double node_interface(double, double);
double node_addforce(int, double, double, int);
node node_update(node ,int);
node node_update_single(node);
node node_phisolver(char *, int, node, node, node, node, node, node,  double, int);
node node_phisolver_single(char *, int, node, node, node, node, node, node, node, node,  double, int);

double node_csolver(char *, int, node, node, node, node, node, node, double, int);
double node_csolver_single(char *, int, node, node, node, node, node, node, node, double, int);

double node_calcflux_x(node, node, node, node, node);
#endif