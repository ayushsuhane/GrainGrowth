#include<stdio.h>
#include<stdlib.h>
#include"Functions.h"
#include"Global.h"
#include"nodefunctions.h"
#include<string.h>

void Output(int t)
{
	
	char FILENAME1[] = "../output/finaltime1.txt";
	node_outputphi_tofile(FILENAME1, 1);
	char FILENAME2[] = "../output/finaltime2.txt";
	node_outputphi_tofile(FILENAME2, 2);
	char COMP_FILENAME[] = "../output/composition.txt";
	node_outputcomp_tofile(COMP_FILENAME);
}
void output_visual(int t, int grainnum)
{
	char FILENAME[100];
	sprintf( FILENAME, "../output/visual%d_%d.txt", grainnum, t );  
	node_visual(FILENAME, grainnum);
}