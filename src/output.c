#include<stdio.h>
#include<stdlib.h>
#include"Functions.h"
#include"Global.h"
#include"nodefunctions.h"
#include<string.h>

void Output(int t)
{
	printf("Inside Output function\n");
	
	char FILENAME1[] = "../output/finaltime1.txt";
	node_outputphi_tofile(FILENAME1, 1);
	char FILENAME2[] = "../output/finaltime2.txt";
	node_outputphi_tofile(FILENAME2, 2);
	char COMP_FILENAME[] = "../output/composition.txt";
	node_outputcomp_tofile(COMP_FILENAME);
	output_initialfile();
}
void output_visual(int t, int grainnum)
{
	char FILENAME[100];
	extern int noutput;
	extern double nd_dt, nd_totaltime;
	sprintf( FILENAME, "../output/visual%d_%0.2e.txt", grainnum, t*nd_dt*real_totaltime/nd_totaltime);  
	node_visual(FILENAME, grainnum);

	char comp_FILENAME[100];
	sprintf( comp_FILENAME, "../output/visualcomp_%0.2e.txt", t*nd_dt*real_totaltime/nd_totaltime);  
	node_outputcomp_tofile(comp_FILENAME);
}

void output_initialfile()
{
	extern int flag_writetofile;
	if(flag_writetofile)
	{
		extern char initialfile[100];
		FILE *restart = fopen(initialfile, "w");
		if(restart == NULL)
		{
			printf("Couldnot open the file to write\n");
			exit(2);
		}
		extern int gridx, gridy;
		extern node *grid;
		int index;
		for (int i=0; i<gridx ; i++)
		{
			for (int j=0; j<gridy; j++)
			{
				index = j*(gridx) + i;
				fprintf(restart, "%d\n", index);
				fprintf(restart, "%d\n", grid[index].nactive);
				for(int k=0; k< grid[index].nactive; k++)
					fprintf(restart, "%d ", grid[index].activegrain[k]);
				fprintf(restart, "\n");
				for(int k=0; k< grid[index].nactive; k++)
					fprintf(restart, "%e ", grid[index].phi[k]);
				fprintf(restart, "\n");
				for(int k=0; k< grid[index].nactive; k++)
					fprintf(restart, "%c ", grid[index].phase[k]);
				fprintf(restart, "\n");
				fprintf(restart, "%e\n", grid[index].comp);
				fprintf(restart, "\n");
			}
			fprintf(restart,"\n");
		}
		fclose(restart);
	}
	else return;
}