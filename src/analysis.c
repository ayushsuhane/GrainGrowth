#include<stdio.h>
#include<stdlib.h>
#include"analysisfunctions.h"
#include<string.h>
#include"nodefunctions.h"
#include"Global.h"

int Analysisinit(char *init)
{
	extern int ana_grain1, ana_grain2;
	extern double ana_phi;
	extern char **ana_keyword;
	ana_keyword = (char **)malloc(10 * sizeof(char *));
	for ( int i=0; i<10 ; i++) ana_keyword[i] = (char *)malloc(10*sizeof(char));
	if (init == NULL)
	{
		ana_grain1 = 0;
		ana_grain2 = 0;
		ana_phi = 0.0;
		ana_keyword[0] = NULL;
		printf("No analysis included\n");
		return 0;
	}
	// Initialize the extern parameters
	char tmp[15], tmp_g1[10], tmp_g2[10], tmp_phi[10];
	int grain1, grain2;
	double phival;

	sscanf(init, "%15s : %10s %10s %10s", tmp, tmp_g1, tmp_g2, tmp_phi);
	
	if(!strcmp(tmp, "track"))
	{
		ana_grain1 = atoi(tmp_g1);
		ana_grain2 = atoi(tmp_g2);
		ana_phi = atof(tmp_phi);
		strcpy(ana_keyword[0], tmp);

	} 
	printf("Initializing tracking of interface between %d and %d \n", ana_grain1, ana_grain2);
	printf("Keyword activated : %s\n", ana_keyword[0]);
	extern FILE *ana;
	ana = fopen("../output/analysis-linear.txt", "w");
	//ana = fopen("//mnt//d//Work//Code//C//GrainGrowth//output//analysis-linear.txt", "w");
	if(ana == NULL)
	{
		printf("Couldnot open file for analysis\n");
		perror("Error");
		exit(2);
	}
	printf("All the analysis is written in analysis.txt\n");

	return 1;		
}

void Analysis(double time)
{
	printf("Inside analysis\n");
	printf("Time = %e\n", time);
	extern char **ana_keyword;
	extern FILE *ana;
	double pos, c;
	if(!strcmp(ana_keyword[0], "track")) 
	{
		extern int ana_grain1, ana_grain2;
		extern double ana_phi;
		extern int gridx;
		//pos = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, (gridx/2), "pos" );
		//c = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, (gridx/2), "comp" );
		pos = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, gridx, "pos" );
		c = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, gridx, "comp" );
		fprintf(ana, "%e %e %e\n", time, pos, c);
	}
	double totalcomp = analysis_totalcomposition();
	printf("Total Composition = %e\n", totalcomp);
	return;
}


double analysis_track(int grain1, int grain2, double phi, int xstart, int xstop, char *key)
{
	/*
	An interface is characterized by two different grains. 
	Given the indices of two grains, this functions tracks the position of 
	the midpoint of the interface, which can be used to evaluate the velocity.
	*/
	printf("Inside track function\n");
	extern int gridsize;
	extern node *grid;
	extern int gridx, gridy;
	
	int y;
	//Choosing the index for 1D and 2D
	if(gridy == 1) y = 0;
	else y = (int)(gridy/2 - 1);

	double big = 1.0, small = 0.0; // Phi Value
	int bigpos =  0, smallpos = 0; // Position wrt index
	double val;
	int index;
	for(int i = xstart; i < xstop; i++)
	{
		index = y*gridx + i;
		if(grid[index].nactive > 1)
		{
			val = node_phival(grid[index], grain1);
			if ((val > phi) && (val < big) ) 
			{
				big = val;
				bigpos = i;
			}
			else if ((val - phi <= 0.0) && (val-small > 0.0))
			{
				small = val;
				smallpos = i;
			}
		}
		else continue;
	}
	//Linear interpolation between small and big to get the position
	printf("Big : %e, position : %d, Small : %e, position : %d\n", big, bigpos, small, smallpos);
	printf("Phi_1-big = %e, grain:%d, Phi_1-small = %e, grain:%d\n", grid[y*gridx + bigpos].phi[0], grid[y*gridx+bigpos].activegrain[0], grid[y*gridx +smallpos].phi[0], grid[y*gridx+smallpos].activegrain[0]);
	double pos,c;
	if(smallpos == bigpos && big == 1.0 && small == 0.0)
	{
		printf("ERROR: Interface is out of the domain. Either change the domain size, reduce the time or change the boundary conditions\n");
		exit(2);
	}
	pos = ((double)bigpos + ((big - small)/(bigpos - smallpos))*(phi - big));
	c = (double)(grid[y*gridx + bigpos].comp + ((grid[y*gridx + bigpos].comp - grid[y*gridx +smallpos].comp)/(bigpos - smallpos))*(pos - bigpos));
	printf("The interface is at %e and composition is %e\n", pos, c);
	if(!strcmp(key, "pos")) return pos;
	else if(!strcmp(key, "comp")) return c;
}

void Analysisend()
{	
	extern FILE *ana;
	fclose(ana);
	extern char **ana_keyword;
	for (int i = 0; i < 10 ; i++ ) free(ana_keyword[i]);
	free(ana_keyword);

	return;
}

double analysis_totalcomposition()
{
    extern node *grid;
    extern int gridsize;
    extern double initcomp;
    double totalcomp=0.0;
    for(int i=0; i<gridsize; i++)
    {
    	totalcomp += grid[i].comp;
    }
    printf("Change in percentage : %lf\n", 100.0*(initcomp - totalcomp/gridsize)/initcomp );
    return totalcomp/gridsize;

}





/**************************************************************/
/**************Single*****************************************/
int Analysisinit_single(char *init)
{
	extern double ana_phi;
	extern char **ana_keyword;
	ana_keyword = (char **)malloc(10 * sizeof(char *));
	for ( int i=0; i<10 ; i++) ana_keyword[i] = (char *)malloc(10*sizeof(char));
	if (init == NULL)
	{
		ana_phi = 0.0;
		ana_keyword[0] = NULL;
		printf("No analysis included\n");
		return 0;
	}
	// Initialize the extern parameters
	char tmp[15], tmp_phi[10];
	double phival;

	sscanf(init, "%15s : %10s", tmp, tmp_phi);
	
	printf("Keyword %s %s\n", tmp, tmp_phi);
	if(!strcmp(tmp, "track"))
	{
		ana_phi = atof(tmp_phi);
		strcpy(ana_keyword[0], tmp);
		printf("Initializing tracking of interface \n");
		printf("Keyword activated : %s\n", ana_keyword[0]);
		printf("Tracking interface at value %lf\n", ana_phi);
	} 
	extern FILE *ana;
	ana = fopen("../output/analysis-linear.txt", "w");
	//ana = fopen("//mnt//d//Work//Code//C//GrainGrowth//output//analysis-linear.txt", "w");
	if(ana == NULL)
	{
		printf("Couldnot open file for analysis\n");
		perror("Error");
		exit(2);
	}
	printf("All the analysis is written in analysis.txt\n");

	return 1;		
}

void Analysis_single(double time)
{
	printf("Inside analysis\n");
	printf("Time = %e\n", time);
	extern char **ana_keyword;
	extern FILE *ana;
	double pos, c, max_c;
	if(!strcmp(ana_keyword[0], "track")) 
	{
		extern double ana_phi;
		extern int gridx;
		//pos = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, (gridx/2), "pos" );
		//c = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, (gridx/2), "comp" );
		pos = analysis_track_single(ana_phi, 0, gridx, "pos" );
		c = analysis_track_single(ana_phi, 0, gridx, "comp" );
		max_c = analysis_interface_maxc_single();
		fprintf(ana, "%e %e %e %e\n", time, pos, c, max_c);
	}
	double totalcomp = analysis_totalcomposition();
	printf("Total Composition = %e\n", totalcomp);
	return;
}

double analysis_interface_maxc_single()
{
	extern int gridx;
	extern node *grid;
	double maxc = 0.0;
	for (int i=0; i< gridx; i++)
	{
		if(grid[i].phi[0] > 0.0 && grid[i].phi[0] < 1.0) 
		{
			if(grid[i].comp > maxc) maxc = grid[i].comp;
		}
	}
	return maxc;
}

double analysis_track_single(double phi, int xstart, int xstop, char *key)
{
	/*
	An interface is characterized by two different grains. 
	Given the indices of two grains, this functions tracks the position of 
	the midpoint of the interface, which can be used to evaluate the velocity.
	*/
	printf("Inside track function\n");
	extern int gridsize;
	extern node *grid;
	extern int gridx, gridy;
	
	int y;
	//Choosing the index for 1D and 2D
	if(gridy == 1) y = 0;
	else y = (int)(gridy/2 - 1);

	double big = 1.0, small = 0.0; // Phi Value
	int bigpos =  0, smallpos = 0; // Position wrt index
	double val;
	int index;
	for(int i = xstart; i < xstop; i++)
	{
		index = y*gridx + i;
		val = grid[index].phi[0];
		//printf("index : %d, Phi : %lf, compare_phi: %lf\n", index, val, phi);
		if ((val > phi) && (val <= big) ) 
		{
			big = val;
			bigpos = i;
		}
		else if ((val - phi <= 0.0) && (val-small > 0.0))
		{
			small = val;
			smallpos = i;
		}
	}
	
	//Linear interpolation between small and big to get the position
	printf("Big : %e, position : %d, Small : %e, position : %d\n", big, bigpos, small, smallpos);
	printf("Phi_1-big = %e, grain:%d, Phi_1-small = %e, grain:%d\n", grid[y*gridx + bigpos].phi[0], grid[y*gridx+bigpos].activegrain[0], grid[y*gridx +smallpos].phi[0], grid[y*gridx+smallpos].activegrain[0]);
	double pos,c;
	if(smallpos == bigpos && big == 1.0 && small == 0.0)
	{
		printf("ERROR: Interface is out of the domain. Either change the domain size, reduce the time or change the boundary conditions\n");
		exit(2);
	}
	pos = ((double)bigpos + ((big - small)/(bigpos - smallpos))*(phi - big));
	c = (double)(grid[y*gridx + bigpos].comp + ((grid[y*gridx + bigpos].comp - grid[y*gridx +smallpos].comp)/(bigpos - smallpos))*(pos - bigpos));
	printf("The interface is at %e and composition is %e\n", pos, c);
	if(!strcmp(key, "pos")) return pos;
	else if(!strcmp(key, "comp")) return c;
	else return 0.0;
}