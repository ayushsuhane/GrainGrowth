#include<stdio.h>
#include<stdlib.h>
#include"analysisfunctions.h"
#include<string.h>
#include"nodefunctions.h"
#include"Global.h"

int Analysisinit(char *init)
{
	extern int ana_grain1, ana_grain2;
	extern float ana_phi;
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
	float phival;

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
	printf("All the analysis is written in analysis.txt\n");

	return 1;		
}

void Analysis(float time)
{
	extern char **ana_keyword;
	extern FILE *ana;
	float pos;
	if(!strcmp(ana_keyword[0], "track")) 
	{
		extern int ana_grain1, ana_grain2;
		extern float ana_phi;
		extern int gridx;
		pos = analysis_track(ana_grain1, ana_grain2, ana_phi, 0, gridx);
		fprintf(ana, "%f %f\n", time, pos);
	}
	return;
}


float analysis_track(int grain1, int grain2, float phi, int xstart, int xstop)
{
	/*
	An interface is characterized by two different grains. 
	Given the indices of two grains, this functions tracks the position of 
	the midpoint of the interface, which can be used to evaluate the velocity.
	*/
	extern int gridsize;
	extern node *grid;
	extern int gridx, gridy;
	
	int y;
	//Choosing the index for 1D and 2D
	if(gridy == 1) y = 0;
	else y = (int)(gridy/2 - 1);

	float big = 1, small = 0; // Phi Value
	int bigpos, smallpos; // Position wrt index
	float val;
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
			if ((val < phi) && ( val > small))
			{
				small = val;
				smallpos = i;
			}
		}
		else continue;
	}
	//Linear interpolation between small and big to get the position
	//printf("Big : %f, position : %d, Small : %f, position : %d\n", big, bigpos, small, smallpos);
	float pos;
	pos = ((float)bigpos + ((big - small)/(bigpos - smallpos))*(phi - big));
	printf("The interface is at %f\n", pos);
	return pos;
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