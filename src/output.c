#include<stdio.h>
#include<stdlib.h>
#include"Functions.h"
#include"Global.h"
#include"nodefunctions.h"

void Output(int t)
{
	
	char FILENAME[] = "../output/finaltime.txt";
	node_outputphi_tofile(FILENAME, 1);

}