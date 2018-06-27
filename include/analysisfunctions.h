#ifndef ANALYSISFUNCTIONS_H
#define ANALYSISFUNCTIONS_H


int ana_grain1, ana_grain2;
float ana_phi;
char **ana_keyword;
FILE *ana;

int Analysisinit(char *);
void Analysis(float );
float analysis_track(int , int , float, int, int);
void Analysisend();
#endif