#ifndef ANALYSISFUNCTIONS_H
#define ANALYSISFUNCTIONS_H


int ana_grain1, ana_grain2;
double ana_phi;
char **ana_keyword;
FILE *ana;

int Analysisinit(char *);
void Analysis(double );
double analysis_track(int , int , double, int, int, char *);
double analysis_totalcomposition();
void Analysisend();

int Analysisinit_single(char *);
void Analysis_single(double );
double analysis_track_single(double, int, int, char *);
double analysis_interface_maxc_single();
#endif