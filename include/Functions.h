#ifndef FUNCTIONS_H
#define FUNCTIONS_H

int Readcontrols(char *);
int Initialize( char *, char *);
void Readinputs();
void NonDimensionalize();
void Parametrize(char *);
void Printsettings();
void Allocate();
void Setup(char *);
void Setup_single(char *);


void Check_initialstate();
void FreeMemory();

void Solver(char *, char *, int);
void Solver_single(char *, char *, int);
void Output(int );
void output_visual(int, int);
void output_initialfile();



void randomTest();
#endif