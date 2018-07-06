#ifndef FUNCTIONS_H
#define FUNCTIONS_H


int Initialize( char *, char *);
void Readinputs();
void NonDimensionalize();
void Parametrize(char *);
void Printsettings();
void Allocate();
void Setup();

void Check_initialstate();
void FreeMemory();

void Solver(char *, char *, int);
void Output(int );
void output_visual(int, int);



void randomTest();
#endif