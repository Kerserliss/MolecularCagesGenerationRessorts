#ifndef __KAMADA_KAWAI_H
#define __KAMADA_KAWAI_H

#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constant.h"
#include <math.h>
#include <lapacke.h>

typedef struct
{
	int* mat;
	int* correspondance;
}ResultDijkstra;

double euclidien_distance(Point_t* p1, Point_t* p2);

int** InitMatInt(int line, int col);
double** InitMatDouble(int line, int col);
ResultDijkstra* Dijkstra(Graph_t *g, int source);
int GetTabIndice(int* tab, int id, int n);

void rdInit(ResultDijkstra* rd);
ResultDijkstra* rdCreate();
void rdDestroy(ResultDijkstra* rd);

// Computing fonction
double* Hessienne(Graph_t*g,Vertex* v,double** matK, double** matL);
double Gradient_Norme(double* matE);
double* MatrixMultiplication(double* Mat1, int line1, int col1, double* Mat2, int line2, int col2);
double* Newton_Raphson(Graph_t* g, Vertex* v, double** matK, double** matL);
double* Energie_Gradient(Graph_t* g, Vertex* v, double** matK, double** matL);


#endif