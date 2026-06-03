#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure.c"
#include "Heap.c"
#include "constant.h"

int InitMat(int line, int col)
{
	int m[line][col];
	for (int i = 0; i < line; i++)
	{
		 for (int j = 0; j < col; j++)
		 {
		 	m[i][j] = 0;
		 }
	}
	return m
}

int Dijkstra(Graph_t *g, int source, int target, int* correspondance)
{
	// Matrice initialization. Heap initialization.
	mdist[g->size];
	mprev[g->size];
	Heap* h = hpCreate();

	for(int i = 0; i<g->size; i++)
	{
		if(source == correspondance[i])
		{
			mdist[i]=0;
			mprev[i]= -1;
			hpAddElement(h,i);
		}
		else
		{ 
			mdist[i][j] = g->size;
			mprev[i][j] = -1;
			hpAddElement(h,i);
		}
	}

	// Start of the algorithm
	while (h->size !=0)
	{
		int id = hpExtractMin(h);
		Vertex* v = gGetVertex()
	}
}	


float total_system_energy(graph* g)
{
	float e = 0;
	return 1.5	
}