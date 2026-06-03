#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure.h"
#include "Heap.h"
#include "constant.h"
#include "Kamada_Kawai.h"

int** InitMat(int line, int col)
{
	int** m = malloc(line*sizeof(int*));
	for (int i = 0; i < line; i++)
	{
		m[i] = malloc(col*sizeof(int));
		 for (int j = 0; j < col; j++)
		 {
		 	m[i][j] = 0;
		 }
	}
	return m;
}

List_t* Dijkstra(Graph_t *g, int source, int target, int* correspondance)
{
	// Matrice initialization. Heap initialization.
	int mdist[g->size];
	int mprev[g->size];
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
			mdist[i]= g->size;
			mprev[i] = -1;
			hpAddElement(h,i);
		}
	}

	// Start of the algorithm
	while (h->size !=0)
	{
		int id = hpExtractMin(h);
		Vertex* v = gGetVertex(g,correspondance[id]);
		for (int i = 0; i<size(v->neighborhood); i++)
		{
			int temp = mdist[id] + 1; // +1 car distance tjr de 1
			if (temp < mdist[i])
			{
				hpIncreaseVal(h,i,temp);
				mdist[i] = temp;
				mprev[i] = id;
			}
		}
	}

	// Find path
	List_t* seq = lstCreate();
	int size = 0;
	int u = target;
	while(mprev[u] != -1 || u != source)
	{
		lstAddElement(seq,u);
		u = mprev[u];
	}
	return seq;

}	


float total_system_energy(Graph_t* g)
{
	float e = 0;
	return 1.5;
}