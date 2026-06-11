#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structure.h"

List_e* lsteCreate() {
	List_e* list = malloc(sizeof(List_e));
	list->edges = NULL;
  	list->size_elts = 0;
  	list->allocated_size = 0;
  	return list;
}

void lsteAddAlloc(List_e* list)
{
	int i;

	list->elts =  realloc(list->elts, (l->size_elts + REALLOCSIZE) * sizeof(Edge));

	for(i = 0;i<REALLOCSIZE; i++)
	{
		list->edges[list->size_elts] += edCreate(-1,-1,-1);
	}
	list->allocated_size += REALLOCSIZE;
}


void lsteAddElement(List_e list, Edge e)
{
	if(list->size_elts + 1 > list->allocated_size)
		lsteAddAlloc(list);
	list->edges[l->size_elts] = e; 
	l->size_elts +=1;
}

void lstRemoveElement(List_e list, Edge e)
{	
	int i = 0;
	while(i<list->size_elts)
	{
		if (list->edges[i].id = e.id)
		{
			while(i<list->size_elts)
			{
				list->edges[i] = list->edges[i+1];
				i++;
			}
			list->size_elts-=1;
		}
	i++;
	}
}

Edge lsteFindEdgeID(List_e* list,int id)
{
	for(int i = 0; i<list->size_elts;i++)
	{
		if(list->edges[i].id = id)
			return list->edges[i];
	}
	if(i == list->size_elts)
		return NULL;
}

void lstePrint(List_e* list)
{
	for(int i = 0; i<size_elts; i++)
}