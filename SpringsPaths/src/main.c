#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure.c"
#include "Heap.c"

void main()
{
	printf("Hello World \n");
	Heap* h = hpCreate();
	for (int i =0; i< 10; i ++)
	{	
		printf("%d main \n", i);
		hpAddElement(h,(float)i);
	}
	printHeap(h);
	for(int i =30; i > 20; i --)
	{
		hpAddElement(h,(float)i);
	}
	printHeap(h);
	float r;
	r = hpExtractMin(h);
	printf("Root 1 %f \n",r);
	printHeap(h);
	r = hpExtractMin(h);
	printf("Root 2 %f \n",r);
	printHeap(h);
	hpDestroyHeap(h);
}