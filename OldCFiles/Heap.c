#include "structure.h"
#include "Heap.h"
#include <math.h>

// HEAP  Implementation of a min heap with tuple in it.
// We defined that number 1 is the id of the vertex and number two is the distance.


void hpInit(Heap* h)
{

	h->array = NULL;
	h->size = 0;
	h->capacity = 0;
}

void swap(int* a, int* b)
{
	int temp = *a; 
	*a = *b;
	*b = temp;
}

Heap* hpCreate()
{
	Heap* h = malloc(sizeof(Heap));
	hpInit(h);
	return h;
}

void hpAddAlloc(Heap* h)
{
	h->array = realloc(h->array, (h->size + REALLOCSIZE)*sizeof(Tuple));
	for(int i = 0; i<REALLOCSIZE; i++)
	{	
		h->array[h->size + i] = tCreate();
	}
	h->capacity += REALLOCSIZE;
}


void hpAddElement(Heap* h, Tuple* tuple)
{
	if (h->size + 1 > h->capacity){
		hpAddAlloc(h);
	}
	int i = h->size;
	h->size++;
	h->array[i] = tuple;
	while(i !=0 && harray(h,((i-1)/2))->nb2 > harray(h,i)->nb2)
	{	
		swap(&harray(h,((i-1)/2))->nb1,&harray(h,i)->nb1); // We swap the id of the vertex
		swap(&harray(h,((i-1)/2))->nb2,&harray(h,i)->nb2); // And the distance.
		i = (i-1)/2;
	}
}

void hpHeapify(Heap* h, int i)
{
	int min = i; // Actual node
	int left = 2*i +2; // Child left
	int right = 2*i+1; // Child right

	if (left < h->size && harray(h,left)->nb2 < harray(h,min)->nb2)
	{
		min = left; // If value at left is smaller, we change our min index
	}
	if (right < h->size && harray(h,right)->nb2 < harray(h,min)->nb2)
	{
		min = right;  // If value at right is smaller, we change our min index
	}
	if (min != i) // If our i is different as min, that mean we need to change.
	{
		swap(&harray(h,i)->nb1,&harray(h,min)->nb1); // We swap the id of the vertex
		swap(&harray(h,i)->nb2,&harray(h,min)->nb2); // And the distance. // We swap value
		hpHeapify(h,min); // And we call again to continue to change
	}

}

void hpIncreaseVal(Heap* h, int id, int newValue)
{

	int index = 0;

	for (;index<h->size; index++)
	{
		if(harray(h,index)->nb1 == id) // We check if the id we search is the id of the tuple.
		{
			break;
		}
	}

    if (index >= h->size || harray(h,index)->nb2 <= newValue)
    {
        printf("Invalid index or new value is not smaller\n");
        return;
    }

    harray(h,index)->nb2 = newValue;
    while (index > 0 && harray(h,((index - 1) / 2))->nb2 > harray(h,index)->nb2)
    {
        swap(&harray(h,index)->nb1, &harray(h,((index - 1) / 2))->nb1);
        swap(&harray(h,index)->nb2, &harray(h,((index - 1) / 2))->nb2);
        index = (index - 1) / 2;
    }
}

Tuple* hpExtractMin(Heap* h)
{
	Tuple* root = h->array[0];
	h->array[0] = h->array[h->size-1];
	h->size--;

	hpHeapify(h,0);
	return root;
}

void hpDestroyHeap(Heap* h)
{
	free(h->array);
	free(h);
}

void printHeap(Heap *h)
{
    for (int i = 0; i < h->size; ++i)
        printf("Id : %d, dist : %d \n", harray(h,i)->nb1,harray(h,i)->nb2);
    printf("\n");
}