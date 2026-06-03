#include "structure.h"
#include "Heap.h"
#include <math.h>

// HEAP  Implementation of a min heap

Heap* hpInit(Heap* h)
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
	h->array = realloc(h->array, (h->size + REALLOCSIZE)*sizeof(float));
	for(int i = 0; i<REALLOCSIZE; i++)
	{	
		h->array[h->size + i] = -1;
	}
	h->capacity += REALLOCSIZE;
}


void hpAddElement(Heap* h,int elt)
{
	if (h->size + 1 > h->capacity){
		hpAddAlloc(h);
	}
	int i = h->size;
	h->size++;
	h->array[i] = elt;
	while(i !=0 && h->array[(i-1)/2] > h->array[i])
	{	
		swap(&h->array[(i-1)/2],&h->array[i]);
		i = (i-1)/2;
	}
}

void hpHeapify(Heap* h, int i)
{
	int min = i; // Actual node
	int left = 2*i +1; // Child left
	int right = 2*i; // Child right

	if (left < h->size && h->array[left] < h->array[min])
	{
		min = left; // If value at left is smaller, we change our min index
	}
	if (right < h->size && h->array[right] < h->array[min])
	{
		min = right;  // If value at right is smaller, we change our min index
	}
	if (min != i) // If our i is different as min, that mean we need to change.
	{
		swap(&h->array[i], &h->array[min]); // We swap value
		hpHeapify(h,min); // And we call again to continue to change
	}

}

void hpIncreaseVal(Heap* h, int element, int newValue)
{

	int index = 0;

	for (index; index<h->size; index++)
	{
		if(h->array[index] == element)
		{
			break;
		}
	}

    if (index >= h->size || h->array[index] <= newValue)
    {
        printf("Invalid index or new value is not smaller\n");
        return;
    }

    h->array[index] = newValue;
    while (index != 0 && h->array[(index - 1) / 2] < h->array[index])
    {
        swap(&h->array[index], &h->array[(index - 1) / 2]);
        index = (index - 1) / 2;
    }
}

int hpExtractMin(Heap* h)
{
	int root = h->array[0];
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

void printHeap(Heap *heap)
{
    for (int i = 0; i < heap->size; ++i)
        printf("%d \n", heap->array[i]);
    printf("\n");
}