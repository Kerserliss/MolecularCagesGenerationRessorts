#ifndef __HEAP_H
#define __HEAP_H

typedef struct 
{
    int *array;
	int size; // Number of element in the heap
	int capacity; // Allocated size
}Heap;

Heap* hpInit(Heap* h);
void swap(int* a, int* b);
Heap* hpCreate();
void hpAddAlloc(Heap* h);
void hpAddElement(Heap* h, int elt);
void hpHeapify(Heap* h, int i);
int hpExtractMin(Heap* h);
void hpDestroyHeap(Heap* h);
void printHeap(Heap *heap);
void hpIncreaseVal(Heap* h, int element, int newValue);

#endif