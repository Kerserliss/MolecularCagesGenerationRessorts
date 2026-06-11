#ifndef __HEAP_H
#define __HEAP_H

typedef struct 
{
    Tuple** array;
	int size; // Number of element in the heap
	int capacity; // Allocated size
}Heap;



void hpInit(Heap* h);
void swap(int* a, int* b);
Heap* hpCreate();
void hpAddAlloc(Heap* h);
void hpAddElement(Heap* h, Tuple* elt);
void hpHeapify(Heap* h, int i);
Tuple* hpExtractMin(Heap* h);
void hpDestroyHeap(Heap* h);
void printHeap(Heap *heap);
void hpIncreaseVal(Heap* h, int id, int newValue);
  
#endif