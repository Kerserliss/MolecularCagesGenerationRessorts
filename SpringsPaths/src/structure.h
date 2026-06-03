#ifndef __STRUCTURE_H
#define __STRUCTURE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REALLOCSIZE 4

#define coords(v) (v)->coords
#define size(l) (l)->size
#define sizealloc(l) (l)->sizealloc
#define elts(l, i) (l)->elts[(i)]

#define gsize(g) (g)->size
#define v_id(v) (v)->id

// LIST

typedef struct
{
	int* elts; // Array of elts
	unsigned size; // Lenght of the list
} List_t;

// Liste
/**
 * @brief Initializes a new list.
 *
 * This function initializes a new list by setting its elements to NULL and size to 0.
 *
 * @param l Pointer to the list to initialize.
 */
void lstInit(List_t *);

/**
 * @brief Adds allocated memory to the list.
 *
 * This function adds allocated memory to the list to accommodate new elements.
 *
 * @param l Pointer to the list to which memory is added.
 */
void lstAddAlloc(List_t *);

/**
 * @brief Returns the number of elements in the list.
 *
 * This function returns the number of elements in the list.
 *
 * @param l Pointer to the list.
 * @return The number of elements in the list.
 */
unsigned lstNbElements(List_t *);

/**
 * @brief Returns the index of the first available free element in the list.
 *
 * This function returns the index of the first available free element in the list.
 * If no free element is available, it adds and allocates memory for more elements.
 *
 * @param l Pointer to the list.
 * @return The index of the first available free element in the list.
 */
unsigned lstGetIndiceFree(List_t *);

/**
 * @brief Returns the index of the element with the specified identifier in the list.
 *
 * This function returns the index of the element with the specified identifier in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to search for.
 * @return The index of the element with the specified identifier, or -1 if not found.
 */
unsigned lstGetIndice(List_t *, unsigned);

/**
 * @brief Checks if an element with the specified identifier exists in the list.
 *
 * This function checks if an element with the specified identifier exists in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to check for.
 * @return 1 if the element with the specified identifier exists, 0 otherwise.
 */
unsigned lstCheck(List_t *, unsigned);

/**
 * @brief Adds an element to the list if it does not exist.
 *
 * This function adds an element to the list if it does not already exist.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be added.
 */
void lstAddElement(List_t *, unsigned);

/**
 * @brief Removes an element with the specified identifier from the list.
 *
 * This function removes an element with the specified identifier from the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be removed.
 */
void lstRemoveElement(List_t *, unsigned);

/**
 * @brief Creates a new integer list.
 *
 * This function allocates memory for a new integer list and initializes its fields.
 *
 * @return Pointer to the newly allocated integer list.
 */
List_t *lstCreate();

/**
 * @brief Creates a copy of an integer list.
 *
 * This function creates a copy of an integer list by copying its elements.
 *
 * @param l Pointer to the integer list to be copied.
 * @return Pointer to the newly created copy of the integer list.
 */
List_t *lstCopy(List_t *);

/**
 * @brief Copies the list with an offset of the numbering of elements
 * according to the value of the array passed as an argument.
 *
 * @param l List to copy.
 * @param shifts Array of offsets.
 * @return (List_t*) Copied list with shifts.
 */
List_t *lstCopyWithShift(List_t *l, int *mod_pos_nei);

/**
 * @brief Merges two integer lists.
 *
 * This function merges two integer lists, adding elements from the second list to the first.
 * The function then deletes both input lists and returns the merged list.
 *
 * @param l1 First integer list.
 * @param l2 Second integer list.
 * @return Merged integer list.
 */
List_t *lstAddList(List_t *, List_t *);

/**
 * @brief Deletes an integer list and frees the associated memory.
 *
 * This function deletes an integer list and frees the memory associated with its elements.
 *
 * @param l Pointer to the integer list to be deleted.
 */
void lstDelete(List_t *);

// VERTEX


//Point_t
typedef struct
{
	double x;
	double y;
	double z;
} Point_t;


Point_t *ptCreate(double x, double y, double z);
void ptDestroy(Point_t *p);

//Vertex
typedef struct
{
	int id; // Name of the vertex
	List_t* neighborhood; // List of neighbors of the vertex
	Point_t* coords; // Coords of the vertex
} Vertex;

Vertex* vtxInit(Vertex* v, int id, double x, double y, double z);
Vertex* vtxCreate(int id, double x, double y, double z);
void vtxAddNeighboord(Vertex* v, int id_v);
void vtxRemoveNeighboord(Vertex* v, int id_v);
void vtxGetSizeNeighboord(Vertex* v);
void vtxDestroy(Vertex* v);
//Graph

typedef struct 
{
	Vertex** vertices; // Array Of vertices
	int size;
	int capacity;
}Graph_t; // Est ce que ça sert ?

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
float hpExtractMin(Heap* h);
void hpDestroyHeap(Heap* h);
void printHeap(Heap *heap);

//Graph
Graph_t* gInit(Graph_t* g);
Graph_t* gCreate();
void gAddVertex(Graph_t* g,Vertex* v);
void gRemoveVertex(Graph_t* g, Vertex* v, int id);
Vertex* gGetVertex(Graph_t*g,int id);
void gDestroy(Graph_t* g);
#endif