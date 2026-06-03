#include "structure.h"
#include <math.h>

/**************************************/
/* LIST *******************************/
/**************************************/
/**
 * @brief Initializes a new list.
 *
 * This function initializes a new list by setting its elements to NULL and size to 0.
 *
 * @param l Pointer to the list to initialize.
 */
void lstInit(List_t *l) {
  l->elts = NULL;
  l->size = 0;
}

/**
 * @brief Adds allocated memory to the list.
 *
 * This function adds allocated memory to the list to accommodate new elements.
 *
 * @param l Pointer to the list to which memory is added.
 */
void lstAddAlloc(List_t *l) {
  int i;

  l->elts = realloc(l->elts, (size(l) + REALLOCSIZE) * sizeof(int));

  for (i = 0; i < REALLOCSIZE; i++)
    l->elts[l->size + i] = -1;

  size(l) += REALLOCSIZE;
}

/**
 * @brief Returns the number of elements in the list.
 *
 * This function returns the number of elements in the list.
 *
 * @param l Pointer to the list.
 * @return The number of elements in the list.
 */
unsigned lstNbElements(List_t *l) {
  int cpt;

  for (cpt = size(l); cpt > 0 && elts(l, cpt - 1) == -1; cpt--)
    ;

  return cpt;
}

/**
 * @brief Returns the index of the first available free element in the list.
 *
 * This function returns the index of the first available free element in the list.
 * If no free element is available, it adds and allocates memory for more elements.
 *
 * @param l Pointer to the list.
 * @return The index of the first available free element in the list.
 */
unsigned lstGetIndiceFree(List_t *l) {
  int i;

  for (i = 0; i < size(l); i++)
    if (elts(l, i) == -1)
      return i;

  lstAddAlloc(l);
  return i;
}

/**
 * @brief Returns the index of the element with the specified identifier in the list.
 *
 * This function returns the index of the element with the specified identifier in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to search for.
 * @return The index of the element with the specified identifier, or -1 if not found.
 */
unsigned lstGetIndice(List_t *l, unsigned id) {
  int i;

  for (i = 0; i < size(l) && elts(l, i) != -1; i++)
    if (elts(l, i) == id)
      return i;

  return -1;
}

/**
 * @brief Checks if an element with the specified identifier exists in the list.
 *
 * This function checks if an element with the specified identifier exists in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to check for.
 * @return 1 if the element with the specified identifier exists, 0 otherwise.
 */
unsigned lstCheck(List_t *l, unsigned id) {

  if (lstGetIndice(l, id) == -1)
    return 0;
  return 1;
}

/**
 * @brief Adds an element to the list if it does not exist.
 *
 * This function adds an element to the list if it does not already exist.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be added.
 */
void lstAddElement(List_t *l, unsigned id) {

  int i;
  if (lstGetIndice(l, id) == -1) {
    i = lstGetIndiceFree(l);
    l->elts[i] = id;
  }
}

/**
 * @brief Removes an element with the specified identifier from the list.
 *
 * This function removes an element with the specified identifier from the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be removed.
 */
void lstRemoveElement(List_t *l, unsigned id) {

  int i = lstGetIndice(l, id);

  if (i != -1) {
    while (i < size(l) - 1 && elts(l, i) != -1) {
      elts(l, i) = elts(l, i + 1);
      i++;
    }
    elts(l, i) = -1;
  }
}

/**
 * @brief Creates a new integer list.
 *
 * This function allocates memory for a new integer list and initializes its fields.
 *
 * @return Pointer to the newly allocated integer list.
 */
List_t *lstCreate() {

  List_t *l = malloc(sizeof(List_t));
  lstInit(l);

  return l;
}

/**
 * @brief Creates a copy of an integer list.
 *
 * This function creates a copy of an integer list by copying its elements.
 *
 * @param l Pointer to the integer list to be copied.
 * @return Pointer to the newly created copy of the integer list.
 */
List_t *lstCopy(List_t *l) {

  int i;
  List_t *copy = lstCreate();

  size(copy) = lstNbElements(l);

  copy->elts = malloc(size(copy) * sizeof(int));

  for (i = 0; i < size(copy); i++)
    elts(copy, i) = elts(l, i);

  return copy;
}

/**
 * @brief Copies the list with an offset of the numbering of elements
 * according to the value of the array passed as an argument.
 *
 * @param l List to copy.
 * @param shifts Array of offsets.
 * @return (List_t*) Copied list with shifts.
 */
List_t *lstCopyWithShift(List_t *l, int *shifts) {

  int i;
  List_t *copy = lstCreate();

  size(copy) = lstNbElements(l);

  copy->elts = malloc(size(copy) * sizeof(int));

  for (i = 0; i < size(copy); i++)
    elts(copy, i) = elts(l, i) - shifts[elts(l, i)];

  return copy;
}

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
List_t *lstAddList(List_t *l1, List_t *l2) {

  int i;

  List_t *out = lstCopy(l1);

  for (i = 0; i < size(l2) && elts(l2, i) != -1; i++)
    lstAddElement(out, elts(l2, i));

  lstDelete(l2);
  lstDelete(l1);
  return out;
}

/**
 * @brief Deletes an integer list and frees the associated memory.
 *
 * This function deletes an integer list and frees the memory associated with its elemencurrentElemts.
 *
 * @param l Pointer to the integer list to be deleted.
 */
void lstDelete(List_t *l) {

  free(l->elts);
  free(l);
}

void lstPrint(List_t *l)
{
	for(int i = 0; i<size(l); i++)
	{
		printf("%d \n",elts(l,i));

	}
	printf("\n");
}

/**************************************/
/* Point *******************************/
/**************************************/

Point_t *ptCreate(double x, double y, double z){
	Point_t *p = malloc(sizeof(Point_t));
	p->x = x;
	p->y = y;
	p->z = z;
	return p;
}

void ptDestroy(Point_t *p){
	free(p);
}

/**************************************/
/* Vertex *******************************/
/**************************************/

Vertex* vtxInit(Vertex* v, int id, double x, double y, double z)
{
	v->id= id;
	v->neighborhood = lstCreate();
	v->coords = ptCreate(x,y,z);
}

Vertex* vtxCreate(int id, double x, double y, double z)
{
	Vertex* v = malloc(sizeof(Vertex));
	vtxInit(v, id, x, y, z);
	return v;
}

void vtxAddNeighboord(Vertex* v, int id_v)
{
	lstAddElement(v->neighborhood,id_v);
}

void vtxRemoveNeighboord(Vertex* v, int id_v)
{
	lstRemoveElement(v->neighborhood,id_v);
}

void vtxGetSizeNeighboord(Vertex* v)
{
	lstNbElements(v->neighborhood);
}

void vtxDestroy(Vertex* v)
{
	lstDelete(v->neighborhood);
	ptDestroy(v->coords);
	free(v);
}

Graph_t* gInit(Graph_t* g)
{
	g->vertices = NULL;
	g->size = 0;
	g->capacity = 0;
}

Graph_t* gCreate()
{
	Graph_t* g = malloc(sizeof(Graph_t));
	gInit(g);
	return g;
}

void gAddAlloc(Graph_t* g)
{
	g->vertices = realloc(g->vertices,(g->capacity + REALLOCSIZE)*sizeof(Vertex));
	for (int i =0; i<REALLOCSIZE; i++)
		g->vertices[i] = vtxCreate(-1,0,0,0);
	g->capacity += REALLOCSIZE;
}

int gGetIndiceFree(Graph_t* g)
{
	int i;
	for (i = 0; i<g->size; i++)
	{
		if(g->vertices[i]->id == -1)
			return i;
	}
	gAddAlloc(g);
	return i;
}

void gAddVertex(Graph_t* g, Vertex* v)
{
	if(gsize(g)+1 > g->capacity)
	{
		gAddAlloc(g);
	}
	int indice = gGetIndiceFree(g);
	g->vertices[indice] = v;
	g->size += 1;
}

void gRemoveVertex(Graph_t* g,Vertex* v,int id)
{
	for (int i = 0; i<gsize(g); i++)
	{
		if(v_id(g->vertices[i]) == id)
		{
			g->vertices[i]->id = -1;
			g->size -= 1;
			break;
		}
	}
}

Vertex* gGetVertex(Graph_t* g,int id)
{
	Vertex* v = NULL;
	for(int i = 0; i<gsize(g);i++)
	{
		if(v_id(g->vertices[i]) == id)
			{
				v = g->vertices[i];
				break;
			}
		}
	return v;
}

void gAddEdge(Graph_t* g, int id_v1, int id_v2)
{
	Vertex* v1 = gGetVertex(g, id_v1);
	Vertex* v2 = gGetVertex(g, id_v2);
	
	vtxAddNeighboord(v1, id_v2);
	vtxAddNeighboord(v2,id_v1);
}

void gRemoveEdge(Graph_t* g, int id_v1, int id_v2)
{
	Vertex* v1 = gGetVertex(g, id_v1);
	Vertex* v2 = gGetVertex(g, id_v2);
	
	vtxRemoveNeighboord(v1, id_v2);
	vtxRemoveNeighboord(v2,id_v1);
}

void gDestroy(Graph_t* g)
{
	free(g->vertices);
	free(g);
}