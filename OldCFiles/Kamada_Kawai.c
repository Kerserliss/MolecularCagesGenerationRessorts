#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <lapacke.h> // For matrix inversion - Version : 3.12.1
#include "structure.h"
#include "Heap.h"
#include "constant.h"
#include "Kamada_Kawai.h"

int** InitMatInt(int line, int col)
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

double** InitMatDouble(int line, int col)
{
	double** m = malloc(line*sizeof(double*));
	for (int i = 0; i < line; i++)
	{
		m[i] = malloc(col*sizeof(double));
		 for (int j = 0; j < col; j++)
		 {
		 	m[i][j] = 0;
		 }
	}
	return m;
}

int GetTabIndice(int* tab, int id, int n)
{
	for(int i = 0; i<n; i++)
	{
		if(tab[i] == id)
		{
			return i;
		}
	}
	return -1;
}

void rdInit(ResultDijkstra* rd)
{
	rd->mat=NULL;
	rd->correspondance=NULL;
}

ResultDijkstra* rdCreate()
{
	ResultDijkstra* rd = malloc(sizeof(ResultDijkstra));
	rdInit(rd);
	return rd;
}
void rdDestroy(ResultDijkstra* rd)
{
	free(rd->correspondance);
	free(rd);
}

double euclidien_distance(Point_t* p1, Point_t* p2)
{
	double dx = p1->x - p2->x;
	double dy = p1->y - p2->y;
	double dz = p1->z - p2->z;
	return sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
}

ResultDijkstra* Dijkstra(Graph_t *g, int source)
{
	// Matrice initialization. Heap initialization.
	int mdist[g->size];
	int mprev[g->size];
	Heap* h = hpCreate();
	int* correspondance = malloc(g->size*sizeof(int));

	// We parcourt our matrice and set it at -1 for previous and g->size for dist 
	for(int i = 0; i<g->size; i++)
	{	
		// We get the id at the i-eme indice
		int vertex_id = v_id(g->vertices[i]);
		correspondance[i] = vertex_id; // We add it to our correspondance tab.

		Tuple* t = tCreate();
		t-> nb1 = vertex_id;
		// We initiliaze both matrice
		mprev[i] = -1;
		if(source == vertex_id) // If it's our source, we set at 0.
		{
			mdist[i]=0;
			t-> nb2 = 0;
		}
		else // If not, we initiliaze it a the maximun distance.
		{
			mdist[i]= g->size +1;
			t->nb2 = g->size +1;
		}
		hpAddElement(h,t); // And we add our element to our heap.
	}
	
	printf("Heap size before algo : %d \n",h->size);

	for(int j = 0; j<g->size; j++)
	{
		printf("Mprev indice %d : %d \n",j, correspondance[j]);
	}

	Tuple* t_min = hpExtractMin(h);
	// Start of the algorithm
	while (h->size >0 && t_min->nb1!= -1)
	{
		printf("Heap size during algo : %d \n",h->size);
		printf("Id in the tuple : %d \n", t_min->nb1);
		Vertex* v = gGetVertex(g,t_min->nb1);
		printf("Treated vertex : %d \n",v->id);
		for (int i = 0; i<size(v->neighborhood) && elts(v->neighborhood,i) != -1; i++)
		{	
			lstPrint(v->neighborhood);
			int id_neighbor = elts(v->neighborhood,i);
			printf("Treated Neighbord : %d \n",id_neighbor);
			int temp = mdist[GetTabIndice(correspondance,t_min->nb1, g->size)] + 1; // +1 car distance tjr de 1
			printf("Temp : %d \n",temp);
			if (temp < mdist[GetTabIndice(correspondance,id_neighbor, g->size)])
			{
				printf("Heap size before increase Val : %d \n",h->size);
				hpIncreaseVal(h,id_neighbor,temp);
				printf("Heap size after increase Val : %d \n",h->size);
				mdist[GetTabIndice(correspondance,id_neighbor,g->size)] = temp;
				mprev[GetTabIndice(correspondance,id_neighbor, g->size)] = v_id(v);
				for(int j = 0; j<g->size; j++)
				{
					printf("Mprev indice %d : %d \n",j, mprev[j]);
				}
			}
			printHeap(h);
		}
		t_min = hpExtractMin(h);
		printHeap(h);
	}

	ResultDijkstra* rd = rdCreate();
	rd-> mat = mdist;
	rd->correspondance = correspondance;
	return rd;

/*	// Find path
	List_t* seq = lstCreate();
	lstAddElement(seq,target);
	int size = 0;
	int u = GetTabIndice(correspondance,target,g->size);
	for(int j=0; j<g->size;j++)
	{
		printf("Mat prev Id : %d Pred : %d \n",j, mprev[j]);
	}
	int source_corr = GetTabIndice(correspondance,source,g->size);
	printf("Source Corr : %d \n",source_corr);
	while(mprev[u] != -1 && u != source_corr)
	{
		printf("U : %d \n",u);
		printf("Prev U corr : %d \n",GetTabIndice(correspondance,mprev[u],g->size));
		lstAddElement(seq,u);
		u = GetTabIndice(correspondance,mprev[u],g->size);

	}
	lstAddElement(seq,u);

	// Get it in the right sway
	reverseList(seq);
	return seq;*/

}	

double* Hessienne(Graph_t* g,Vertex* v,double** matK, double** matL)
{
	double Ex2;
	double Ey2;
	double Ez2;
	double Exy;
	double Exz;
	double Eyz;
	int indice_v = gGetVertexIndice(g,v_id(v));
	for (int i =0; i < g->size; i++)
	{
		if(indice_v == i) // We ignore if it the same vertex
		{
			continue;
		}

		// Values Initialization
		Vertex* v_i = g->vertices[i];
		double euclidien_d = euclidien_distance(coords(v),coords(g->vertices[i]));
		double dx = coords(v)->x - coords(v_i)->x;
		double dy = coords(v)->y - coords(v_i)->y;
		double dz = coords(v)->z - coords(v_i)->z;

		// Compute diagonale of hessian matrice
		Ex2 += matK[indice_v][i]*(1-((matL[indice_v][i]*(pow(dy,2)+pow(dz,2)))/pow(euclidien_d,3)));
		Ey2 += matK[indice_v][i]*(1-((matL[indice_v][i]*(pow(dx,2)+pow(dz,2)))/pow(euclidien_d,3)));
		Ez2 += matK[indice_v][i]*(1-((matL[indice_v][i]*(pow(dx,2)+pow(dy,2)))/pow(euclidien_d,3)));

		// Compute non diagnonale of hessian matrice
		Exy += matK[indice_v][i] *(((dx*dy)*matL[indice_v][i])/pow(euclidien_d,3));
		Exz += matK[indice_v][i] *(((dx*dz)*matL[indice_v][i])/pow(euclidien_d,3));
		Eyz += matK[indice_v][i] *(((dy*dz)*matL[indice_v][i])/pow(euclidien_d,3));
	}


	double mat_temp[3][3] = {{Ex2,Exy,Exz},{Exy,Ey2,Eyz},{Exz,Eyz,Ez2}};
	double* mat = malloc(3*3*sizeof(double));
	memcpy(mat, mat_temp,3*3*sizeof(double));
	return mat;
}

double* Energie_Gradient(Graph_t* g, Vertex* v, double** matK, double** matL)
{
	double Ex;
	double Ey;
	double Ez;
	int indice_v = gGetVertexIndice(g,v_id(v));
	double* matE = malloc(3*sizeof(double));
	for (int i =0; i < g->size; i++)
	{
		if(indice_v == i) // We ignore if it the same vertex
		{
			continue;
		}

		// Values Initialization
		Vertex* v_i = g->vertices[i];
		double euclidien_d = euclidien_distance(coords(v),coords(g->vertices[i]));
		double dx = coords(v)->x - coords(v_i)->x;
		double dy = coords(v)->y - coords(v_i)->y;
		double dz = coords(v)->z - coords(v_i)->z;

		//Compute Energie Gradient
		Ex += matK[indice_v][i]*(1-(matL[indice_v][i]/euclidien_d))*(dx);
		Ey += matK[indice_v][i]*(1-(matL[indice_v][i]/euclidien_d))*(dy);
		Ez += matK[indice_v][i]*(1-(matL[indice_v][i]/euclidien_d))*(dz);
	}
	matE[0] = Ex;
	matE[1] = Ey;
	matE[2] = Ez;
	return matE;
}

double Gradient_Norme(double* matE)
{
	return sqrt(pow(matE[0],2)+pow(matE[1],2)+ pow(matE[2],2));
}

double* MatrixMultiplication(double* Mat1, int line1, int col1, double* Mat2, int line2, int col2)
{
	double* result = malloc(line1*line2*sizeof(double));
	for(int i = 0; i<line1; i++)
	{
		for(int j=0;j<col2; j++)
		{
			result[i*line1+j] =0;
			for(int k=0;k<line2; k++)
			{
				result[i*line1+j] += Mat1[i*line1+k] * Mat2[k*line1+j];
			}
		}
	}
	return result;
}

double* Newton_Raphson(Graph_t* g, Vertex* v, double** matK, double** matL)
{
	double* matNewton;
	double* matEnergie = Energie_Gradient(g,v,matK,matL);
	double* matHessienne = Hessienne(g,v,matK,matL);

	// Inversion of Hessienne Matrice
	lapack_int ipiv[3]; 
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, matHessienne, 3, ipiv); // Decomposition in factor P, L and U
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, matHessienne, 3, ipiv); // Inversion.
	for(int i = 0; i<3;i++)
	{
		for (int j =0; j<1;j++)
		{
			matEnergie[i*3+j]*=-1;
		}
	}

	matNewton = MatrixMultiplication(matHessienne,3,3,matEnergie,3,1);
	return matNewton;


}

void Kamada_Kawai(Graph_t* g)
{
	// Compute Distance between each vertex of the graph.
	int** matDistance = InitMatInt(g->size,g->size);
	for(int i =0; i<g->size; i++)
	{	
		Vertex* v = g->vertices[i];
		ResultDijkstra* rd = Dijkstra(g,v_id(v));
		matDistance[i] = rd->mat;
	}

	// Compute Lij and Kij for each vertices
	double** matK = malloc(g->size*sizeof(double*));
	double** matL = malloc(g->size*sizeof(double*));
	for (int i = 0; i<g->size; i++)
	{
		matK[i] = malloc(g->size*sizeof(double));
		matK[i] = malloc(g->size*sizeof(double));
		for(int j = 0; j<g->size; j++)
		{
			if(i==j)
			{
				continue;
			}
			matK[i][j] = K_RIGIDITY /(matDistance[i][j]*matDistance[i][j]);
			matL[i][j] = matDistance[i][j] * LENGHT;
		}
	}
	
	// Compute the norme for each vectex. Implentation of a max heap


}