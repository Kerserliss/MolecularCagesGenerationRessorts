#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "structure.h"
#include "Heap.h"
#include "Kamada_Kawai.h"


Graph_t* Parse_file()
{
	FILE *file = fopen("graph2.txt", "r");
    if (file == NULL)
    {
        perror("Error opening file");
    }
    Graph_t* g = gCreate();
    int id;
    int pos[3];
    int neighboords[100];
    char line[256];

    while(fgets(line,sizeof(line),file) != NULL)
    {
    	line[strcspn(line, "\n")] = '\0';
    	Vertex* v;

    	char* f_id = strtok(line,":");
    	char *f_pos = strtok(NULL, ":");
        char *f_nghb = strtok(NULL, ":");

        if (f_id == NULL || f_pos == NULL || f_nghb == NULL)
            perror("Invalid format");

        id = atoi(f_id);
        char* token = strtok(f_pos, ",");
        pos[0] = atoi(token);
        token = strtok(NULL, ",");
        pos[1] = atoi(token);
        token = strtok(NULL, ",");
        pos[2] = atoi(token);

        token = strtok(f_nghb, ",");

        v = vtxCreate(id,pos[0], pos[1], pos[2]);
        while (token != NULL)
        {
            vtxAddNeighboord(v, atoi(token));
            token = strtok(NULL, ",");
        }
        gAddVertex(g,v);
    }
    fclose(file);
    return g;
}

void printMat1(double** mat, int n, int m)
{
    for (int i = 0; i<n;i++)
    {
        for(int j= 0;j<m;j++)
        {
            printf(" %f ",mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printMatInt(int** mat, int n, int m)
{
    for (int i = 0; i<n;i++)
    {
        for(int j= 0;j<m;j++)
        {
            printf(" %d ",mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
void printMat2(double* mat, int n, int m)
{
    for (int i = 0; i<n;i++)
    {
        for(int j= 0;j<m;j++)
        {
            printf(" %f ",mat[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void main()
{
	Graph_t* g=Parse_file();
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
    printMatInt(matDistance,g->size,g->size);
    printMat1(matK,g->size,g->size);
    printMat1(matL,g->size,g->size);
}