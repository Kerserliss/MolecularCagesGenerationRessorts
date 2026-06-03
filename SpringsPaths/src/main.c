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
	FILE *file = fopen("graph.txt", "r");
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

void main()
{
	Graph_t* g = Parse_file();
	int comp[] = {1,2,3,4,5,6};
	List_t* seq1 = Dijkstra(g,1,2,comp);
	lstPrint(seq1);
}