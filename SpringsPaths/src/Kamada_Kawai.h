#ifndef __KAMADA_KAWAI_H
#define __KAMADA_KAWAI_H

#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constant.h"

int** InitMat(int line, int col);
List_t* Dijkstra(Graph_t *g, int source, int target, int* correspondance);

#endif