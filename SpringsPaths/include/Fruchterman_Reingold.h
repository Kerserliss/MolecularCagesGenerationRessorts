#ifndef _FRUCHTEMAN_REINGOLD_H
#define _FRUCHTEMAN_REINGOLD_H

#include "structure.h"
double Repulsion_force(double x,double k);
double Attraction_Force(double x,double k);
double RSMD_Cage_dist(int n, int m, int** mat, Cage_t* s);
double RSMD_Cage_angle(int n, int m, int** mat, Cage_t* s);
void ptPrint(Point_t p);
double ptSum(Point_t p);
void Add_Path(Cage_t* s);
double LinearCoolingFunction(double temperature);
double CoolingFunction(double temperature);
void ComputeEdgeMat(Cage_t*s,int n, int m,int **mat);
void Fruchterman_Reingold(Cage_t* s);

#endif
