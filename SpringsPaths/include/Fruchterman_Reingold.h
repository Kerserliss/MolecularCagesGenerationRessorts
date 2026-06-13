#ifndef _FRUCHTEMAN_REINGOLD_H
#define _FRUCHTEMAN_REINGOLD_H

double Repulsion_force(double x,double k);
double Attraction_Force(double x,double k);

void ptPrint(Point_t p);
void Add_Path(Cage_t* s);
double CoolingFunction(double temperature);
void ComputeEdgeMat(Cage_t*s,int n, int m,int **mat);
void Fruchterman_Reingold(Cage_t* s);

#endif