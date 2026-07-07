#ifndef _FRUCHTEMAN_REINGOLD_H
#define _FRUCHTEMAN_REINGOLD_H

#include "main.h"
#include "structure.h"
#include "substrat.h"
#include "interconnection.h"
double Repulsion_force(double x,double k);
double Attraction_Force(double x,double k);
double Attraction_Force_quatradic(double x, double k);
double RSMD_Cage_dist(int n, int m, int** mat, Cage_t* s);
double RSMD_Cage_angle(int n, int m, int** mat, Cage_t* s);
void ptPrint(Point_t p);
double ptSum(Point_t p);
void Add_Path(Cage_t* s, int source, int target, int nb_atom_to_place);
double LinearCoolingFunction(double temperature);
double CoolingFunction(double temperature);
void ComputeEdgeMat(Cage_t*s,int n, int m,int **mat);
void Fruchterman_Reingold_file(Cage_t* s, double k_attraction1,double k_attraction2, double k_repulsion, double** substrat, int nb_atom_substrat, char* name1);
void Fruchterman_Reingold(Cage_t* s, double k_attraction1,double k_attraction2, double k_repulsion, GridSubstrat* grid_t, double step, double* mat_RMSD,int** edge_mat,Options_t options);

#endif
