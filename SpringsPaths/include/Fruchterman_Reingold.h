#ifndef _FRUCHTEMAN_REINGOLD_H
#define _FRUCHTEMAN_REINGOLD_H

double Repulsion_force(double x);
double Attraction_Force(double x);

void Add_Path(Cage_t* s);
void Fruchterman_Reinglod(Cage_t* s);

#endif