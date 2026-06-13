#include "structure.h"
#include "Fruchterman_Reingold.h"
#include "constant.h"
#include "distance.h"
#include "util.h"

#include <math.h>
#include <float.h>

double Repulsion_force(double x,double k)
{
	return pow(k,2)/x;
}


double Attraction_Force(double x,double k)
{
	return pow(x,2)/k;
}

double CoolingFunction(double temperature)
{
	return temperature * COOLING_RATE;
}

void Add_Path(Cage_t* s)
{	
	// Creation of a list to have all Linkable in a list.
	List_t* flag_list = lstCreate();
	for (int i = 0; i<s->size;i++)
	{
		if(flag(atom(s,i))==LINKABLE_F && cageNbNeighborhood(atom(s,i))<2) // If our linkable have already too much neighbord we don't consider it.
		{	
			lstAddElement(flag_list,i);
		}
	}
	for (int i = 0; i< lstNbElements(flag_list)&& cageNbNeighborhood(atom(s,elts(flag_list,i)))<2;i++) // Same check here.
	{	
		AtomCage_t* atom_i = atom(s,elts(flag_list,i));

		double min[2] = {-1,DBL_MAX};

		for (int j =0; j<lstNbElements(flag_list)&& cageNbNeighborhood(atom(s,elts(flag_list,j)))<2;j++)
		{	
			if (i!=j)
			{	
				AtomCage_t* atom_j = atom(s,elts(flag_list,j));
				double distance_euclidian = dist(coords(atom_i),coords(atom_j));

					if (min[1]>distance_euclidian)
					{
						min[0] = (double)elts(flag_list,j);
						min[1] = distance_euclidian;
					}
			}
		}
		if (min[0] != -1)
		{		
			AtomCage_t* atom_select = atom(s,(int)min[0]);
			Point_t vector_dir = vector(coords(atom_i),coords(atom_select));

			Point_t norm = normalization(vector_dir,1.22);

			int nb_atoms_to_place = (int)round(min[1]/DIST_SIMPLE_PATTERN)-1;

			int id_pred = elts(flag_list,i);
			for (int k =0; k<nb_atoms_to_place;k++)
				{
					Point_t coords = ptAdd(coords(atom(s,id_pred)),norm);
					int id_new_atom = cageAddAtom(s,coords);
					flag(atom(s,id_new_atom)) = SPRING_PATH_F;
					cageAddEdge(s,id_new_atom,id_pred);
				id_pred = id_new_atom;
				}
			cageAddEdge(s,id_pred,(int)(min[0]));
		}
	}

}

void ComputeEdgeMat(Cage_t*s,int n, int m,int** mat)
{
	for(int i = 0; i<n; i++)
	{
		AtomCage_t* Current_atom = atom(s,i);
		for(int k = 0; k<m;k++)
		{
			mat[i][k] = -1;
		}
		for(int j= 0; j<cageNbNeighborhood(Current_atom);j++)
		{
			int id_neighbor = elts(Current_atom->neighborhood,j);
			mat[i][id_neighbor] = 1;
			AtomCage_t* Neighbor_atom = atom(s,id_neighbor);
			for(int l = 0; l<cageNbNeighborhood(Neighbor_atom);l++)
			{
				int id_neighbor_2 = elts(Neighbor_atom->neighborhood,l);
				if(id_neighbor_2!=i)
				{
					mat[i][id_neighbor_2] = 2;
				}
			}
		}
	}
}

void ptPrint(Point_t p)
{
  printf("\tCoords x : %f\n", p.x);
  printf("\tCoords y : %f\n", p.y);
  printf("\tCoords z : %f\n", p.z);
}

void Fruchterman_Reingold(Cage_t* s)
{
	double temperature = TEMPERATURE;
	Point_t* Vector_disp_list = malloc(s->size*sizeof(Point_t));
	Point_t delta;
	double distance;
	int** mat = (int**)malloc(s->size*sizeof(int*));
	for(int i = 0; i<s->size;i++)
	{
		mat[i] = (int*)malloc(s->size*sizeof(int));
	}

	ComputeEdgeMat(s,s->size,s->size,mat);
	for (int i =0; i<ITERATIONS; i++)
	{
		printf("ITERATIONS : %d \n",i);
		for (int j = 0; j<s->size; j++)
			{
				Vector_disp_list[j] = ptInit(0);
			}
		// Computing the repulsive force.
		for(int j = 0; j<s->size; j++)
		{	
			if(flag(atom(s,j))==SPRING_PATH_F)
			{			
				for(int k = 0; k<s->size;k++)
				{
					if (k != j)
					{
						delta = vector(coords(atom(s,k)),coords(atom(s,j)));
						distance = dist(coords(atom(s,k)),coords(atom(s,j)));
						Vector_disp_list[j] = ptAdd(Vector_disp_list[j],ptMul(normalization(delta,1),Repulsion_force(fabs(distance),1.8)));
					}
				}	
			}	
		}
		// Computing the attractive force between neighbord.
		for(int j = 0; j<s->size;j++)
		{	
			AtomCage_t* atom_v = atom(s,j);
			for (int k = j; k<s->size;k++)
			{	
				if(mat[j][k] != -1)
				{
					AtomCage_t* atom_u = atom(s,k);
					if (flag(atom_v)==SPRING_PATH_F)
					{
						delta = vector(coords(atom_u),coords(atom_v));
						distance = dist(coords(atom_u),coords(atom_v));
						if (flag(atom_u)==SPRING_PATH_F)
						{
							if (mat[j][k] == 1)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),1.5)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),1.5)));
							}
							else
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),2.44)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),2.44)));	
							}
						}
						else
						{
							if (mat[j][k] == 1)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),1.5)));
							}
							else
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),2.44)));
							}
						}
					}
					else if(flag(atom_u) == SPRING_PATH_F)
					{
						delta = vector(coords(atom_u),coords(atom_v));
						distance = dist(coords(atom_u),coords(atom_v));
							if (mat[j][k] == 1)
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),1.5)));
							}
							else
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),2.44)));	
							}
					}
				}		
			}
		}
		printf("Before adding to the atom \n");
		for (int j = 0; j<s->size; j++)
			{
				printf("Vector %d : \n",j);
				ptPrint(Vector_disp_list[j]);
				
			}

		// Adding the displacement to the point.
		for (int j = 0; j<s->size; j++)
		{
			if(!ptEqual(ptInit(0),Vector_disp_list[j]))
			{
				AtomCage_t* Current_atom = atom(s,j);
				Point_t calcul = ptMul(normalization(Vector_disp_list[j],1),fmin(temperature,dist(ptInit(0),Vector_disp_list[j])));
				coords(Current_atom) = ptAdd(coords(Current_atom),calcul);
			}
		}
		for (int j = 0; j<s->size; j++)
			{
				printf("Atom %d : \n",j);
				ptPrint(coords(atom(s,j)));
				
			}
		// Cooling the temperature.
		temperature = CoolingFunction(temperature);
	}	
}