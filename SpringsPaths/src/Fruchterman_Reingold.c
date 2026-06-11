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

void ComputeEdgeMat(Cage_t*s,int n, int m,int** mat[n][m])
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

void Fruchterman_Reinglod(Cage_t* s)
{
	Point_t* Point_list_dips = malloc(s->size*sizeof(Point_t));
	**int mat [s->size][s->size] = malloc(s->size*sizeof(int*));
	for(int i = 0; i<s->size;i++)
	{
		mat[i] = malloc(s->size*sizeof(int));
	}

	ComputeEdgeMat(s,s->size,s->size,mat);
	for (int i =0; i<ITERATIONS; i++)
	{
		// Computing the repulsive force.
		for(int j = 0; j<s->size; j++)
		{
			Point_list[j] = ptInit(0);
			if(flag(atom(s,j))==SPRING_PATH_F)
			{			
				for(int k = 0; k<s->size;k++)
				{
					if (k != j)
					{
						Point_t delta = vector(coords(atom(s,k)),coords(atom(s,j)));
						double dist = dist(coords(atom(s,k)),coords(atom(s,j)));
						disp = ptAdd(disp,ptMul(normalization(delta,1),Repulsion_force(fabs(dist),1.8)));
						Point_list[j] = disp;
					}
				}	
			}	
		}
		// Computing the attractive force between neighbord.
		for(int j = 0; j<s->size;j++)
		{	
			if(flag(atom(s,j))==SPRING_PATH_F)
			{
				AtomCage_t* atom_j = atom(s,j);
				for (int k = 0; k<neighborhoodSize(atom_j))
				{	
					AtomCage_t* atom_neighbord = atom(s,elts(atom_j->neighborhood,k));
					delta = vector(coords(atom_neighbord),coords(atom_j));
					dist = dist(coords(atom_neighbord),coords(atom_j));
					Point_list[j] = ptSub(Point_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(dist),1.5)));
					if (flag(atom_neighbord) ==  SPRING_PATH_F)
					{
						Point_list[elts(atom_j->neighborhood,k)] = ptAdd(Point_list[elts(atom_j->neighborhood,k)],ptMul(normalization(delta,1),Attraction_Force(fabs(dist),1.5)));
						for (int l = 0; l<neighborhoodSize(atom_neighbord))
						{
							if (elts(atom_neighbord->neighborhood,l) =! j && atom(s,elts(atom_neighbord->neighborhood,l)) == SPRING_PATH_F)
							{
								AtomCage_t* atom_neighbord_k = atom(s,elts(atom_neighbord->neighborhood,l));
								delta = vector(coords(atom_neighbord_k),coords(atom_neighbord));
								dist = dist(coords(atom_neighbord_k),coords(atom_neighbord));
								Point_list[elts(atom_j->neighborhood,k)] =  
							}
						}
					}
				}
			}
		}
	}	
}