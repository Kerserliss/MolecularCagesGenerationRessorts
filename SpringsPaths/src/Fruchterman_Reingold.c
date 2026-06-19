#include "structure.h"
#include "Fruchterman_Reingold.h"
#include "constant.h"
#include "distance.h"
#include "util.h"
#include "interconnection.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

double Repulsion_force(double x,double k)
{
	return pow(k,2)/x;
}

double RSMD_Cage_dist(int n, int m, int** mat, Cage_t* s)
{
    double RMSD = 0;
    double distance = 0;
    for (int i = 0; i<n; i++)
    {
        AtomCage_t* atom_i = atom(s,i);
        for (int j = i; j<m ;j++)
        {
            if (mat[i][j] ==1)
            {
                AtomCage_t* atom_j = atom(s,j);
                distance = dist(coords(atom_i),coords(atom_j))-1.5;
                RMSD += pow(fabs(distance),2);
            }
        }
    }
    RMSD = sqrt(RMSD/s->size);
    return RMSD;
}

double RSMD_Cage_angle(int n, int m, int** mat, Cage_t* s)
{
    double RMSD = 0;
    double angle_c = 0;
    for (int i = 0; i<n; i++)
    {
        AtomCage_t* atom_b = atom(s,i);
        for (int j = i; j<m ;j++)
        {
            if (mat[i][j] ==2)
            {
                AtomCage_t* atom_c = atom(s,j);
                for (int k = 0; k<cageNbNeighborhood(atom_c);k++)
                {
                    if (lstCheck(atom_b->neighborhood,elts(atom_c->neighborhood,k))) // Find the atom in commun.
                    {
                        AtomCage_t* atom_a = atom(s,elts(atom_c->neighborhood,k));
                        angle_c = angle(coords(atom_a),coords(atom_b), coords(atom_c));
                        RMSD += pow(fabs(angle_c - 109.5),2);
                    }
                }
            }
        }
    }
    RMSD = sqrt(RMSD/s->size);
    return RMSD;
}

double Attraction_Force(double x,double k)
{
	return pow(x,2)/k;
}

double Attraction_Force_quatradic(double x, double k)
{
	return k*log(x);
}

double LinearCoolingFunction(double temperature)
{
    return temperature - COOLING_RATE;
}

double CoolingFunction(double temperature)
{
	return temperature * COOLING_RATE;
}


void Add_Path(Cage_t* s)
{
	// Creation of a list to have all Linkable in a list.
	List_t* flag_list = lstCreate();
	List_t* already_linked = lstCreate();
	for (int i = 0; i<s->size;i++)
	{
	    printf("Id Atom %d flag : %d\n",i,flag(atom(s,i)));
		if(flag(atom(s,i))==LINKABLE_F) // If our linkable have already too much neighbord we don't consider it.
		{
			lstAddElement(flag_list,i);
		}
	}
	// Search connexe componants in the graph.

	for (int i = 0; i< lstNbElements(flag_list);i++) // Same check here.
	{
	    if(!(lstCheck(already_linked,i)))
	    {
    		AtomCage_t* atom_i = atom(s,elts(flag_list,i));

    		double min[2] = {-1,DBL_MAX};

    		for (int j =0; j<lstNbElements(flag_list);j++)
    		{
    			if (i!=j && !(lstCheck(already_linked,j)))
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

    			double dist_atom = Al_kashi_therorem(1.5,1.5,1.911136)/2;
    			Point_t norm = normalization(vector_dir,dist_atom);

    			int nb_atoms_to_place = (int)round(min[1]/dist_atom)-1;

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
                lstAddElement(already_linked, i);
                lstAddElement(already_linked, (int)min[0]);
    		}
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

double ptSum(Point_t p)
{
    double sum = p.x + p.y + p.z;
    return sum;
}

void Fruchterman_Reingold(Cage_t* s, double k_attraction1,double k_attraction2, double k_repulsion,char* name1)
{
	FILE* fp1 = fopen(name1,"a");
	

	double RSMD_dist;
	double RSMD_angle;

	double temperature = TEMPERATURE;
	Point_t* Vector_disp_list = malloc(s->size*sizeof(Point_t));
	Point_t delta;
	double distance;
	double dist_al_kashi = Al_kashi_therorem(1.5,1.5,1.911136)/2;
	int** mat = (int**)malloc(s->size*sizeof(int*));
	for(int i = 0; i<s->size;i++)
	{
		mat[i] = (int*)malloc(s->size*sizeof(int));
	}


	int i = 0;
	ComputeEdgeMat(s,s->size,s->size,mat);
	while(i<ITERATIONS && temperature >0)
	{
	    if(i%10000 == 0)
			printf("Iterations : %d \n",i);
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
					if (k != j&&dist(coords(atom(s,k)),coords(atom(s,j)))<3)
					{
						delta = vector(coords(atom(s,k)),coords(atom(s,j)));
						distance = dist(coords(atom(s,k)),coords(atom(s,j)));
						Vector_disp_list[j] = ptAdd(Vector_disp_list[j],ptMul(normalization(delta,1),Repulsion_force(fabs(distance),k_repulsion)));
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
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction1)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction1)));
							}
							else
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction2)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction2)));
							}
						}
						else
						{
							if (mat[j][k] == 1)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction1)));
							}
							else
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction2)));
							}
						}
					}
					else if(flag(atom_u) == SPRING_PATH_F)
					{
						delta = vector(coords(atom_u),coords(atom_v));
						distance = dist(coords(atom_u),coords(atom_v));
							if (mat[j][k] == 1)
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction1)));
							}
							else
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force(fabs(distance),k_attraction2)));
							}
					}
				}
			}
		}
		// Point_t max = ptInit(0);
		Point_t pt0 = ptInit(0);
		for (int j = 0; j<s->size; j++)
		{
			if(!ptEqual(pt0,Vector_disp_list[j]))
			{
				AtomCage_t* Current_atom = atom(s,j);
				Point_t calcul = ptMul(normalization(Vector_disp_list[j],1),fmin(temperature,dist(ptInit(0),Vector_disp_list[j])));
				coords(Current_atom) = ptAdd(coords(Current_atom),calcul);
			}
		}
		temperature = CoolingFunction(temperature);
		i += 1;

	}
	printf("I : %d , Temp : %f \n",i,temperature);
	RSMD_dist = RSMD_Cage_dist(s->size, s->size, mat, s);
	RSMD_angle = RSMD_Cage_angle(s->size, s->size, mat, s);
	printf("RMSD dist : %f \n",RSMD_dist);
	printf("RMSD angle : %f \n", RSMD_angle);

	fprintf(fp1,"%f,%f,%f\n",k_attraction2,RSMD_dist,RSMD_angle);
	fclose(fp1);

}
