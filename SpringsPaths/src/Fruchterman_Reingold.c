#include "structure.h"
#include "substrat.h"
#include "Fruchterman_Reingold.h"
#include "constant.h"
#include "distance.h"
#include "util.h"
#include "output.h"

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
    int compteur = 0;
    for (int i = 0; i<n; i++)
    {
        AtomCage_t* atom_i = atom(s,i);
        if(flag(atom_i) == SPRING_PATH_F)
        {
            compteur+= 1;
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
    }
    RMSD = sqrt(RMSD/compteur);
    return RMSD;
}

double RSMD_Cage_angle(int n, int m, int** mat, Cage_t* s)
{
    double RMSD = 0;
    double angle_c = 0;
    int compteur = 0;
    for (int i = 0; i<n; i++)
    {
        AtomCage_t* atom_b = atom(s,i);
        {
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
                            if(flag(atom_b)==SPRING_PATH_F ||flag(atom_c)==SPRING_PATH_F ||flag(atom_a)==SPRING_PATH_F)
                            {
                                compteur +=1;
                                angle_c = angle(coords(atom_a),coords(atom_b), coords(atom_c));
                                RMSD += pow(fabs(angle_c - 109.5),2);
                            }
                        }
                    }
                }
            }
        }
    }
    RMSD = sqrt(RMSD/compteur);
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


void Add_Path(Cage_t* s, int source, int target, int nb_atom_to_place)
{
	// Recuperation of each atom : source and target.
	AtomCage_t* atom_source = atom(s,source);
	AtomCage_t* atom_target = atom(s,target);

	// Computing vector for normalisation
    Point_t vector_dir = vector(coords(atom_source),coords(atom_target));
    double dist_vector = dist(coords(atom_source),coords(atom_target))/nb_atom_to_place;
    Point_t norm = normalization(vector_dir,dist_vector);

    // Adding atom
    int id_pred = source;
    for (int k =0; k<nb_atom_to_place;k++)
    {
        Point_t coords = ptAdd(coords(atom(s,id_pred)),norm);
        int id_new_atom = cageAddAtom(s,coords);
        flag(atom(s,id_new_atom)) = SPRING_PATH_F;
        cageAddEdge(s,id_new_atom,id_pred);
        id_pred = id_new_atom;
    }
    // Checking if our last atom is collinding with our atom target.
    if(ptEqual(coords(atom(s,id_pred)), coords(atom(s,target))))
    {
        // If yes, we "merge" the last atom and the atom target.
        for (int i = 0; i<cageNbNeighborhood(atom(s,id_pred));i++)
        {
            cageAddEdge(s,target,elts(neighborhood(atom(s,id_pred)),i));
        }
        cageRemoveAtom(s,id_pred);
    }
    else {
        // Else we add an edge between the last atom and the atom target.
        cageAddEdge(s,id_pred,target);
    }
}

double RMSD_angle_path(Cage_t* s, int source, int target)
{
    double RMSD = 0;
    int current = source;
    int id_prev = -1;
    double angle_current = 0;

    while(current!=target)
    {
        for (int i = 0; i<cageNbNeighborhood(atom(s,current)); i++)
        {
            int neigbord1 = elts(atom(s,current)->neighborhood,i);
            if(neigbord1!=id_prev)
            {
                for(int j = 0;j<cageNbNeighborhood(atom(s,neigbord1));j++)
                {
                    int neigbord2 = elts(atom(s,neigbord2)->neighborhood,j);
                    if(neigbord2!= current)
                    {
                        angle_current = angle(coords(atom(s,neigbord1)),coords(atom(s,current)),coords(atom(s,neigbord2)));
                        RMSD += pow(fabs(angle_current - 109.5),2);
                    }
                }
            }
        }
    }
    RMSD = sqrt(RMSD/s->size);
    return RMSD;
}

void ComputeEdgeMat(Cage_t*s,int n, int m,int** mat)
{
    // For each atom
	for(int i = 0; i<n; i++)
	{
		AtomCage_t* Current_atom = atom(s,i);
		printf("i : %d \n",i);
		// We initialize each case of our matrix.
		for(int k = 0; k<m;k++)
		{
		    printf("k : %d \n",k);
			mat[i][k] = -1;
		}
		// For each neighbord, we change the corresponding case, replacing -1 by 1.
		for(int j= 0; j<cageNbNeighborhood(Current_atom);j++)
		{

			int id_neighbor = elts(Current_atom->neighborhood,j);
			mat[i][id_neighbor] = 1;
			AtomCage_t* Neighbor_atom = atom(s,id_neighbor);
			// For each neigborh of our neighbord, we replace -1 by 2.
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

void Fruchterman_Reingold_file(Cage_t* s, double k_attraction1,double k_attraction2, double k_repulsion, double** substrat, int nb_atom_substrat, char* name1)
{
	double RSMD_dist;
	double RSMD_angle;

	FILE* fp1 = fopen(name1,"a");

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
		// Computing the repulsive force intermol and with the substrat.
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
				for (int k = 0; k<nb_atom_substrat; k++)
				{
                Point_t coords_sub = ptInit(0);
                coords_sub.x = substrat[k][0];
                coords_sub.y = substrat[k][1];
                coords_sub.z = substrat[k][2];
				    if (dist(coords_sub, coords(atom(s,j)))<3)
				    {
						delta = vector(coords_sub,coords(atom(s,j)));
						distance = dist(coords_sub,coords(atom(s,j)));
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
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
							}
							else
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
							}
						}
						else
						{
							if (mat[j][k] == 1)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
							}
							else
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
							}
						}
					}
					else if(flag(atom_u) == SPRING_PATH_F)
					{
						delta = vector(coords(atom_u),coords(atom_v));
						distance = dist(coords(atom_u),coords(atom_v));
							if (mat[j][k] == 1)
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
							}
							else
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
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
				Point_t calcul = ptMul(normalization(Vector_disp_list[j],1),fmin(temperature,vectorLength(Vector_disp_list[j])));

				coords(Current_atom) = ptAdd(coords(Current_atom),calcul);
			}
		}
		temperature = CoolingFunction(temperature);
		i += 1;

	}
	printf("I : %d , Temp : %f \n",i,temperature);
	RSMD_dist = RSMD_Cage_dist(s->size, s->size, mat, s);
	RSMD_angle = RSMD_Cage_angle(s->size, s->size, mat, s);

	fprintf(fp1, "%f;%f;%f\n",k_attraction2,RSMD_dist,RSMD_angle);
	fclose(fp1);
	printf("RMSD dist : %f \n",RSMD_dist);
	printf("RMSD angle : %f \n", RSMD_angle);

}

void Fruchterman_Reingold(Cage_t* s, double k_attraction1,double k_attraction2, double k_repulsion, GridSubstrat* grid_t, double step, double* mat_RMSD,int** edge_mat,Options_t options)
{
    // if(options.verbose)
    // {
    //     printf("In Fruchterman\n");
    // }

    // Initialization of our variable.
	double RMSD_dist;
	double RMSD_angle;

	double temperature = TEMPERATURE;
	Point_t* Vector_disp_list = malloc(s->size*sizeof(Point_t));
	Point_t delta;
	double distance;
	int i = 0;


	while(i<ITERATIONS && temperature >0.0001)
	{
	    if(i%10000 == 0 && options.verbose)
	    {
			// printf("Iterations : %d \n",i);
			// printf("Temperature : %f\n",temperature);
		}
		// We initialize our array of vector at 0.
		for (int j = 0; j<s->size; j++)
		{
			Vector_disp_list[j] = ptInit(0);
		}
		// Computing the repulsive force intermol and with the substrat.
		for(int j = 0; j<s->size; j++)
		{
			if(flag(atom(s,j))==SPRING_PATH_F)
			{
                AtomCage_t* atom_j = atom(s,j);
				for(int k = 0; k<s->size;k++)
				{
				    AtomCage_t* atom_k = atom(s,k);
					// If our atom is too close, we compute repulsive force.
					if (k != j&&dist(coords(atom_k),coords(atom_j))<3)
					{

						delta = vector(coords(atom_k),coords(atom_j));
						distance = dist(coords(atom_k),coords(atom_j));
						Vector_disp_list[j] = ptAdd(Vector_disp_list[j],ptMul(normalization(delta,1),Repulsion_force(fabs(distance),k_repulsion)));
						// if (Vector_disp_list[j].x != Vector_disp_list[j].x ||Vector_disp_list[j].y != Vector_disp_list[j].y ||Vector_disp_list[j].z != Vector_disp_list[j].z )
						// {
      //             		    if (j == 25)
      //             		    {
      //                           printf("Coords between atom\n");
      //                           printf("Atom %d \n",j);
      //                           ptPrint(coords(atom_j));
      //                           printf("Atom %d \n",k);
      //                           ptPrint(coords(atom_k));
      //            			}
						// 	if (i% 1000 == 0)
						// 	{
						// 	    printf("Nan in repulsion between vector \n");
						// 		ptPrint(Vector_disp_list[j]);
						// 		printf("Iteration :  %d  Atom 1  : %d Atom 2 : %d dist : %f\n",i,j,k,distance);
						// 		printf("Delta \n");
						// 		ptPrint(delta);
						// 		cageWriteMol2_Spring("Path_atom.mol2",s);
						// 	}
						// }
					}
				}
				// We check collision with the substrat.
				double* collision_adresss = checkGridCollisionSubstratPointTReturnPoint(coords(atom(s,j)), grid_t, step);
				if (collision_adresss != NULL)
				{
                    Point_t coords_sub = ptInit(0);
                    coords_sub.x = collision_adresss[0];
                    coords_sub.y = collision_adresss[1];
                    coords_sub.z = collision_adresss[2];

                    delta = vector(coords_sub,coords(atom(s,j)));
					distance = dist(coords_sub,coords(atom(s,j)));
					Vector_disp_list[j] = ptAdd(Vector_disp_list[j],ptMul(normalization(delta,1),Repulsion_force(fabs(distance),k_repulsion)));
					// if (Vector_disp_list[j].x != Vector_disp_list[j].x ||Vector_disp_list[j].y != Vector_disp_list[j].y ||Vector_disp_list[j].z != Vector_disp_list[j].z )
					// {
					// 	if (i% 1000 == 0)
					// 	{
					// 	    printf("Nan in repulsion between substrat \n");
					// 	}
					// }
				}
			}
		}


		// Computing the attractive force between neighbord.
		for(int j = 0; j<s->size;j++)
		{
			AtomCage_t* atom_v = atom(s,j);
			for (int k = j; k<s->size;k++)
			{
			    // If we have a link with one of the atom, we compute our force.
				if(edge_mat[j][k] != -1)
				{
					AtomCage_t* atom_u = atom(s,k);
					if (flag(atom_v)==SPRING_PATH_F) // We check if the atom who have a link with the first one is allowed to move.
					{
					    // If it's the case, we can compute force.
						delta = vector(coords(atom_u),coords(atom_v));
						distance = dist(coords(atom_u),coords(atom_v));
						if (flag(atom_u)==SPRING_PATH_F) // If the first can also move, we add force on both of the atom with different parameters deping if mat[j][k] = 1 or 2.
						{
							if (edge_mat[j][k] == 1&&distance > 1.5)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
							}
							else if (edge_mat[j][k] == 2 && distance > 2.44)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
							}
						}
						else // if the first can not move, we compute only on the second.
						{
							if (edge_mat[j][k] == 1 && distance > 1.5)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
							}
							else if (edge_mat[j][k] == 2 && distance > 2.44)
							{
								Vector_disp_list[j] = ptSub(Vector_disp_list[j],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
							}
						}
					}
					else if(flag(atom_u) == SPRING_PATH_F) // And if it's the first who can only move, we add force on only the first.
					{
						delta = vector(coords(atom_u),coords(atom_v));
						distance = dist(coords(atom_u),coords(atom_v));
							if (edge_mat[j][k] == 1 && distance > 1.5)
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction1)));
							}
							else if (edge_mat[j][k] == 2 && distance > 2.44)
							{
								Vector_disp_list[k] = ptAdd(Vector_disp_list[k],ptMul(normalization(delta,1),Attraction_Force_quatradic(fabs(distance),k_attraction2)));
							}
					}
					// if (Vector_disp_list[j].x != Vector_disp_list[j].x ||Vector_disp_list[j].y != Vector_disp_list[j].y ||Vector_disp_list[j].z != Vector_disp_list[j].z )
					// {
					// 	if (i% 1000 == 0)
					// 	{
					// 	    printf("Nan in attraction between vector \n");
					// 	}
					// }
				}
			}
		}
		// Point_t max = ptInit(0);
		Point_t pt0 = ptInit(0);
		// Now we add our force to each atom.
		for (int j = 0; j<s->size; j++)
		{
		    AtomCage_t* Current_atom = atom(s,j);
			// if (coords(Current_atom).x != coords(Current_atom).x || coords(Current_atom).y != coords(Current_atom).y || coords(Current_atom).z != coords(Current_atom).z)
			// {
			//     if (i %1000 == 0)
			//         printf("Nan , iteration : %d, atom : %d \n",i,j);
			// }
			// If our displacement vector is null, it's mean that our atom cannot move, so we don't add.
			if(!ptEqual(pt0,Vector_disp_list[j]))
			{
				Point_t calcul = ptMul(normalization(Vector_disp_list[j],1),fmin(temperature,vectorLength(Vector_disp_list[j])));
				// if (calcul.x != calcul.x || calcul.y != calcul.y || calcul.z != calcul.z)
				// {
				//     if (i %1000 == 0)
				//     {
				// 		ptPrint(Vector_disp_list[j]);
    //                     printf("Nan calcul : ");
    //                     ptPrint(calcul);
				// 	}
				// }

				coords(Current_atom) = ptAdd(coords(Current_atom),calcul);
			}

		}
		// We cool our temperature.
		temperature = CoolingFunction(temperature);
		i += 1;

	}
	//printf("I : %d , Temp : %f \n",i,temperature);
	if(options.verbose)
        printf("Finished computing \n");

	// We compute the RMSD on our molecular cage and we stock it in our matrix.
	RMSD_dist = RSMD_Cage_dist(s->size, s->size, edge_mat, s);
	RMSD_angle = RSMD_Cage_angle(s->size, s->size, edge_mat, s);
	mat_RMSD[0] = RMSD_dist;
	mat_RMSD[1] = RMSD_angle;

	//printf("RMSD dist : %f \n",RMSD_dist);
	//printf("RMSD angle : %f \n", RMSD_angle);
}
