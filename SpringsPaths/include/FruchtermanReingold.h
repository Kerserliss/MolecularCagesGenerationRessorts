#ifndef _FRUCHTEMAN_REINGOLD_H
#define _FRUCHTEMAN_REINGOLD_H

#include "main.h"
#include "structure.h"
#include "substrat.h"
#include "interconnection.h"
/**
 * @brief Compute the repulsion force.
 */
double RepulsionForce(double x,double k);

/**
 * @brief Compute the attractive force.
 */
double AttractionForce(double x,double k);

/**
 * @brief Compute the quadratic attractive force.
 */
double AttractionForceQuadratic(double x, double k);

/**
 * @brief Compute the distance RMSD for the entire cage.
 * @param n Number of line in the matrix.
 * @param m Number of colunms in the matrix.
 * @param mat Edge matrix pre-calculated.
 * @param s Molecular cage.
 */
double RMSDCageDist(int n, int m, int** mat, Cage_t* s);

/**
 * @brief Compute the angle RMSD for the entire cage.
 * @param n Number of line in the matrix.
 * @param m Number of colunms in the matrix.
 * @param mat Edge matrix pre-calculated.
 * @param s Molecular cage.
 */
double RMSDCageAngle(int n, int m, int** mat, Cage_t* s);

/**
 * @brief Print a point used for debug.
 * @param p Point to print
 */
void ptPrint(Point_t p);

/**
 * @brief Sum each coordinate, used for debug.
 * @param p Point to sum
 */
double ptSum(Point_t p);

/**
 * @brief Add a linear line composed of atom linked between them between a source atom and a target atom.
 * The number of atom to place is determined before, when the function is called.
 * @param s Molecular cage.
 * @param source Id of the source atom.
 * @param target Id of the target atom.
 * @param nb_atom_to_place Number of atom to place.
 */
void AddPath(Cage_t* s, int source, int target, int nb_atom_to_place);

/**
 * @brief Linear cooling function
 */
double LinearCoolingFunction(double temperature);

/**
 * @brief Exponential cooling function.
 */
double CoolingFunction(double temperature);

/**
 * @brief Compute an edge matrix with 1 if the two atoms are immediate neighbord, 2 if they are at distance 2 in the graph, -1 else.
 * @param s Molecular cage.
 * @param n Number of line
 * @param m Number of columns.
 * @param mat double pointer to the matrix.
 */
void ComputeEdgeMat(Cage_t*s,int n, int m,int **mat);

/**
 * @brief Execute the FruchtermanReingold algorithm on a molecular cage.
 *
 * @details
 * ### Key Steps of the FruchtermanReingold Algorithm :
 * 1. **Initialization**: Initialization of variable and our displacement matrix. This displacement matrix will be reset at each iteration
 * and it will be used to calculate how our atoms will be moved
 * 2. **Force calculation** :
 *    - Atoms can only move if they have the Spring Path flag.
 *    - For the Repulsive Force :
 *          - Calculated only if we are closer than 3 with the euclidian distance.
 *          - If we have a collision with the substrat, calculation too.
 *
 *    - For the Attractive force :
 *          If atoms are too far away, and they have a relation in the edge matrix ( 1 or 2 ), calculation of the atractive force. Else, we do nothing.
 * 3. **Movement**
 *    - We add our displacement nd we cool our temperature.
 *    - Repeat the three step until our temperature is inferior to 0.0001 or i > ITERATION ( constant ).
 *
 * 4. **Ending**
 *    - Calculation of RMSD and stock it in the RMSD matrix.
 *
 * @param s Molecular cage.
 * @param k_attraction 1 Attraction force for immediate neighbord.
 * @param k_attraction 2 Attraction force for distance 2 neighbord in the graph.
 * @param k_repulsion Repulsion force.
 * @param grid_t Grid Substrat
 * @param step Step of the grid substrat ( need removal, because it's in constant )
 * @param mat_RMSD To stock RMSD result.
 * @param edge_mat Edge matrix.
 * @param options Options for verbose.
 */
void FruchtermanReingold(Cage_t* s, double k_attraction1,double k_attraction2, double k_repulsion, GridSubstrat* grid_t, double step, double* mat_RMSD,int** edge_mat,Options_t options);

#endif
