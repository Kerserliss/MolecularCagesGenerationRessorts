#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "main.h"
#include "structure.h"

/**
 * @file output.h
 * @brief Output Header File
 *
 * This file contains function for handling various output-related operations.
 */

extern int cpt_result_cage_egv; // Extern global variable to count the number of cages generated

/**
 * @brief Create a directory with the given input name in the "../results" folder.
 *
 * @param input The input filename.
 * @return A dynamically allocated string representing the created directory name.
 */
char *createDir(char *);

/**
 * @brief Write the contents of a list to the console.
 *
 * @param l The list to be written.
 */
void lstWrite(List_t *);

/**
 * @brief Write the information of a cage to the console.
 *
 * @param s The cage to be written.
 */
void cageWrite(Cage_t *);

/**
 * @brief Writes the contents of a cage structure to a .mol2 file format.
 *
 * This function generates a .mol2 file representing the molecular structure of the cage, including atoms,
 * paths, and bonds. The .mol2 format is commonly used for molecular modeling and visualization.
 *
 * @param output The path to the output .mol2 file.
 * @param s A pointer to the `Cage_t` structure representing the molecular system.
 * @param paths A pointer to the `Paths_t` structure containing path information (e.g., atom paths and their positions).
 * @param interTree An array representing the interconnection tree, where each pair of integers defines the start
 *                  and end points of a connection.
 * @param nb_atoms_paths The number of additional atoms in paths to be included in the file.
 * @param cumul_RMSD The cumulative Normalized Root Mean Squared Deviation (RMSD) for all paths, used for analysis.
 *
 * @details
 * ### Function Workflow:
 * 1. **File Initialization**:
 *    - The function opens the specified file for writing. If the file cannot be opened, an error message is printed,
 *      and the function exits.
 *
 * 2. **Molecule Section**:
 *    - Writes the general header of the .mol2 file, including the number of atoms and bonds.
 *    - Assumes the "SMALL" and "GASTEIGER" molecule types in the .mol2 format.
 *
 * 3. **Atom Section**:
 *    - Iterates over all vertices in the cage and writes atom information such as ID, type, and coordinates.
 *    - Determines the type of each atom based on its flag (`CYCLE_F`, `HYDRO_PATTERN_F`, etc.).
 *    - For atoms in paths, their positions are derived from the `Paths_t` structure, and additional atoms (e.g.,
 * hydrogens) are handled specifically.
 *
 * 4. **Bond Section**:
 *    - Iterates over the neighborhood of each atom in the cage to define bonds between connected atoms.
 *    - Includes bonds for paths, connecting the atoms defined in the `interTree` array and additional atoms in
 * `Paths_t`.
 *
 * 5. **Error Handling**:
 *    - If any writing operation fails, an error message is printed, and the program exits with an error code.
 *
 * 6. **Memory Management**:
 *    - Allocates memory for indexing cage vertices to atom IDs.
 *    - Ensures all dynamically allocated memory (e.g., `index`) is freed before the function exits.
 *
 * ### Notes:
 * - The function assumes that all vertices and paths in the input structures are properly initialized.
 * - The atom types (e.g., `H`, `C`, `O`, etc.) and flags (`CYCLE_F`, `HYDRO_PATTERN_F`) must conform to the expected
 * format.
 * - The generated file is in the .mol2 format, compatible with most molecular modeling tools.
 *
 * @see CAGE_nbAtom
 * @see CAGE_nbEdges
 * @see Paths_t
 */
void cageWriteMol2(char *output, Cage_t *s, Paths_t *paths, int *interTree, int nbAtomPaths, double cumul_MSD);

/**
 * @brief Writes output files for a given cage structure .
 *
 * This function generates and writes a `.mol2` file for the given cage structure,
 * including all atoms and edges. It's also write the number of collisions and the execution time fir the cage.
 * It updates atom flags as needed for writing and restores them afterward.
 *
 * @param cage A pointer to the `Cage_t` structure representing the molecular system.
 * @param options The `Options_t` structure containing input file name, output directory,
 *                and the maximum number of results to be written.
 * @param sp A pointer to the molecular cage generated with associated RMSD
 * @param start Time of start of the exectution, from the main
 * @param collision Number of time the moc collision with the substrat ( TODO : check how it's done )
 *
 * @details
 * ### Function Workflow:
 * 1. **Prepare Output Directory**:
 *    - Creates a file based on the ouput file name.
 *
 * 2. **Write `.mol2` File**:
 *    - Constructs the output file name using the input file name and result index.
 *    - Write the differents variables ( time, collision ).
 *    - Write atom and edge given their flag.
 *
 */
void cageWriteMol2_Spring(char *output, SpringPath_t* sp,time_t start, int collision, Options_t options);

/**
 * @brief Write the output files for a given cage.
 *
 * @param options The options with input file name and number of results.
 * @param s The cage to be written.
 * @param paths The paths of the cage.
 */
void writeCageOutput(Cage_t *cage, Paths_t *paths, int *interTree, Options_t options);

/**
 * @brief Write the execution time and statistics to a file.
 *
 * This function calculates the execution time of the program and writes it to a specified output file.
 * It also includes the number of interconnection trees generated during the execution.
 *
 * @param outputname The name of the output file where the execution time and statistics will be written.
 *
 * @details
 * ### Function Workflow:
 * 1. **Open File**:
 *    - Opens the specified file in write mode. If the file cannot be opened, an error message is printed.
 *
 * 2. **Calculate Execution Time**:
 *    - Calculates the total execution time in seconds using `difftime`.
 *    - Calculates CPU time used using `clock()`.
 *
 * 3. **Write Execution Time**:
 *    - Writes the formatted execution time to the file, including hours, minutes, and seconds.
 *
 * 4. **Write Interconnection Trees**:
 *    - Writes the number of interconnection trees generated during execution to the file.
 *
 * 5. **Close File**:
 *    - Closes the file after writing all information.
 *
 */
void writeTime(const char *outputname);

/**
 * @brief Write the parameters of execution.
 *
 * @param options Options_t type.
 */
void writeParameters(Options_t options);

/**
 * @brief Emit partial statistics when the solver receives a termination signal.
 */
void flushStatsOnSignal(void);

#ifdef ENABLE_STATS
/**
 * @brief print the stats of execution.
 */
void printStats();

/**
 * @brief write the stats of execution.
 *
 * @param outpuname the path of the file to print.
 */
void writeStats(const char *outputname);
#endif

#endif
