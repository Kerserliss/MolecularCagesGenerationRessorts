/**
 * @file main.c
 * @brief Main program for generating cages for a specific substrate.
 *
 * This program reads a substrate molecule from an input .moc2 file, reads a
 * partial cage from an input .moc2 and generates whole cages.
 */



#include <getopt.h>
#include <math.h>
//#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "constant.h"
#include "distance.h"
#include "interconnection.h"
#include "main.h"
#include "output.h"
#include "structure.h"
#include "substrat.h"
#include "util.h"
#include "Fruchterman_Reingold.h"

// for stats
#include "assembly.h"

struct timespec start_time_ts;
clock_t start_clock;

/**
 * Main function to run the cage generation process.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return Exit status of the program.
 */
int main(int argc, char** argv)
{
    srand((unsigned)time(NULL));
    time_t start = time(&start);
    /******************************** Options *****/
    int opt;
    Options_t options = {NULL,
                         NULL,
                         NULL,
                         DEFLT_SIZEMAX,
                         DEFLT_MAX_RESULTS,
                         DEFLT_BANNED_EDGES,
                         DEFLT_ONE_CAGE_BY_INTERCONNECTION_TREE,
                         DEFLT_PATH_BOUNDARY,
                         DEFLT_DYNAMIC_PATH_LIMIT,
                         DEFLT_SORT_INTERCONNECTION_TREES,
                         DEFLT_SPRING_PATH,
                         DEFLT_VERBOSE};
    options.enablePathBoundary = is_path_boundary_filter_enabled();
    int boundary_cli_override = 0;
    //printf("After default option\n");

    while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
        switch (opt) {
        case 'i':
            options.input = optarg;
            break;
        case 'n':
            options.numMoc = optarg;
            break;
        case 'o':
            options.output = optarg;
            break;
        case 's':
            options.sizeMaxPath = atoi(optarg);
            break;
        case 'r':
            options.maxResults = atoi(optarg);
            break;
        case 'b':
            options.isBannedEdges = atoi(optarg);
            break;
        case 't':
            options.oneCageByInterconnectionTree = atoi(optarg);
            break;
        case 'p':
            options.enablePathBoundary = atoi(optarg);
            set_path_boundary_filter_enabled(options.enablePathBoundary != 0);
            boundary_cli_override = 1;
            break;
        case 'l':
            options.enableDynamicPathLimit = atoi(optarg);
            break;
        case 'g':
            options.sortInterTreesBeforePaths = atoi(optarg);
            break;
        case 'c' :
            options.springPath = atoi(optarg);
            break;
        case 'v' :
            options.verbose = atoi(optarg);
            break;

        case 'h':
        default:
            usage();
            break;

      }
    }
    if (options.output == NULL) {
        options.output = DEFLT_OUTPUTPATH;
    }
  // printf("In main before opening file\n");
  // char name1[100];
  // sprintf(name1, "Result_k_repulsion_%f_k_attraction1_%f_quadratic.csv",k_r,k1_a);
  // FILE* fp1 = fopen(name1,"w");
  // fprintf(fp1,"K_attraction2;RMSD_dist;RMSD_angle\n");
  // fclose(fp1);
  // printf("In main after opening file\n");

  // Chargement du substrat :
  double **substrat_t = NULL;
  GridSubstrat grid_sub = importSubstratToGrid(options.input, &substrat_t);

  Cage_t *import_cage = cageImport(options.input, options.numMoc,options);
  if (options.verbose)
      printf("Loaded\n");
  findInterconnection(import_cage, &grid_sub, &substrat_t, start, options);
  freeGridSubstrat(grid_sub);
  free2DDouble(substrat_t, grid_sub.substratSize);
  cageDelete(import_cage);
}

/**
 * Print the usage of the program.
 */
void usage() {
        fprintf(stderr, USAGE_FMT, DEFLT_SIZEMAX, DEFLT_MAX_RESULTS, DEFLT_BANNED_EDGES,
          DEFLT_ONE_CAGE_BY_INTERCONNECTION_TREE, DEFLT_PATH_BOUNDARY, DEFLT_DYNAMIC_PATH_LIMIT,
          DEFLT_SORT_INTERCONNECTION_TREES);
  exit(EXIT_FAILURE);
}

/**
 * @brief Handles termination signals (SIGTERM/SIGINT).
 *
 * This function is invoked when a timeout signal is received. It prints
 * the values of global variables `cptResultCage_EGV` and `cptInterTree_EGV`
 * (representing the number of results and interconnection trees, respectively)
 * to `stderr`, and then gracefully terminates the program.
 *
 * @param signum The signal number that triggered the handler (e.g., SIGTERM).
 *
 * @note This function uses the `exit(1)` call to terminate the program.
 *       Ensure that any required cleanup has been completed before invoking
 * this handler.
 */
void timeoutHandler(int signum) {
  fprintf(stderr, "Termination signal (%d) received!\n", signum);
  flushStatsOnSignal();
  fprintf(stderr, "Global Variable Values:\n");
  fprintf(stderr, "  Cages Results: %d\n", cpt_result_cage_egv);
  fprintf(stderr, "  Interconnection Trees: %d\n", cpt_inter_tree_egv);
  exit(EXIT_FAILURE); // Exit the program
}
