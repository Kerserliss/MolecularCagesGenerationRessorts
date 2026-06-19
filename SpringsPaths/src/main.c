/**
 * @file main.c
 * @brief Main program for generating cages for a specific substrate.
 *
 * This program reads a substrate molecule from an input .moc2 file, reads a
 * partial cage from an input .moc2 and generates whole cages.
 */

#include <getopt.h>
#include <math.h>
#include <signal.h>
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
int main()
{
  int i = 0;
  double k1_a = 1.34;
  double k2_a = 1.22;
  double k_r = 1.7;
  char name1[100];
	sprintf(name1, "Result_k_repulsion_%f_k_attraction1_%f.csv",k_r,k1_a);
  FILE* fp1 = fopen(name1,"a");
  fprintf(fp1,"K_attractivité1,RMSD_dist,RMSD_angle\n");
  fclose(fp1);
  
  // while(i<100)
  // {
  Cage_t* s = cageImport("src/demos_test/TestNoSub1path","0");
  Add_Path(s);
    
  Fruchterman_Reingold(s,k1_a,k2_a,k_r,name1);
  cageDelete(s);
  cageWriteMol2_Spring("./Result_strange.mol2", s);
  //  }
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
