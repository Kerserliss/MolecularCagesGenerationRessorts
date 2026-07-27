# SpringPath

This repository contain the work on generating moleculars path between binding patterns with forces and springs.
This is supposed to work with PathsGeneration, but it was not test.

## Repository Layout
- **Include** and **source** contain all the code
- **run_demos** permit a full run on all demos files with euclidian distance and astar distance . Its take time, take 8 CPU heart ( can modify ).
- **run_plot** scrapt result from a complete run from run_demos and generate differents plots.

### Method
 - **Step 1 :** Load all files and options.
 - **Step 2 :** Find interconnexion tree.
 - **Step 3 :** ( with springpath chosen, else refere to the read me in PathsGeneration) Simulated anneling.
 - **Step 4 :** Produting a molecular cage.
 - **Step 5 :** Save file under the option name chosen.

This work was done by me , during a summer intership at Laboratoire DAVID under the supervision of Noe Demange, Yann Strozeki and Sandrine Vial.
It's based on the work of Marie Bricage, and the thesis of Noe Demange.
