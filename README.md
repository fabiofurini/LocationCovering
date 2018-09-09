# LocationCovering
This site contains accompanying material to the paper:

Benders Decomposition for Very Large Scale Partial Set Covering and Maximal Covering Problems

by Fabio Furini, Ivana Ljubic and Jean-Francois.

Solver for Maximal Covering Location and  partial set covering location problems

The Branch-and-Benders-cut algorithm (described in the publication mentioned above) can be downloaded by clicking here. To use the solver, a license file ..

The source code of the exact algorithm is provided and it needs CPLEX libraries to be compiled

The solver can be run as, e.g.,

./CFLU GRID_MCLP_n100_m10000_d1_100_f1_1_s1.dat 2 600  6.25 10

For the instance format, see below (or the file readme.md in the provided zipfile). The software is for academic purposes only, see also the file license.md in the provided zipfile.

