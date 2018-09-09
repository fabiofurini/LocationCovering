# LocationCovering
This site contains accompanying material to the paper:

> Benders Decomposition for Very Large Scale Partial Set Covering and Maximal Covering Location Problems  
Fabio Furini, Ivana Ljubic and Jean-Francois.


The Branch-and-Benders-cut algorithm for Maximal Covering Location and  partial set covering location problems (described in the publication mentioned above) can be downloaded by clicking here (see the files above). To use the solver, a license file is provided. To compile and run the code  CPLEX libraries are required.

The solver can be run as, e.g.,

./[exeNAME] GRID_MCLP_n100_m10000_d1_100_f1_1_s1.dat 2 600  6.25 10

The parameters are: i) the instance name ii) 1 for the compact model 2 for the Branch-and-Benders-cut iii) the radius of coverage iiii) for the MCLP the Covering Demand and for the PSCLP the Budget

As far as the instance format is concerned: the first line reports the number of potential facility locations and the number of customers. Then for each potential facility location the file reports the cohordinates (x,y) and the  fixed cost. Finally, for each customer the file reports the cohordinates (x,y) and the demand.    

The software is for academic purposes only, see also the file license.md  provided.

