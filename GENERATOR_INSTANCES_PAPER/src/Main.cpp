#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <iomanip>

using namespace std;

typedef struct {


	int n_locations;
	int n_clients;
	int *demand;
	int *fixed_cost;
	double param_sparsity;
	double param_demand;
	int seed;


	bool cohordinates_loaded;
	double *x_location;
	double *y_location;

	double *x_client;
	double *y_client;

	int min_demand,max_demand;
	int min_cost,max_cost;

	int size_grid;


} instance;

// return a integer random value in range min-max
/*****************************************************************/
int randomBETWEEN(int min,int max)
/*****************************************************************/
{
	return rand() % (max - min +1) +min;
}

// return random value in range min-max
/*****************************************************************/
double randomBETWEEN_double(int min,int max)
/*****************************************************************/
{
	return (rand()/(double)RAND_MAX)*(max-min) + min;
}

/*****************************************************************/
void create_demands_and_costs(instance *inst)
/*****************************************************************/
{

	for ( int j = 0; j < inst->n_locations; j++ )
	{
		inst->fixed_cost[j]=randomBETWEEN(inst->min_cost,inst->max_cost);
	}

	for ( int i = 0; i < inst->n_clients; i++ ){
		inst->demand[i]=randomBETWEEN(inst->min_demand,inst->max_demand);
	}

}




/*****************************************************************/
void random_generator_grid(instance *inst)
/*****************************************************************/
{


	cout << "***LOCATIONS\t" << inst->n_locations << endl;
	cout << "***CLIENTS\t" << inst->n_clients << endl;

	inst->fixed_cost = (int *) calloc(inst->n_locations, sizeof(int));
	inst->demand = (int *) calloc(inst->n_clients, sizeof(int));

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	inst->x_location= (double *) calloc(inst->n_locations, sizeof(double));
	inst->y_location= (double *) calloc(inst->n_locations, sizeof(double));

	inst->x_client= (double *) calloc(inst->n_clients, sizeof(double));
	inst->y_client= (double *) calloc(inst->n_clients, sizeof(double));

	inst->cohordinates_loaded=true;


	cout << "SAMPLING POINT IN THE GRID\t" << 0 << "\t" << inst->size_grid << endl;

	for ( int j = 0; j < inst->n_locations; j++ )
	{


		double x=randomBETWEEN_double(0,inst->size_grid);
		double y=randomBETWEEN_double(0,inst->size_grid);
		inst->x_location[j]=x;
		inst->y_location[j]=y;

	}

	for ( int i = 0; i < inst->n_clients; i++ )
	{

		double x=randomBETWEEN_double(0,inst->size_grid);
		double y=randomBETWEEN_double(0,inst->size_grid);;
		inst->x_client[i]=x;
		inst->y_client[i]=y;
	}


	///////////////////////////////////
	create_demands_and_costs(inst);
	///////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	char FILE_NAME[10000];
	ofstream output_file;

	sprintf(FILE_NAME,"INSTANCES_NEW/GRID_INST_n%d_m%d_d%d_%d_f%d_%d_s%d.dat",inst->n_locations,inst->n_clients,inst->min_demand,inst->max_demand,inst->min_cost,inst->max_cost,inst->seed);
	cout << "FILE:\t\t" << FILE_NAME << endl;
	output_file.open(FILE_NAME);
	output_file << fixed
			<< inst->n_locations << "\t"
			<< inst->n_clients << "\t"
			<< endl;

	for ( int j = 0; j < inst->n_locations; j++ )
	{
		output_file << "F\t" <<  j << "\t" << inst->x_location[j] << "\t" << inst->y_location[j] << "\t" << inst->fixed_cost[j] << endl;
	}

	for ( int i = 0; i < inst->n_clients; i++ )
	{
		output_file << "C\t" << i << "\t" << inst->x_client[i]  << "\t" << inst->y_client[i] << "\t" << inst->demand[i]<< endl;
	}
	output_file.close();
	exit(-1);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////



}


/**************************************************************/
int main(int argc, char** argv)
/**************************************************************/
{

	instance inst;

	if (argc == 9)
	{
		inst.n_locations=atoi(argv[1]);
		inst.n_clients=atoi(argv[2]);
		inst.min_demand=atoi(argv[3]);
		inst.max_demand=atoi(argv[4]);
		inst.min_cost=atoi(argv[5]);
		inst.max_cost=atoi(argv[6]);
		inst.seed=atoi(argv[7]);
		inst.size_grid=atoi(argv[8]);

		srand(inst.seed);

	}
	else
	{
		cout << "ERROR NUMBER STANDARD PARAMETERS -->>> n";
		exit(-1);
	}


	cout << "n_locations\t" << inst.n_locations << endl;
	cout << "n_clients\t" << inst.n_clients << endl;
	cout << "min_demand\t" << inst.min_demand << endl;
	cout << "max_demand\t" << inst.max_demand << endl;
	cout << "min_cost\t" << inst.min_cost << endl;
	cout << "max_cost\t" << inst.max_cost << endl;
	cout << "seed\t" << inst.seed << endl;

	///////////////////////////////
	random_generator_grid(&inst);
	///////////////////////////////

	printf("\n\nDONE!\n");

	return 1;
}



