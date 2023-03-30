

#include "global_functions.h"

#define load_forward_star

//#define cover_uncovered_clients

//#define print_point
//#define print_neighbourhoods
//#define print_cost
//#define print_demands
//#define print_scenarios
//#define print_distances

//#define print_stars





/*****************************************************************/
double compute_single_location_coverage(instance *inst, int location)
/*****************************************************************/
{
	double total_demand=0;

	for (int  k = inst->NFS[location]; k < inst->NFS[location+1]; k++ )
	{
		total_demand+=inst->demand[inst->AFS[k]];
	}

	return total_demand;

}

/*****************************************************************/
double compute_subset_location_coverage(instance *inst,double *y)
/*****************************************************************/
{

	double total_demand=0;

	bool *reached_clients=new bool[inst->n_clients];

	for ( int j = 0; j < inst->n_clients; j++){
		reached_clients[j]=false;
	}

	for ( int j = 0; j < inst->n_locations; j++){

		if(y[j]< 0.5){continue;}

		for (int  k = inst->NFS[j]; k < inst->NFS[j+1]; k++ )
		{
			if(reached_clients[inst->AFS[k]]==false)
			{
				total_demand+=inst->demand[inst->AFS[k]];
				reached_clients[inst->AFS[k]]=true;
			}
		}
	}

	delete []reached_clients;

	return total_demand;
}

/*****************************************************************/
int compute_subset_location_coverage_clients(instance *inst,double *y,bool *reached_clients)
/*****************************************************************/
{

	int number=0;

	for ( int j = 0; j < inst->n_clients; j++){
		reached_clients[j]=false;
	}

	for ( int j = 0; j < inst->n_locations; j++){

		if(y[j]< 0.5){continue;}

		for (int  k = inst->NFS[j]; k < inst->NFS[j+1]; k++ )
		{
			if(reached_clients[inst->AFS[k]]==false){
				reached_clients[inst->AFS[k]]=true;
				number++;
			}
		}
	}

	return number;
}


/*****************************************************************/
void compute_location_coverage(instance *inst)
/*****************************************************************/
{

	double *location_coverage_ratio=new double[inst->n_locations];
	double *location_coverage=new double[inst->n_locations];

	for ( int j = 0; j < inst->n_locations; j++){

		//	cout << "LOCATION\t"<< j << "\tCLIENTS\t" ;
		location_coverage[j]=0;
		location_coverage_ratio[j]=0;
		for (int  k = inst->NFS[j]; k < inst->NFS[j+1]; k++ )
		{
			//		cout << inst->AFS[k] << "\t";
			location_coverage[j]+=inst->demand[inst->AFS[k]];
		}
		//	cout << endl;
		location_coverage_ratio[j]=inst->fixed_cost[j]/location_coverage[j];
	}


	//	cout << "COVERAGE\n";
	//	for ( int j = 0; j < inst->n_locations; j++){
	//		cout << "location\t" << j << "\tfixed cost\t" << inst->fixed_cost[j] << "\tlocation_coverage\t" << location_coverage[j] << "\tRATIO\t" << location_coverage_ratio[j] << endl;
	//	}

	delete []location_coverage_ratio;
	delete []location_coverage;

	//exit(-1);
}

/*****************************************************************/
void build_data_structure(instance *inst,vector < vector < int > > dummy)
/*****************************************************************/
{



#ifdef print_cost
	cout << "COSTS:\n";
	for ( int j = 0; j < inst->n_locations; j++ )
	{
		cout << inst->fixed_cost[j] << "\n";

	}
	cout << endl;
#endif


#ifdef print_demands
	cout << "DEMANDS:\n";
	for ( int i = 0; i < inst->n_clients; i++ ){
		cout << inst->demand[i] << "\n";
	}
	cout << endl;
#endif



#ifdef print_neighbourhoods
	for ( int i = 0; i < inst->n_clients; i++ ){

		cout << "Client\t" << i << "\t number of facilities \t " << dummy[i].size() << "\t\t";

		for ( int j = 0; j <  dummy[i].size(); j++){
			cout << dummy[i][j] << " ";
		}

		cout << endl;
	}
	cout << endl;
#endif


	inst->cliets_OK=new bool[inst->n_clients];

	//////////////////////////////////////////////
	for ( int i = 0; i < inst->n_clients; i++ )
	{
		inst->cliets_OK[i]=true;
		if((int)(dummy[i].size())==0)
		{
			inst->cliets_OK[i]=false;
		}
	}
	//////////////////////////////////////////////


	inst->client_not_covered=0;
	inst->client_single_covered=0;
	for ( int i = 0; i < inst->n_clients; i++ ){
		if((int)(dummy[i].size())==0) {
			inst->client_not_covered++;
		}
		if((int)(dummy[i].size())==1) {
			inst->client_single_covered++;
		}
	}
	cout << "client_not_covered\t" << inst->client_not_covered << endl;
	cout << "client_single_covered\t" << inst->client_single_covered << endl;
	//	cin.get();

#ifdef cover_uncovered_clients
	//check if all the clients have at least a location
	bool cliets_OK=true;
	for ( int i = 0; i < inst->n_clients; i++ ){
		if((int)(dummy[i].size())==0) {
			int facility=randomBETWEEN(0,inst->n_locations);
			cout << "clients\t" << i << "\tadding\t" <<  facility << endl;
			dummy[i].push_back(facility);
			cliets_OK=false;
		}
	}
	if(cliets_OK){
		cout << "***CONNECTIONS ARE OK\n";
	}
#endif

	inst->NFS=new int[inst->n_locations+1];
	inst->AFS=new int[MAX_CONNECTIONS];
	inst->DP=new int[inst->n_locations];
	inst->NBS=new int[inst->n_clients+1];
	inst->ABS=new int[MAX_CONNECTIONS];
	inst->DM=new int[inst->n_clients];

	inst->counter_c=0;

	inst->NBS[0]=0;
	for ( int i = 0; i < inst->n_clients; i++ ){

		inst->DM[i]=dummy[i].size();

		inst->NBS[i+1]=inst->NBS[i]+dummy[i].size();

		for ( int j = 0; j <  dummy[i].size(); j++){
			//			cout << counter_c << "\t" << dummy[i][j] << endl;
			//			cin.get();
			inst->ABS[inst->counter_c++]=dummy[i][j];
		}
	}

	cout << "\n***REAL_CONNECTIONS\t" << inst->counter_c << endl;
	//cin.get();


#ifdef 	load_forward_star

	cout << "\n***BUILDING FORWARD STARS\t" << endl;
	inst->counter_l=0;

	inst->NFS[0]=0;
	for ( int j = 0; j < inst->n_locations; j++)
	{

		int local_counter=0;

		for ( int i = 0; i < inst->n_clients; i++ )
		{
			bool found=false;
			for ( int jj = 0; jj <  dummy[i].size() && !found; jj++)
			{
				if(j==dummy[i][jj])
				{
					inst->AFS[inst->counter_l++]=i;
					local_counter++;
					found=true;
				}
			}
		}


		inst->DP[j]=local_counter;

		inst->NFS[j+1]=inst->NFS[j]+local_counter;
	}

//	cout << "\n***REAL_CONNECTIONS\t" << inst->counter_l << endl;

	compute_location_coverage(inst);

#endif


#ifdef print_stars

	for(int i=0;i<inst->n_clients;i++){
		cout << "Backward star of client \t" << i << "\t size of his neighbourhood\t" << inst->DM[i] << endl;
		for (int  k = inst->NBS[i]; k < inst->NBS[i+1]; k++ )
		{
			cout << "location\t" << inst->ABS[k]  << endl;
		}
		cin.get();
	}

	for ( int j = 0; j < inst->n_locations; j++){
		cout << "Forward star of location \t" << j << "\t size of his neighbourhood\t" << inst->DP[j] << endl;
		for (int  k = inst->NFS[j]; k < inst->NFS[j+1]; k++ )
		{
			cout << "client\t" << inst->AFS[k] << endl;
		}
		cin.get();
	}


#endif





}




/*****************************************************************/
double distance(double x1,double y1,double x2, double y2)
/*****************************************************************/
{


	double distancex = (x2 - x1)*(x2 - x1);
	double distancey = (y2 - y1)*(y2 - y1);

	return   sqrt(distancex + distancey);
}

/*****************************************************************/
double distance_manatthan(double x1,double y1,double x2, double y2)
/*****************************************************************/
{


	double distancex;
	if(x2 >= x1){
		distancex=x2 - x1;
	}
	else{
		distancex=x1 - x2;
	}

	double distancey;
	if(y2 >= y1){
		distancey=y2 - y1;
	}
	else{
		distancey=y1 - y2;
	}

	return   (distancex + distancey);
}



/*****************************************************************/
void READ_NEW_FILE(instance *inst)
/*****************************************************************/
{

	cout << "INSTANCE\t" << inst->input_file << endl;


	ifstream in(inst->input_file);
	if(!in)
	{
		ofstream err("Error.log",ios::app);
		cout << "File could not be opened. " << endl;
		exit(1);
	}

	in >> inst->n_locations;
	in >> inst->n_clients;

	cout << "\n\n***LOCATIONS\t" << inst->n_locations << endl;
	cout << "***CLIENTS\t" << inst->n_clients << endl;

	inst->fixed_cost = (int *) calloc(inst->n_locations, sizeof(int));
	inst->demand = (int *) calloc(inst->n_clients, sizeof(int));

	inst->x_location= (double *) calloc(inst->n_locations, sizeof(double));
	inst->y_location= (double *) calloc(inst->n_locations, sizeof(double));

	inst->x_client= (double *) calloc(inst->n_clients, sizeof(double));
	inst->y_client= (double *) calloc(inst->n_clients, sizeof(double));

	inst->cohordinates_loaded=true;


	for ( int j = 0; j < inst->n_locations; j++ )
	{
		char dummy_char;
		int dummy_int;

		in >>dummy_char;
		in >>dummy_int;

//		cout << dummy_char << "\t" <<dummy_int << endl;
//		cin.get();

		in >> inst->x_location[j];
		in >> inst->y_location[j];
		in >> inst->fixed_cost[j];
	}

	for ( int i = 0; i < inst->n_clients; i++ )
	{

		char dummy_char;
		int dummy_int;

		in >>dummy_char;
		in >>dummy_int;

//		cout << dummy_char << "\t" <<dummy_int << endl;
//		cin.get();

		in >> inst->x_client[i];
		in >> inst->y_client[i];
		in >> inst->demand[i];
	}


#ifdef	print_point
	for ( int j = 0; j < inst->n_locations; j++ )
	{
		cout << "Location\t" <<  j << "\t(x,y)\t" << inst->x_location[j] << "\t" << inst->y_location[j] << endl;
	}

	for ( int i = 0; i < inst->n_clients; i++ )
	{

		cout << "Client\t" << i << "\t(x,y)\t" << inst->x_client[i]  << "\t" << inst->y_client[i] << endl;
	}
#endif



#ifdef	print_distances
	cout << "DISTANCES\n";
	for ( int j = 0; j < inst->n_locations; j++ )
	{
		cout << "location\t" << j << endl;
		for ( int i = 0; i < inst->n_clients; i++ )
		{
			cout << "distance to client\t" << i << "\t" << distance_manatthan(inst->x_location[j],inst->y_location[j],inst->x_client[i],inst->y_client[i]) << "\t" << distance(inst->x_location[j],inst->y_location[j],inst->x_client[i],inst->y_client[i])<< endl;

		}
	}
	cin.get();
#endif



	cout << "\n\nBUILDING neighbourhoods\n";
	vector < vector < int > > dummy;
	for ( int i = 0; i < inst->n_clients; i++ ){

		vector < int > local_dummy;

		for ( int j = 0; j < inst->n_locations; j++)
		{

			if(distance(inst->x_client[i],inst->y_client[i],inst->x_location[j],inst->y_location[j]) < inst->RADIUS)
			{
				local_dummy.push_back(j);
			}
		}
		dummy.push_back(local_dummy);
	}

	///////////////////////////////////
	build_data_structure(inst,dummy);
	///////////////////////////////////

}




/*****************************************************************/
void free_data(instance *inst)
/*****************************************************************/
{

	free(inst->fixed_cost);
	free(inst->demand);


	delete []inst->NFS;
	delete []inst->AFS;
	delete []inst->DP;
	delete []inst->NBS;
	delete []inst->ABS;
	delete []inst->DM;



	if(inst->cohordinates_loaded==true){
		free(inst->x_location);
		free(inst->y_location);
		free(inst->x_client);
		free(inst->y_client);
	}

	delete[]inst->cliets_OK;

}

/*****************************************************************/
void  compute_super_rho_DFL(instance *inst,double *DFL_super_rho,double *DFL_SET)
/*****************************************************************/
{

	for (int i = 0; i < inst->n_locations; i++)
	{
		DFL_SET[i]=1.0;
	}

	double F_I=compute_subset_location_coverage(inst,DFL_SET);

	//	cout << "F_I\t" << F_I <<  endl;

	for (int i = 0; i < inst->n_locations; i++)
	{


		DFL_SET[i]=0.0;

		DFL_super_rho[i]=F_I-compute_subset_location_coverage(inst,DFL_SET);

		DFL_SET[i]=1.0;

	}


	//	cout << "***cut_DFL_super_rho\n";
	//	for (int i = 0; i < inst->n_locations; i++)
	//	{
	//		cout << "location\t" << i <<  "\tsuper_rho\t" << DFL_super_rho[i]  << endl;
	//	}
	//
	//	cin.get();

}

/*****************************************************************/
void  compute_super_rho_CFL(instance *inst,double *CFL_super_rho,double *CFL_SET)
/*****************************************************************/
{

	for (int i = 0; i < inst->n_locations; i++)
	{
		CFL_SET[i]=1.0;
	}

	double F_I=compute_subset_location_coverage(inst,CFL_SET);

	//	cout << "F_I\t" << F_I <<  endl;

	for (int i = 0; i < inst->n_locations; i++)
	{


		CFL_SET[i]=0.0;

		CFL_super_rho[i]=F_I-compute_subset_location_coverage(inst,CFL_SET);

		CFL_SET[i]=1.0;

	}


//		cout << "***cut_CFL_super_rho\n";
//		for (int i = 0; i < inst->n_locations; i++)
//		{
//			cout << "location\t" << i <<  "\tsuper_rho\t" << CFL_super_rho[i]  << endl;
//		}
//
//		cin.get();
}

/*****************************************************************/
void  compute_single_rho_CFL(instance *inst,double *CFL_single_rho,double *CFL_SET)
/*****************************************************************/
{

	for (int i = 0; i < inst->n_locations; i++)
	{
		CFL_SET[i]=0.0;
	}



	for (int i = 0; i < inst->n_locations; i++)
	{

		CFL_SET[i]=1.0;

		CFL_single_rho[i]=compute_subset_location_coverage(inst,CFL_SET);

		CFL_SET[i]=0.0;

	}


	//	cout << "***cut_CFL_single_rho\n";
	//	for (int i = 0; i < inst->n_locations; i++)
	//	{
	//		cout << "location\t" << i <<  "\tsingle_rho\t" << CFL_single_rho[i]  << endl;
	//	}
	//
	//	cin.get();
}


/*****************************************************************/
void  compute_single_rho_DFL(instance *inst,double *DFL_single_rho,double *DFL_SET)
/*****************************************************************/
{

	for (int i = 0; i < inst->n_locations; i++)
	{
		DFL_SET[i]=0.0;
	}



	for (int i = 0; i < inst->n_locations; i++)
	{

		DFL_SET[i]=1.0;

		DFL_single_rho[i]=compute_subset_location_coverage(inst,DFL_SET);

		DFL_SET[i]=0.0;

	}


	//	cout << "***cut_DFL_single_rho\n";
	//	for (int i = 0; i < inst->n_locations; i++)
	//	{
	//		cout << "location\t" << i <<  "\tsingle_rho\t" << DFL_single_rho[i]  << endl;
	//	}
	//
	//	cin.get();
}

/*****************************************************************/
void load_I_TILDE_CFL(instance *inst, bool rounding,double *I_TILDE_CFL,double *CFL_BEN_2_Y)
/*****************************************************************/
{


	for ( int j = 0; j < inst->n_clients; j++){

		I_TILDE_CFL[j]=0.0;

		for (int  k = inst->NBS[j]; k < inst->NBS[j+1]; k++ )
		{

			I_TILDE_CFL[j]=I_TILDE_CFL[j]+CFL_BEN_2_Y[inst->ABS[k]];
		}

		if(rounding){
			I_TILDE_CFL[j]=(int)(I_TILDE_CFL[j]+0.0001);
		}
	}


	//	cout << "\n\nI_TILDE_BEN_2\n";
	//	for ( int j = 0; j < inst->n_clients; j++){
	//		cout << "location\t" << j << "\t"<<I_TILDE_CFL[j] << endl;
	//	}
	//	cin.get();

}

/*****************************************************************/
void load_I_TILDE_DFL(instance *inst, bool rounding,double *I_TILDE_DFL,double *DFL_BEN_2_Y)
/*****************************************************************/
{


	for ( int j = 0; j < inst->n_clients; j++){

		I_TILDE_DFL[j]=0.0;

		for (int  k = inst->NBS[j]; k < inst->NBS[j+1]; k++ )
		{

			I_TILDE_DFL[j]=I_TILDE_DFL[j]+DFL_BEN_2_Y[inst->ABS[k]];
		}

		if(rounding){
			I_TILDE_DFL[j]=(int)(I_TILDE_DFL[j]+0.0001);
		}
	}


	//				cout << "\n\nI_TILDE_BEN_2_DFL\n";
	//				for ( int j = 0; j < inst->n_clients; j++){
	//					cout << "location\t" << j << "\t"<<I_TILDE_DFL[j] << endl;
	//				}
	//				cin.get();

}


/*****************************************************************/
void comb_solve_model_CFL_BEN_2_AUX_1(instance *inst,double *I_TILDE,double *AUX_SOL)
/*****************************************************************/
{


	for ( int i = 0; i < inst->n_clients; i++ )
	{

		if(I_TILDE[i] < 1 - 0.0001 ||  inst->DM[i]==1)
		{
			//pi
			AUX_SOL[i]=inst->demand[i];
			//sigma
			AUX_SOL[i+inst->n_clients]=0.0;
		}
		else
		{
			//pi
			AUX_SOL[i]=0.0;
			//sigma
			AUX_SOL[i+inst->n_clients]=inst->demand[i];
		}
	}

}


/*****************************************************************/
void comb_solve_model_CFL_BEN_2_AUX_2(instance *inst,double *I_TILDE,double *AUX_SOL)
/*****************************************************************/
{


	for ( int i = 0; i < inst->n_clients; i++ )
	{


		if(I_TILDE[i] <= 1 + 0.0001)
		{
			//pi
			AUX_SOL[i]=inst->demand[i];
			//sigma
			AUX_SOL[i+inst->n_clients]=0.0;
		}
		else
		{
			//pi
			AUX_SOL[i]=0.0;
			//sigma
			AUX_SOL[i+inst->n_clients]=inst->demand[i];
		}
	}

}

/*****************************************************************/
void comb_solve_model_DFL_BEN_2_AUX_1(instance *inst,double *I_TILDE,double *AUX_SOL)
/*****************************************************************/
{



	for ( int i = 0; i < inst->n_clients; i++ )
	{

		if(I_TILDE[i] < 1 - 0.0001 ||  inst->DM[i]==1)
		{

			//pi
			AUX_SOL[i]=inst->demand[i];
			//sigma
			AUX_SOL[i+inst->n_clients]=0.0;
		}
		else
		{
			//pi
			AUX_SOL[i]=0.0;
			//sigma
			AUX_SOL[i+inst->n_clients]=inst->demand[i];
		}
	}


}

/*****************************************************************/
void comb_solve_model_DFL_BEN_2_AUX_2(instance *inst,double *I_TILDE,double *AUX_SOL)
/*****************************************************************/
{



	for ( int i = 0; i < inst->n_clients; i++ )
	{

		if(I_TILDE[i] <= 1 + 0.0001)
		{

			//pi
			AUX_SOL[i]=inst->demand[i];
			//sigma
			AUX_SOL[i+inst->n_clients]=0.0;
		}
		else
		{
			//pi
			AUX_SOL[i]=0.0;
			//sigma
			AUX_SOL[i+inst->n_clients]=inst->demand[i];
		}
	}


}





