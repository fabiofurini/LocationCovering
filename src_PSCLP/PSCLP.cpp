

#include "PSCLP.h"


//#define write_prob

//#define solve_LP

//#define print_solution
//#define print_solution_lp




int CPXPUBLIC mycutcallback_CFL_FAKE(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
{

	(*useraction_p)=CPX_CALLBACK_DEFAULT;

	return 0;
}


/***********************************************************************************/
int position_y_CFL(instance *inst,int location)
/***********************************************************************************/
{
	return location;
}

/***********************************************************************************/
int position_z_CFL(instance *inst,int client)
/***********************************************************************************/
{
	return inst->n_locations + client;
}

/*****************************************************************/
void build_model_CFL(instance *inst)
/*****************************************************************/
{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	inst->env_CFL=CPXopenCPLEX(&(inst->status));
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	inst->lp_CFL=CPXcreateprob(inst->env_CFL,&(inst->status),"CFL");
	if(inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables *
	inst->ccnt=inst->n_locations+inst->n_clients;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->c_type=(char*) calloc(inst->ccnt,sizeof(char));


	inst->colname=(char**) calloc(inst->ccnt,sizeof(char*));
	for(int i=0;i<inst->ccnt;i++){inst->colname[i]=(char*) calloc(2000,sizeof(char));}

	int counter=0;
	for ( int j = 0; j < inst->n_locations; j++){

		inst->obj[counter]=inst->fixed_cost[j];
		inst->lb[counter]=0.0;
		inst->ub[counter]=1.0;
		inst->c_type[counter]='B';
		//cout << "POSITION Y\t" << position_y_CFL(inst,j) << endl;
		sprintf(inst->colname[counter], "y%d",j);
		counter++;

	}

	cout << "\n*** CONTINUOUS Z VARIABLES\n";
	for ( int i = 0; i < inst->n_clients; i++ ){

		inst->obj[counter]=0.0;
		inst->lb[counter]=0.0;
		inst->ub[counter]=1.0;
		inst->c_type[counter]='C';
		//cout << "POSITION Z\t" << position_z_CFL(inst,i) << endl;
		sprintf(inst->colname[counter], "z%d",i);
		counter++;
	}

	inst->status=CPXnewcols(inst->env_CFL,inst->lp_CFL,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->c_type,inst->colname);
	if(inst->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	for(int i=0;i<inst->ccnt;i++){free(inst->colname[i]);}
	free(inst->colname);


	// * setting the objective function in the minimization form
	CPXchgobjsen(inst->env_CFL,inst->lp_CFL,CPX_MIN);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	for ( int i = 0; i < inst->n_clients; i++ )
	{
		// * creating the covering  constraint *
		inst->rcnt=1;
		inst->nzcnt=inst->DM[i]+1;

		//		cout << "size\t" << inst->DM[i] << "\tfrom\t" << inst->NBS[i] << "\tto\t" <<  inst->NBS[i+1] << endl;
		//		cin.get();

		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

		inst->rhs[0]=0.0;
		inst->sense[0]='G';


		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		int counter=0;
		for (int  k = inst->NBS[i]; k < inst->NBS[i+1]; k++ )
		{
			//			cout << "location\t" << inst->ABS[k] << endl;
			inst->rmatval[counter]=1.0;
			inst->rmatind[counter++]=position_y_CFL(inst,inst->ABS[k]);
			//			cout << "location\t" << inst->ABS[k]  << endl;
		}


		inst->rmatval[counter]=-1.0;
		inst->rmatind[counter++]=position_z_CFL(inst,i);

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_CFL,inst->lp_CFL,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
		if(inst->status!=0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	// * creating the demand constraint *
	inst->rcnt=1;
	inst->nzcnt=inst->n_clients;
	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

	inst->rhs[0]=inst->COVERING_DEMAND;
	inst->sense[0]='G';


	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	for ( int i = 0; i < inst->n_clients; i++ )
	{
		inst->rmatval[i]=inst->demand[i];
		inst->rmatind[i]=position_z_CFL(inst,i);

	}

	inst->rmatbeg[0]=0;

	inst->status=CPXaddrows(inst->env_CFL,inst->lp_CFL,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
	if(inst->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * writing the created ILP model on a file *
	//	char dummy_name[1000];
	//	sprintf(dummy_name,"./LPS_PSCP/PSCPn%dm%dR%.2fD%.2fs%d.lp",inst->n_locations,inst->n_clients,inst->param_sparsity,inst->param_demand,inst->seed);
	//	inst->status=CPXwriteprob(inst->env_CFL,inst->lp_CFL,dummy_name,NULL);
	//	if(inst->status!=0)
	//	{
	//		printf("error in CPXwriteprob\n");
	//		exit(-1);
	//	}
	//	exit(-1);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status=CPXwriteprob(inst->env_CFL,inst->lp_CFL,"CFL.lp",NULL);
	if(inst->status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif



}



/*****************************************************************/
void solve_model_CFL(instance *inst)
/*****************************************************************/
{



	CPXsetintparam (inst->env_CFL, CPX_PARAM_SCRIND, CPX_ON);


	//	// * Set relative tolerance *
	//	inst->status = CPXsetdblparam (inst->env_CFL, CPX_PARAM_EPAGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}
	//
	//	// * Set a tolerance *
	//	inst->status = CPXsetdblparam (inst->env_CFL, CPX_PARAM_EPGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPGAP\n");
	//	}
	//
	//	// * Set mip tolerances integrality *
	//	inst->status = CPXsetdblparam (inst->env_CFL, CPX_PARAM_EPINT, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//
	//	// * Set Feasibility tolerance *
	//	inst->status = CPXsetdblparam (inst->env_CFL, CPX_PARAM_EPRHS, 1e-9);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPRHS\n");
	//	}

	// * Set number of CPU*
	inst->status = CPXsetintparam (inst->env_CFL, CPX_PARAM_THREADS, inst->number_of_CPU);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	// * Set time limit *
	inst->status = CPXsetdblparam (inst->env_CFL, CPX_PARAM_TILIM,inst->timelimit);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	//CPXsetintparam(inst->env_CFL, CPX_PARAM_PREIND, CPX_OFF);

	//	if(inst->option==1000)
	//	{
	//		solve_LP_CFL(inst);
	//	}


	//	if(inst->option==2)
	//	{
	//
	//		//this is the only only one necessary to avoid the removal of all continuous variables
	//		CPXsetintparam(inst->env_CFL, CPX_PARAM_PREIND, CPX_OFF);
	//
	//
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_AGGIND, CPX_OFF);
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_BNDSTRENIND, CPX_OFF);
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_COEREDIND, CPX_OFF);
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_RELAXPREIND, CPX_OFF);
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_REDUCE, CPX_OFF);
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_PREPASS, CPX_OFF);
	//		//		CPXsetintparam(inst->env_CFL, CPX_PARAM_REPEATPRESOLVE, CPX_OFF);
	//
	//		cout << "***********\n\n AUTOMATIC BENDER'S DECOMPOSITION\n\n";
	//
	//		inst->status = CPXsetintparam (inst->env_CFL, CPXPARAM_Benders_Strategy, 3);
	//		if (inst->status)
	//		{
	//			printf ("error for CPX_PARAM_EPRHS\n");
	//		}
	//	}

	//	if(inst->option==1)
	//	{
	//OPEN FAKE CALLBACK FOR A FAIR COMPARISON
	//		cout << "OPEN FAKE CALLBACK FOR FAIR COMPARISONS\n";
	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_CFL,mycutcallback_CFL_FAKE,inst);
	if (inst->status)
	{
		printf ("error for CPXsetlazyconstraintcallbackfunc\n");
	}
	//	}


	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();

	cout << "\nCPXmipopt:\n";
	inst->status=CPXmipopt(inst->env_CFL,inst->lp_CFL);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");
		//exit(-1);
	}

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	///////////////////////////////////////////////////////////////////////////////////


	bool sol_found=true;

	// * getting the solution
	inst->x=(double*) calloc(inst->n_locations+inst->n_clients,sizeof(double));


	inst->status=CPXgetmipx(inst->env_CFL,inst->lp_CFL,inst->x,0,inst->n_locations+inst->n_clients-1);
	if(inst->status!=0)
	{
		sol_found=false;
		printf("error in CPXgetmipx\n");
	}

	inst->objval=-1;
	inst->status=CPXgetmipobjval(inst->env_CFL,inst->lp_CFL,&(inst->objval));
	if(inst->status!=0)
	{
		sol_found=false;
		printf("error in CPXgetmipobjval\n");
	}

	printf("\n\nMIP solution value ->\t\%f",inst->objval);


	int open_facilities=-1;
	int satisfied_clients=-1;

	if(sol_found){

		open_facilities=0;
		satisfied_clients=0;

		for ( int j = 0; j < inst->n_locations; j++){

			if( (int)(inst->x[position_y_CFL(inst,j)]+0.5) ==1){
				open_facilities++;
			}

		}

		for ( int i = 0; i < inst->n_clients; i++ ){

			if( (int)(inst->x[position_z_CFL(inst,i)]+0.5)==1){
				satisfied_clients++;
			}
		}




	}

#ifdef print_solution
	printf("\n\nSolution\n");
	for ( int j = 0; j < inst->n_locations; j++){

		cout << "Location\t" << j << "\tval:\t" << (int)(inst->x[position_y_CFL(inst,j)]+0.5) << endl;

	}
	printf("\n");
	for ( int i = 0; i < inst->n_clients; i++ ){

		cout << "Client\t" << i << "\tval:\t" << (int)(inst->x[position_z_CFL(inst,i)]+0.5) << endl;
	}
	printf("\n");
#endif


	inst->bestobjval=-1;
	inst->status=CPXgetbestobjval(inst->env_CFL,inst->lp_CFL,&(inst->bestobjval));
	if(inst->status!=0)
	{
		printf("error in CPXgetbestobjval\n");
	}

	inst->lpstat=CPXgetstat(inst->env_CFL,inst->lp_CFL);
	inst->nodecount = CPXgetnodecnt(inst->env_CFL, inst->lp_CFL);

	cout << "\n\nlpstat\t" << inst->lpstat << endl;

	cout << "\n***open_facilities\t" << open_facilities << endl;
	cout << "***satisfied_clients\t" << satisfied_clients << endl;


	///////////////////////////////////////////////////////////////////////////////////
	/* linear programming relaxation*/

	CPXsetintparam (inst->env_CFL, CPX_PARAM_SCRIND, CPX_OFF);

	double solution_time_lp=0;
	double cplex_lp=-1;



	///////////////////////////////////////////////////////////////////////////////////

	int cur_numcols=CPXgetnumcols(inst->env_CFL,inst->lp_CFL);
	int cur_numrows=CPXgetnumrows(inst->env_CFL,inst->lp_CFL);

	cout << "\n\nSTAT:\tobjval\t" << inst->objval << "\tbestobjval\t" << inst->bestobjval << "\tlpstat\t" << inst->lpstat << "\topen_facilities\t" << open_facilities << "\tsatisfied_clients\t" << satisfied_clients << "\ttime\t"<< solution_time<< endl << endl;


	//	printf("\nnumcols\t%d\n",cur_numcols);
	//	printf("\nnumrows\t%d\n",cur_numrows);

	//	ofstream compact_file;
	//	compact_file.open("info_CFL.txt", ios::app);
	//	compact_file << fixed
	//			<< inst->input_file << "\t"
	//			<< inst->n_locations << "\t"
	//			<< inst->n_clients << "\t"
	//			<< inst->RADIUS << "\t"
	//			<< inst->COVERING_DEMAND << "\t"
	//			<<  inst->objval<< "\t"
	//			<<  inst->bestobjval<< "\t"
	//			<<  inst->lpstat<< "\t"
	//			<<   solution_time << "\t"
	//			<<  inst->nodecount<< "\t"
	//			<<  cplex_lp<< "\t"
	//			<<  solution_time_lp<< "\t"
	//			<<  open_facilities<< "\t"
	//			<<  satisfied_clients<< "\t"
	//			<<  inst->algorithm<< "\t"
	////			<<  inst->option<< "\t"
	//			<<  inst->seed<< "\t"
	//			<<  inst->client_not_covered<< "\t"
	//			<<  inst->client_single_covered<< "\t"
	//			<< endl;
	//	compact_file.close();


	free(inst->x);




}


/*****************************************************************/
void clean_model_CFL(instance *inst)
/*****************************************************************/
{

	inst->status=CPXfreeprob(inst->env_CFL,&(inst->lp_CFL));
	if(inst->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	inst->status=CPXcloseCPLEX(&(inst->env_CFL));
	if(inst->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

}
