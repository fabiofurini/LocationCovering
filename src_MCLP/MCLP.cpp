

#include "MCLP.h"

//#define solve_LP

//#define write_prob

//#define print_solution
//#define print_solution_lp



int CPXPUBLIC mycutcallback_DFL_FAKE(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
{

	(*useraction_p)=CPX_CALLBACK_DEFAULT;

	return 0;
}


/***********************************************************************************/
int position_y_DFL(instance *inst,int location)
/***********************************************************************************/
{
	return location;
}

/***********************************************************************************/
int position_z_DFL(instance *inst,int client)
/***********************************************************************************/
{
	return inst->n_locations + client;
}

/*****************************************************************/
void build_model_DFL(instance *inst)
/*****************************************************************/
{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	inst->env_DFL=CPXopenCPLEX(&(inst->status));
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	inst->lp_DFL=CPXcreateprob(inst->env_DFL,&(inst->status),"DFL");
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

		inst->obj[counter]=0.0;
		inst->lb[counter]=0.0;
		inst->ub[counter]=1.0;
		inst->c_type[counter]='B';
		//cout << "POSITION Y\t" << position_y_DFL(inst,j) << endl;
		sprintf(inst->colname[counter], "y%d",j);
		counter++;

	}

	cout << "\n*** CONTINUOUS Z VARIABLES\n";
	for ( int i = 0; i < inst->n_clients; i++ ){

		inst->obj[counter]=inst->demand[i];
		inst->lb[counter]=0.0;
		inst->ub[counter]=1.0;

		//		if(!inst->cliets_OK[i])
		//		{
		//			cout << "CLIENT\t" << i << "\tUNCOVERED\n";
		//			inst->ub[counter]=0.0;
		//		}

		inst->c_type[counter]='C';
		//cout << "POSITION Z\t" << position_z_DFL(inst,i) << endl;
		sprintf(inst->colname[counter], "z%d",i);
		counter++;
	}

	inst->status=CPXnewcols(inst->env_DFL,inst->lp_DFL,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->c_type,inst->colname);
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
	CPXchgobjsen(inst->env_DFL,inst->lp_DFL,CPX_MAX);
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
			inst->rmatind[counter++]=position_y_DFL(inst,inst->ABS[k]);
			//			cout << "location\t" << inst->ABS[k]  << endl;
		}


		inst->rmatval[counter]=-1.0;
		inst->rmatind[counter++]=position_z_DFL(inst,i);

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_DFL,inst->lp_DFL,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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


	// * creating the budget constraint *
	inst->rcnt=1;
	inst->nzcnt=inst->n_locations;
	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

	inst->rhs[0]=inst->BUDGET;
	inst->sense[0]='L';


	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	for ( int i = 0; i < inst->n_locations; i++ )
	{
		inst->rmatval[i]=inst->fixed_cost[i];
		inst->rmatind[i]=position_y_DFL(inst,i);

	}

	inst->rmatbeg[0]=0;

	inst->status=CPXaddrows(inst->env_DFL,inst->lp_DFL,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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


#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status=CPXwriteprob(inst->env_DFL,inst->lp_DFL,"DFL.lp",NULL);
	if(inst->status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	exit(-1);
#endif



}



/*****************************************************************/
void solve_model_DFL(instance *inst)
/*****************************************************************/
{



	CPXsetintparam (inst->env_DFL, CPX_PARAM_SCRIND, CPX_ON);


	//	// * Set relative tolerance *
	//	inst->status = CPXsetdblparam (inst->env_DFL, CPX_PARAM_EPAGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}
	//
	//	// * Set a tolerance *
	//	inst->status = CPXsetdblparam (inst->env_DFL, CPX_PARAM_EPGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPGAP\n");
	//	}
	//
	//	// * Set mip tolerances integrality *
	//	inst->status = CPXsetdblparam (inst->env_DFL, CPX_PARAM_EPINT, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//
	//	// * Set Feasibility tolerance *
	//	inst->status = CPXsetdblparam (inst->env_DFL, CPX_PARAM_EPRHS, 1e-9);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPRHS\n");
	//	}

	// * Set number of CPU*
	inst->status = CPXsetintparam (inst->env_DFL, CPX_PARAM_THREADS, inst->number_of_CPU);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	// * Set time limit *
	inst->status = CPXsetdblparam (inst->env_DFL, CPX_PARAM_TILIM,inst->timelimit);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	//	if(inst->option==1000)
	//	{
	//		solve_LP_DFL(inst);
	//	}
	//
	//	if(inst->option==2)
	//	{
	//
	//		//this is the only only one necessary to avoid the removal of all continuous variables
	//		CPXsetintparam(inst->env_DFL, CPX_PARAM_PREIND, CPX_OFF);
	//
	//
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_AGGIND, CPX_OFF);
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_BNDSTRENIND, CPX_OFF);
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_COEREDIND, CPX_OFF);
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_RELAXPREIND, CPX_OFF);
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_REDUCE, CPX_OFF);
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_PREPASS, CPX_OFF);
	////				CPXsetintparam(inst->env_DFL, CPX_PARAM_REPEATPRESOLVE, CPX_OFF);
	//
	//		cout << "***********\n\n AUTOMATIC BENDER'S DECOMPOSITION\n\n";
	//
	//		inst->status = CPXsetintparam (inst->env_DFL, CPXPARAM_Benders_Strategy, 3);
	//		if (inst->status)
	//		{
	//			printf ("error for CPX_PARAM_EPRHS\n");
	//		}
	//	}


	//	if(inst->option==1)
	//	{

	//OPEN FAKE CALLBACK FOR A FAIR COMPARISON WITH BENDERS (SAME BRANCHING STRATEGY!!!!)
	//cout << "OPEN FAKE CALLBACK FOR FAIR COMPARISONS\n";
	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_DFL,mycutcallback_DFL_FAKE,inst);
	if (inst->status)
	{
		printf ("error for CPXsetlazyconstraintcallbackfunc\n");
	}
	//	}

	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();

	cout << "\nCPXmipopt:\n";
	inst->status=CPXmipopt(inst->env_DFL,inst->lp_DFL);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;
	///////////////////////////////////////////////////////////////////////////////////

	bool sol_found=true;

	// * getting the solution

	inst->x=(double*) calloc(inst->n_locations+inst->n_clients,sizeof(double));


	inst->status=CPXgetmipx(inst->env_DFL,inst->lp_DFL,inst->x,0,inst->n_locations+inst->n_clients-1);
	if(inst->status!=0)
	{
		sol_found=false;
		printf("error in CPXgetmipx\n");
	}

	inst->status=CPXgetmipobjval(inst->env_DFL,inst->lp_DFL,&(inst->objval));
	if(inst->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
	}

	printf("\n\nMIP solution value ->\t\%f",inst->objval);

	int open_facilities=-1;
	int satisfied_clients=-1;

	if(sol_found){


		open_facilities=0;
		satisfied_clients=0;

		for ( int j = 0; j < inst->n_locations; j++){

			if( (int)(inst->x[position_y_DFL(inst,j)]+0.5) ==1){
				open_facilities++;
			}

		}

		for ( int i = 0; i < inst->n_clients; i++ ){

			if( (int)(inst->x[position_z_DFL(inst,i)]+0.5)==1){
				satisfied_clients++;
			}
		}

//		/////////////////////////////////////////////////////////////
//		if(inst->cohordinates_loaded==true){
//			int *zz=new int[inst->n_clients];
//			int *yy=new int[inst->n_locations];
//			for ( int j = 0; j < inst->n_locations; j++){
//				yy[j]= (int)(inst->x[position_y_DFL(inst,j)]+0.5);
//			}
//			for ( int i = 0; i < inst->n_clients; i++ ){
//				zz[i]=(int)(inst->x[position_z_DFL(inst,i)]+0.5);
//			}
//
//			draw_grid_sol(inst,yy,zz);
//			delete[]zz;
//			delete[]yy;
//		}
//		/////////////////////////////////////////////////////////////
	}

#ifdef print_solution
	printf("\n\nSolution\n");
	for ( int j = 0; j < inst->n_locations; j++){

		cout << "Location\t" << j << "\tval:\t" << (int)(inst->x[position_y_DFL(inst,j)]+0.5) << endl;

	}
	printf("\n");
	for ( int i = 0; i < inst->n_clients; i++ ){

		cout << "Client\t" << i << "\tval:\t" << (int)(inst->x[position_z_DFL(inst,i)]+0.5) << endl;
	}
	printf("\n");
#endif

	inst->bestobjval=-1;
	inst->status=CPXgetbestobjval(inst->env_DFL,inst->lp_DFL,&(inst->bestobjval));
	if(inst->status!=0)
	{
		printf("error in CPXgetbestobjval\n");
	}

	inst->lpstat=CPXgetstat(inst->env_DFL,inst->lp_DFL);
	inst->nodecount = CPXgetnodecnt(inst->env_DFL, inst->lp_DFL);

	cout << "\n\nlpstat\t" << inst->lpstat << endl;


	cout << "\n***open_facilities\t" << open_facilities << endl;
	cout << "***satisfied_clients\t" << satisfied_clients << endl;


	cout << "\n\nSTAT:\tobjval\t" << inst->objval << "\tbestobjval\t" << inst->bestobjval << "\tlpstat\t" << inst->lpstat << "\topen_facilities\t" << open_facilities << "\tsatisfied_clients\t" << satisfied_clients << "\ttime\t"<< solution_time<< endl << endl;


	int cur_numcols=CPXgetnumcols(inst->env_DFL,inst->lp_DFL);
	int cur_numrows=CPXgetnumrows(inst->env_DFL,inst->lp_DFL);


	//	printf("\nnumcols\t%d\n",cur_numcols);
	//	printf("\nnumrows\t%d\n",cur_numrows);

//	ofstream compact_file;
//	compact_file.open("info_DFL.txt", ios::app);
//	compact_file << fixed
//			<< inst->input_file << "\t"
//			<< inst->n_locations << "\t"
//			<< inst->n_clients << "\t"
//			<< inst->RADIUS << "\t"
//			<< inst->BUDGET << "\t"
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
//			<<  inst->seed<< "\t"
//			<< endl;
//	compact_file.close();


	free(inst->x);




}


/*****************************************************************/
void clean_model_DFL(instance *inst)
/*****************************************************************/
{

	inst->status=CPXfreeprob(inst->env_DFL,&(inst->lp_DFL));
	if(inst->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	inst->status=CPXcloseCPLEX(&(inst->env_DFL));
	if(inst->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

}
