

#include "MCLP_BEN.h"

//#define USE_CUT_LOOP

#define lambda_update 0.9

//#define DEBUG
//#define DEBUG_plus

#define EPSILON_cut_DFL_BEN_2 0.0001
#define MAX_SIZE 10000000

double _cut_DFL_BEN_2_rmatval[MAX_SIZE];
double _cut_DFL_BEN_2_Y[MAX_SIZE];
double _cut_DFL_BEN_2_SET[MAX_SIZE];
int    _cut_DFL_BEN_2_rmatind[MAX_SIZE];
double _cut_DFL_BEN_2_obj_value;
double _cut_DFL_BEN_2_RHS;
double _cut_DFL_BEN_2_value;


#define TOLL_FRAC_BEN 10000

//#define print_solution
//#define print_solution_lp

//#define local_cut_DFL_BEN_2
#define global_cut_DFL_BEN_2
//#define global_cut_DFL_BEN_2_purgeable

//#define write_prob
//#define write_prob_BEN_2_AUX

//#define solve_LP

//#define print_solution
//#define print_solution_lp


double *AUX_SOL_BEN_2_DFL;
double *I_TILDE_BEN_2_DFL;


int n_cuts_BEN_DFL_1=0;
int n_cuts_BEN_FRAC_DFL_1=0;
int n_cuts_BEN_DFL_2=0;
int n_cuts_BEN_FRAC_DFL_2=0;


double _cut_DFL_BEN_2_local_rho[MAX_SIZE];
double _cut_DFL_BEN_2_super_rho[MAX_SIZE];
double _cut_DFL_BEN_2_single_rho[MAX_SIZE];

/*****************************************************************/
int CPXPUBLIC mycutcallback_DFL_BEN_2(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
/*****************************************************************/
{


	(*useraction_p)=CPX_CALLBACK_DEFAULT;

	instance *inst=(instance *) cbhandle;

	int num_variables=inst->n_locations+1;
	int status;
	double nu;
	double val_comparison;

	status=CPXgetcallbacknodex(env,cbdata,wherefrom,_cut_DFL_BEN_2_Y,0,num_variables-1);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}

	status=CPXgetcallbacknodeobjval(env,cbdata,wherefrom,&_cut_DFL_BEN_2_obj_value);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}

	double _BEST_INTEGER;
	status=CPXgetcallbackinfo(env,cbdata,wherefrom,CPX_CALLBACK_INFO_BEST_INTEGER,&_BEST_INTEGER);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}
	double _BEST_REMAINING;
	status=CPXgetcallbackinfo(env,cbdata,wherefrom,CPX_CALLBACK_INFO_BEST_REMAINING,&_BEST_REMAINING);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}
	int _NODE;
	status=CPXgetcallbackinfo(env,cbdata,wherefrom,CPX_CALLBACK_INFO_NODE_COUNT,&_NODE);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}


	////////////////////////////////////////////////////////////
	load_I_TILDE_DFL(inst,true,I_TILDE_BEN_2_DFL,_cut_DFL_BEN_2_Y);
	comb_solve_model_DFL_BEN_2_AUX_1(inst,I_TILDE_BEN_2_DFL,AUX_SOL_BEN_2_DFL);
	////////////////////////////////////////////////////////////

#ifdef DEBUG

	cout << "\n\n******************************************************\n\n";

	cout << _BEST_INTEGER << " _BEST_INTEGER " << _BEST_REMAINING  << " _BEST_REMAINING " << " _NODE " << _NODE <<  endl;

	cout << "_cut_DFL_BEN_2_obj_value\t" << _cut_DFL_BEN_2_obj_value << endl;

	printf("\n\nsolution_DFL_BEN_2\n");
	for (int i = 0; i < inst->n_locations; i++)
	{
		printf("location %d -> %f\n",i,_cut_DFL_BEN_2_Y[i]);
	}


	for ( int i = 0; i < inst->n_clients; i++ )
	{
		cout << i <<	"\tpi\t" << AUX_SOL_BEN_2_DFL[i] << "\tsigma\t" << AUX_SOL_BEN_2_DFL[i+inst->n_clients] << endl;
	}
	cin.get();
#endif


	int rcnt=1;
	int nzcnt=inst->n_locations+1;
	double first_part=0;


	_cut_DFL_BEN_2_RHS=0.0;
	//sum of the sigma variable values
	for ( int j = 0; j < inst->n_clients; j++)
	{
		_cut_DFL_BEN_2_RHS +=	AUX_SOL_BEN_2_DFL[inst->n_clients+j];
	}



	for (int i = 0; i <inst->n_locations; i++)
	{

		_cut_DFL_BEN_2_rmatval[i]=0.0;

		for (int  k = inst->NFS[i]; k < inst->NFS[i+1]; k++ )
		{

			_cut_DFL_BEN_2_rmatval[i]=_cut_DFL_BEN_2_rmatval[i] - AUX_SOL_BEN_2_DFL[inst->AFS[k]];

		}

		first_part = first_part - (_cut_DFL_BEN_2_rmatval[i] *_cut_DFL_BEN_2_Y[i]);

		_cut_DFL_BEN_2_rmatind[i]=i;
	}

	_cut_DFL_BEN_2_rmatval[inst->n_locations]=1.0;
	_cut_DFL_BEN_2_rmatind[inst->n_locations]=inst->n_locations;

	//	cout << "current teta\t" << _cut_DFL_BEN_2_Y[inst->n_locations] << endl;
	//    cout << "val check\t" << _cut_DFL_BEN_2_RHS + first_part << endl;
	//	cout << "RHS\t" << _cut_DFL_BEN_2_RHS << endl;
	//	cin.get();

	if(_cut_DFL_BEN_2_Y[inst->n_locations]>_cut_DFL_BEN_2_RHS + first_part + 0.0001)
	{

#ifdef DEBUG_plus
		cout << "FIND CUT\tRHS\t" << _cut_DFL_BEN_2_RHS << endl;
#endif

#ifdef DEBUG_plus
		cout << "_cut_DFL_BEN_2_rmatval\n";
		for(int j=0;j<inst->n_locations;j++){
			cout << _cut_DFL_BEN_2_rmatval[j] << "\t";
		}
		cout << endl;

		cout << "_cut_DFL_BEN_2_rmatind\n";
		for(int j=0;j<inst->n_locations;j++){
			cout << _cut_DFL_BEN_2_rmatind[j] << "\t";
		}
		cout << endl;
#endif

#ifdef local_cut_DFL_BEN_2
		status=CPXcutcallbackaddlocal (env,cbdata,wherefrom,nzcnt,_cut_DFL_BEN_2_RHS,'L',_cut_DFL_BEN_2_rmatind,_cut_DFL_BEN_2_rmatval);
		if(status!=0){
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
#endif

#ifdef global_cut_DFL_BEN_2
		status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,_cut_DFL_BEN_2_RHS,'L',_cut_DFL_BEN_2_rmatind,_cut_DFL_BEN_2_rmatval,0);
		if(status!=0){
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
#endif

#ifdef global_cut_DFL_BEN_2_purgeable
		status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,_cut_DFL_BEN_2_RHS,'L',_cut_DFL_BEN_2_rmatind,_cut_DFL_BEN_2_rmatval,1);
		if(status!=0){
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
#endif

		n_cuts_BEN_DFL_1++;

		//inst->number_of_cuts=inst->number_of_cuts+1;

		//		cout << "******************ADDED\n";

#ifdef write_prob_call_back
		double per_test;
		CPXLPptr lp_test;
		CPXgetcallbacknodelp(env,cbdata,wherefrom,&lp_test);
		CPXwriteprob(env,lp_test,"prob.lp",NULL);
		cout << "NODE LP WRITTEN\n\n\n";
		cin.get();
#endif

		(*useraction_p)=CPX_CALLBACK_SET;

	}




#ifdef DEBUG
	cout << "\n\n******************************************************\n\n";
#endif

	//cout << "DONE\n";

	return 0;
}




/*****************************************************************/
int CPXPUBLIC myusercutcallback_DFL_BEN_2(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
/*****************************************************************/
{

	(*useraction_p)=CPX_CALLBACK_DEFAULT;

	instance *inst=(instance *) cbhandle;

	int num_variables=inst->n_locations+1;
	int status;
	double nu;
	double val_comparison;

	status=CPXgetcallbacknodex(env,cbdata,wherefrom,_cut_DFL_BEN_2_Y,0,num_variables-1);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}

	status=CPXgetcallbacknodeobjval(env,cbdata,wherefrom,&_cut_DFL_BEN_2_obj_value);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}

	double _BEST_INTEGER;
	status=CPXgetcallbackinfo(env,cbdata,wherefrom,CPX_CALLBACK_INFO_BEST_INTEGER,&_BEST_INTEGER);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}
	double _BEST_REMAINING;
	status=CPXgetcallbackinfo(env,cbdata,wherefrom,CPX_CALLBACK_INFO_BEST_REMAINING,&_BEST_REMAINING);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}
	int _NODE;
	status=CPXgetcallbackinfo(env,cbdata,wherefrom,CPX_CALLBACK_INFO_NODE_COUNT,&_NODE);
	if(status!=0){
		printf("cannot get the x\n");
		exit(-1);
	}

	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	//cout << "NODES\t" << _NODE << endl;
	//	if(_NODE>0)
	//	{
	//		return 0;
	//	}
	//
	//	if(n_cuts_BEN_FRAC_DFL_1>=500)
	//	{
	//		return 0;
	//	}
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////
	load_I_TILDE_DFL(inst,false,I_TILDE_BEN_2_DFL,_cut_DFL_BEN_2_Y);
	comb_solve_model_DFL_BEN_2_AUX_1(inst,I_TILDE_BEN_2_DFL,AUX_SOL_BEN_2_DFL);
	////////////////////////////////////////////////////////////


#ifdef DEBUG

	cout << "\n\n******************************************************\n\n";

	cout << _BEST_INTEGER << " _BEST_INTEGER " << _BEST_REMAINING  << " _BEST_REMAINING " << " _NODE " << _NODE <<  endl;

	cout << "_cut_DFL_BEN_2_obj_value\t" << _cut_DFL_BEN_2_obj_value << endl;

	printf("\n\nsolution_DFL_BEN_2\n");
	for (int i = 0; i < inst->n_locations; i++)
	{
		printf("location %d -> %f\n",i,_cut_DFL_BEN_2_Y[i]);
	}


	for ( int i = 0; i < inst->n_clients; i++ )
	{
		cout << i <<	"\tpi\t" << AUX_SOL_BEN_2_DFL[i] << "\tsigma\t" << AUX_SOL_BEN_2_DFL[i+inst->n_clients] << endl;
	}
	cin.get();
#endif


	int rcnt=1;
	int nzcnt=inst->n_locations+1;


	_cut_DFL_BEN_2_RHS=0.0;
	//sum of the sigma variable values
	for ( int j = 0; j < inst->n_clients; j++)
	{
		_cut_DFL_BEN_2_RHS+=	AUX_SOL_BEN_2_DFL[inst->n_clients+j];
	}


	double first_part=0;
	for (int i = 0; i <inst->n_locations; i++)
	{

		_cut_DFL_BEN_2_rmatval[i]=0.0;

		for (int  k = inst->NFS[i]; k < inst->NFS[i+1]; k++ )
		{

			_cut_DFL_BEN_2_rmatval[i]=_cut_DFL_BEN_2_rmatval[i] - AUX_SOL_BEN_2_DFL[inst->AFS[k]];

		}

		first_part = first_part - (_cut_DFL_BEN_2_rmatval[i] *_cut_DFL_BEN_2_Y[i] );

		_cut_DFL_BEN_2_rmatind[i]=i;
	}

	_cut_DFL_BEN_2_rmatval[inst->n_locations]=1.0;
	_cut_DFL_BEN_2_rmatind[inst->n_locations]=inst->n_locations;

	//	cout << "current teta\t" << _cut_DFL_BEN_2_Y[inst->n_locations] << endl;
	//	cout << "val check\t" << _cut_DFL_BEN_2_RHS + first_part << endl;
	//	cout << "RHS\t" << _cut_DFL_BEN_2_RHS << endl;
	//	cin.get();

	//	cout << "\nVIOL:\t"<<_cut_DFL_BEN_2_Y[inst->n_locations]- _cut_DFL_BEN_2_RHS + first_part << endl;

	double viol= (_cut_DFL_BEN_2_Y[inst->n_locations] - (_cut_DFL_BEN_2_RHS + first_part))/_cut_DFL_BEN_2_Y[inst->n_locations]*100 ;
	viol=max(viol,0.0001);
	//	cout << "viol\t" << viol << "\t"<<_cut_DFL_BEN_2_Y[inst->n_locations] - (_cut_DFL_BEN_2_RHS + first_part) << endl;
	//	cin.get();


	if(viol>1.0)
	{
		//	if(_cut_DFL_BEN_2_Y[inst->n_locations]>_cut_DFL_BEN_2_RHS + first_part + TOLL_FRAC_BEN)
		//	{

#ifdef DEBUG_plus
		cout << "FIND CUT\tRHS\t" << _cut_DFL_BEN_2_RHS << endl;
#endif

#ifdef DEBUG_plus
		cout << "_cut_DFL_BEN_2_rmatval\n";
		for(int j=0;j<inst->n_locations;j++){
			cout << _cut_DFL_BEN_2_rmatval[j] << "\t";
		}
		cout << endl;

		cout << "_cut_DFL_BEN_2_rmatind\n";
		for(int j=0;j<inst->n_locations;j++){
			cout << _cut_DFL_BEN_2_rmatind[j] << "\t";
		}
		cout << endl;
#endif

#ifdef local_cut_DFL_BEN_2
		status=CPXcutcallbackaddlocal (env,cbdata,wherefrom,nzcnt,_cut_DFL_BEN_2_RHS,'L',_cut_DFL_BEN_2_rmatind,_cut_DFL_BEN_2_rmatval);
		if(status!=0){
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
#endif

#ifdef global_cut_DFL_BEN_2
		status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,_cut_DFL_BEN_2_RHS,'L',_cut_DFL_BEN_2_rmatind,_cut_DFL_BEN_2_rmatval,0);
		if(status!=0){
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
#endif

#ifdef global_cut_DFL_BEN_2_purgeable
		status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,_cut_DFL_BEN_2_RHS,'L',_cut_DFL_BEN_2_rmatind,_cut_DFL_BEN_2_rmatval,1);
		if(status!=0){
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
#endif


		//inst->number_of_cuts=inst->number_of_cuts+1;

		n_cuts_BEN_FRAC_DFL_1++;

		//cout << "ADDED\n";

#ifdef write_prob_call_back
		double per_test;
		CPXLPptr lp_test;
		CPXgetcallbacknodelp(env,cbdata,wherefrom,&lp_test);
		CPXwriteprob(env,lp_test,"prob.lp",NULL);
		cout << "NODE LP WRITTEN\n\n\n";
		cin.get();
#endif

		(*useraction_p)=CPX_CALLBACK_SET;

	}


#ifdef DEBUG
	cout << "\n\n******************************************************\n\n";
#endif

	//cout << "DONE\n";

	return 0;
}



/*****************************************************************/
void build_model_DFL_BEN_2(instance *inst)
/*****************************************************************/
{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	inst->env_DFL_BEN_2=CPXopenCPLEX(&(inst->status));
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	inst->lp_DFL_BEN_2=CPXcreateprob(inst->env_DFL_BEN_2,&(inst->status),"DFL_BEN_2");
	if(inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables *
	inst->ccnt=inst->n_locations+1;
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
		//cout << "POSITION Y\t" << position_y_DFL_BEN_2(inst,j) << endl;
		sprintf(inst->colname[counter], "y%d",j);
		counter++;

	}

	inst->obj[counter]=1;
	inst->lb[counter]=0.0;


	double total_demand=0;

	for(int i=0;i<inst->n_clients;i++){
		total_demand+=inst->demand[i];
	}

	inst->ub[counter]=total_demand;
	//inst->ub[counter]=CPX_INFBOUND;


	inst->c_type[counter]='C';
	//cout << "POSITION Y\t" << position_y_DFL_BEN_2(inst,j) << endl;
	sprintf(inst->colname[counter], "tt");
	counter++;


	inst->status=CPXnewcols(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->c_type,inst->colname);
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
	CPXchgobjsen(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,CPX_MAX);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	// * creating the covering  constraint at least one facility *
	//	cout << "******ADDING AT CONSTRAINTS AT LEAST ONE LOCATION\n\n";
	//	inst->rcnt=1;
	//	inst->nzcnt=inst->n_locations;
	//
	//	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	//	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));
	//
	//	inst->rhs[0]=1.0;
	//	inst->sense[0]='G';
	//
	//
	//	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	//	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	//	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));
	//
	//	for ( int i = 0; i < inst->n_locations; i++ )
	//	{
	//		inst->rmatval[i]=1.0;
	//		inst->rmatind[i]=i;
	//	}
	//
	//	inst->rmatbeg[0]=0;
	//
	//	inst->status=CPXaddrows(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
	//	if(inst->status!=0)
	//	{
	//		printf("error in CPXaddrows\n");
	//		exit(-1);
	//	}
	//
	//	free(inst->rmatbeg);
	//	free(inst->rmatval);
	//	free(inst->rmatind);
	//	free(inst->rhs);
	//	free(inst->sense);


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
		inst->rmatind[i]=i;

	}

	inst->rmatbeg[0]=0;

	inst->status=CPXaddrows(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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


	AUX_SOL_BEN_2_DFL=(double*) calloc(2*inst->n_clients+1,sizeof(double));
	I_TILDE_BEN_2_DFL=(double*) calloc(inst->n_clients,sizeof(double));



	cout << "****compute compute_super_rho_DFL\n";
	compute_super_rho_DFL(inst,_cut_DFL_BEN_2_super_rho,_cut_DFL_BEN_2_SET);
	cout << "****compute compute_single_rho_DFL\n";
	compute_single_rho_DFL(inst,_cut_DFL_BEN_2_single_rho,_cut_DFL_BEN_2_SET);



#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status=CPXwriteprob(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,"DFL_BEN_2.lp",NULL);
	if(inst->status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	cout << "write DFL_BEN_2.lp\n\n";
	exit(-1);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

}





/*****************************************************************/
void solve_model_DFL_BEN_2(instance *inst)
/*****************************************************************/
{



	CPXsetintparam (inst->env_DFL_BEN_2, CPX_PARAM_SCRIND, CPX_ON);


	//	// * Set relative tolerance *
	//	inst->status = CPXsetdblparam (inst->env_DFL_BEN_2, CPX_PARAM_EPAGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}
	//
	//	// * Set a tolerance *
	//	inst->status = CPXsetdblparam (inst->env_DFL_BEN_2, CPX_PARAM_EPGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPGAP\n");
	//	}
	//
	//	// * Set mip tolerances integrality *
	//	inst->status = CPXsetdblparam (inst->env_DFL_BEN_2, CPX_PARAM_EPINT, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//
	//	// * Set Feasibility tolerance *
	//	inst->status = CPXsetdblparam (inst->env_DFL_BEN_2, CPX_PARAM_EPRHS, 1e-9);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPRHS\n");
	//	}

	// * Set number of CPU*
	inst->status = CPXsetintparam (inst->env_DFL_BEN_2, CPX_PARAM_THREADS, inst->number_of_CPU);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}

	// * Set time limit *
	inst->status = CPXsetdblparam (inst->env_DFL_BEN_2, CPX_PARAM_TILIM,inst->timelimit);
	if (inst->status)
	{
		printf ("error for CPX_PARAM_EPRHS\n");
	}



	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();


	CPXsetintparam(inst->env_DFL_BEN_2, CPX_PARAM_MIPCBREDLP, CPX_OFF);        // let MIP callbacks work on the original model
	CPXsetintparam(inst->env_DFL_BEN_2, CPX_PARAM_PRELINEAR, CPX_OFF);              // assure linear mappings between the presolved and original models
	CPXsetintparam(inst->env_DFL_BEN_2, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);


	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_DFL_BEN_2,mycutcallback_DFL_BEN_2,inst);
	if (inst->status)
	{
		printf ("error for CPXsetlazyconstraintcallbackfunc\n");
	}

	//	if(inst->option>1)
	//	{
	cout << "******FRACTIONAL SEPARATION******\n\n";
	inst->status = CPXsetusercutcallbackfunc(inst->env_DFL_BEN_2,myusercutcallback_DFL_BEN_2,inst);
	if (inst->status)
	{
		printf ("error for CPXsetuserconstraintcallbackfunc\n");
	}
	//	}

	cout << "\nCPXmipopt:\n";
	inst->status=CPXmipopt(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2);
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
	inst->x=(double*) calloc(inst->n_locations,sizeof(double));


	inst->status=CPXgetmipx(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,inst->x,0,inst->n_locations-1);
	if(inst->status!=0)
	{
		sol_found=false;
		printf("error in CPXgetmipx\n");
	}

	inst->objval=-1;
	inst->status=CPXgetmipobjval(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,&(inst->objval));
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

		for ( int j = 0; j < inst->n_locations; j++)
		{

			if( (int)(inst->x[j]+0.5) ==1)
			{
				open_facilities++;
			}
		}
	}

#ifdef print_solution
	printf("\n\nSolution\n");
	for ( int j = 0; j < inst->n_locations; j++){

		cout << "Location\t" << j << "\tval:\t" << (int)(inst->x[position_y_DFL_BEN_2(inst,j)]+0.5) << endl;

	}
	printf("\n");

	printf("\n");
#endif


	inst->bestobjval=-1;
	inst->status=CPXgetbestobjval(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2,&(inst->bestobjval));
	if(inst->status!=0)
	{
		printf("error in CPXgetbestobjval\n");
	}

	inst->lpstat=CPXgetstat(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2);
	inst->nodecount = CPXgetnodecnt(inst->env_DFL_BEN_2, inst->lp_DFL_BEN_2);

	cout << "\n\nlpstat\t" << inst->lpstat << endl;

	///////////////////////////////////////////////////////////////////////////////////////////
	bool *reached_clients=new bool[inst->n_clients];
	satisfied_clients=compute_subset_location_coverage_clients(inst,inst->x,reached_clients);
	delete[] reached_clients;
	///////////////////////////////////////////////////////////////////////////////////////////

//	///////////////////////////////////////////////////////////////////////////////////////////
//	bool *reached_clients=new bool[inst->n_clients];
//	satisfied_clients=compute_subset_location_coverage_clients(inst,inst->x,reached_clients);
//
//	if(inst->cohordinates_loaded==true){
//		int *zz=new int[inst->n_clients];
//		int *yy=new int[inst->n_locations];
//		for ( int j = 0; j < inst->n_locations; j++)
//		{
//			yy[j]= (int)(inst->x[j]+0.5);
//		}
//		for ( int i = 0; i < inst->n_clients; i++ ){
//			if(reached_clients[i])
//			{
//				zz[i]=1;
//			}
//			else
//			{
//				zz[i]=0;
//			}
//		}
//
//		//draw_grid_sol(inst,yy,zz);
//		delete[]zz;
//		delete[]yy;
//	}
//	delete[] reached_clients;
//	///////////////////////////////////////////////////////////////////////////////////////////

	cout << "\n***open_facilities\t" << open_facilities << endl;
	cout << "***satisfied_clients\t" << satisfied_clients << endl;


	cout << "\n\nSTAT:\tobjval\t" << inst->objval << "\tbestobjval\t" << inst->bestobjval << "\tlpstat\t" << inst->lpstat << "\topen_facilities\t" << open_facilities << "\tsatisfied_clients\t" << satisfied_clients << "\ttime\t"<< solution_time<< endl << endl;




	int cur_numcols=CPXgetnumcols(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2);
	int cur_numrows=CPXgetnumrows(inst->env_DFL_BEN_2,inst->lp_DFL_BEN_2);

	cout <<"n_cuts_BEN_DFL_1\t" <<  n_cuts_BEN_DFL_1 << "\n";
	cout <<"n_cuts_BEN_FRAC_DFL_1\t"<<  n_cuts_BEN_FRAC_DFL_1 << "\n";
	cout <<"n_cuts_BEN_DFL_2\t" <<  n_cuts_BEN_DFL_2 << "\n";
	cout <<"n_cuts_BEN_FRAC_DFL_2\t"<<  n_cuts_BEN_FRAC_DFL_2 << "\n";

	//	printf("\nnumcols\t%d\n",cur_numcols);
	//	printf("\nnumrows\t%d\n",cur_numrows);

	//	ofstream compact_file;
	//	compact_file.open("info_DFL_BEN_2.txt", ios::app);
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
	////			<<  inst->option<< "\t"
	//			<<  inst->seed<< "\t"
	//			<<  inst->client_not_covered<< "\t"
	//			<<  inst->client_single_covered<< "\t"
	//			<<  n_cuts_BEN_DFL_1 << "\t"
	//			<<  n_cuts_BEN_FRAC_DFL_1 << "\t"
	//			<<  n_cuts_BEN_DFL_2 << "\t"
	//			<<  n_cuts_BEN_FRAC_DFL_2 << "\t"
	//			<<  inst->counter_c << "\t"
	//			<< endl;
	//	compact_file.close();



	free(inst->x);

}


/*****************************************************************/
void clean_model_DFL_BEN_2(instance *inst)
/*****************************************************************/
{

	inst->status=CPXfreeprob(inst->env_DFL_BEN_2,&(inst->lp_DFL_BEN_2));
	if(inst->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	inst->status=CPXcloseCPLEX(&(inst->env_DFL_BEN_2));
	if(inst->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

	free(AUX_SOL_BEN_2_DFL);
	free(I_TILDE_BEN_2_DFL);

}
