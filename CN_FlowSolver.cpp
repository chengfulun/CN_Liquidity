#include "CN_CplexSolver.h"



#include "ilcplex/cplex.h"
#include <stdlib.h>
#include <string.h>

static int
buildNetwork (CPXENVptr env, CPXNETptr net, WidgetGraph* widgetNet, int nnodes, int narcs,
	double * &supply, int* &head, int* &tail, double* &obj, double* &ub, double* &lb);

int CplexSolver::solve(CplexConverter& converter)
{
	/* Declare variables and arrays for retrieving problem data and
	  solution information later on. */
	// cout << "lp solver" << endl;	
	int      narcs;
	int      nnodes;
	int      solstat;
	double   objval;

	CPXENVptr env = NULL;
	CPXNETptr net = NULL;
	int       status;
	int status1;
	int       i, j;
	int cnt = 0;

	/* Initialize the CPLEX environment */

	env = CPXopenCPLEX (&status);

	/* If an error occurs, the status value indicates the reason for
	  failure.  A call to CPXgeterrorstring will produce the text of
	  the error message.  Note that CPXopenCPLEX produces no
	  output, so the only way to see the cause of the error is to use
	  CPXgeterrorstring.  For other CPLEX routines, the errors will
	  be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

	if ( env == NULL ) {
	  char  errmsg[CPXMESSAGEBUFSIZE];
	  fprintf (stderr, "Could not open CPLEX environment.\n");
	  CPXgeterrorstring (env, status, errmsg);
	  fprintf (stderr, "%s", errmsg);
	  goto TERMINATE;
	}

	/* Turn on output to the screen */

	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( status ) {
	  fprintf (stderr, 
				"Failure to turn on screen indicator, error %d.\n", status);
	  goto TERMINATE;
	}

	/* Create the problem. */

	net = CPXNETcreateprob (env, &status, "netex1");

	/* A returned pointer of NULL may mean that not enough memory
	  was available or there was some other problem.  In the case of 
	  failure, an error message will have been written to the error 
	  channel from inside CPLEX.  In this example, the setting of
	  the parameter CPXPARAM_ScreenOutput causes the error message to
	  appear on stdout.  */

	if ( net == NULL ) {
	  fprintf (stderr, "Failed to create network object.\n");
	  goto TERMINATE;
	}

	/* Fill in the data for the problem.  Note that since the space for
	  the data already exists in local variables, we pass the arrays
	  directly to the routine to fill in the data structures.  */

	status = buildNetwork (env, net, converter.widgetNet, 
		converter.nnodes, converter.narcs,
		converter.supply, converter.head, converter.tail, 
		converter.obj, converter.ub, converter.lb);

	// if(supply != NULL){
	// 	cout << "before not clear" << endl;
	// }

	if ( status ) {
	  fprintf (stderr, "Failed to build network problem.\n");
	  goto TERMINATE;
	}


	/* Optimize the problem and obtain solution. */

	status = CPXNETprimopt (env, net);
	if ( status ) {
	  fprintf (stderr, "Failed to optimize network.\n");
	  goto TERMINATE;
	}

	/* get network dimensions */

	narcs  = CPXNETgetnumarcs  (env, net);
	nnodes = CPXNETgetnumnodes (env, net);

	status = CPXNETsolution (env, net, &solstat, &objval, 
		converter.x, converter.pi, converter.slack, converter.dj);
	// cout << "status: " << status << endl;
	if ( status ) {
	  // fprintf (stderr, "Failed to obtain solution.\n");
	  goto TERMINATE;
	}

	/* Write the output to the screen. */

	// printf ("\nSolution status = %d\n", solstat);
	if (solstat != 1 && solstat != 6 && solstat != 14){
		// cout << "status" << solstat << endl; 
		status1 = -1;
		goto TERMINATE;
	}
	// printf ("Solution value  = %f\n\n", objval);

	status1 = 0;
	
	
TERMINATE:

	/* Free up the problem as allocated by CPXNETcreateprob, if necessary */

	if ( net != NULL ) {
	  status = CPXNETfreeprob (env, &net);
	  if ( status ) {
		 fprintf (stderr, "CPXNETfreeprob failed, error code %d.\n", status);
	  }
	}

	if ( env != NULL ) {
	  status = CPXcloseCPLEX (&env);

	  if ( status ) {
	  char  errmsg[CPXMESSAGEBUFSIZE];
		 fprintf (stderr, "Could not close CPLEX environment.\n");
		 CPXgeterrorstring (env, status, errmsg);
		 fprintf (stderr, "%s", errmsg);
	  }
	}

	return (status1);

}



static int
buildNetwork (CPXENVptr env, CPXNETptr net, WidgetGraph* widgetNet, int nnodes, int narcs,
	double * &supply, int* &head, int* &tail, double* &obj, double* &ub, double* &lb)
{
	int status = 0;

	/* definitions to improve readability */

#  define inf    CPX_INFBOUND

	if ( CPXNETgetnumnodes (env, net) > 0 ) {
	  status = CPXNETdelnodes (env, net, 0,
								CPXNETgetnumnodes (env, net)-1);
	  if ( status ) goto TERMINATE;
	}

	/* Set optimization sense */

	status = CPXNETchgobjsen (env, net, CPX_MAX);
	if ( status ) goto TERMINATE;

	/* Add nodes to network along with their supply values,
	  but without any names. */

	status = CPXNETaddnodes (env, net, nnodes, supply, NULL);
	if ( status ) goto TERMINATE;

	/* Add arcs to network along with their objective values and
	  bounds, but without any names. */

	status = CPXNETaddarcs (env, net, narcs, tail, head, lb, ub, obj, NULL);
	if ( status ) goto TERMINATE;

	// if (supply != NULL){
	// 	cout << "before exit build, supply not empty" << endl;
	// }

TERMINATE:

	return (status);

}
