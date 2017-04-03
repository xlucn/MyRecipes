#include "NR.h"
/**
 * @brief a struct type to contain the solution of a system of ODE.
 */
struct _SODEsol{
	int step; /**< the steps used in the integration. It is useful in variable step size methods.*/
	double *t; /**< the pointer to the independent virable list */
	double **y; /**< the pointer to the dependent variables list */
};

/**
 * @brief a struct type to contain the solution of an ODE.
 */
struct _ODEsol{
	int step; /**< the steps used in the integration. It is useful in variable step size methods. */
	double *t; /**< the pointer to the independent virable list */
	double *y; /**< the pointer to the dependent variable list */
};

SODEsol newSODEsol(int step, double *t, double **y)
{
	SODEsol s = (SODEsol)malloc_s(sizeof(struct _SODEsol));
	s->step = step;
	s->t = t;
	s->y = y;
	return s;
}

ODEsol newODEsol(int step, double *t, double *y)
{
	ODEsol s = (ODEsol)malloc_s(sizeof(struct _ODEsol));
	s->step = step;
	s->t = t;
	s->y = y;
	return s;
}

int ODEsolGetStep(ODEsol sol)
{
	return sol->step;
}

double *ODEsolGetT(ODEsol sol)
{
	return sol->t;
}

double *ODEsolGetY(ODEsol sol)
{
	return sol->y;
}

/**
 * @brief free a struct type ODEsol
 * @param sol a struct type ODEsol
 */
void delODEsol(ODEsol sol)
{
	free(sol->t);
	free(sol->y);
	free(sol);
}

int SODEsolGetStep(SODEsol sol)
{
	return sol->step;
}

double *SODEsolGetT(SODEsol sol)
{
	return sol->t;
}

double **SODEsolGetY(SODEsol sol)
{
	return sol->y;
}

/**
 * @brief free a struct type SODEsol
 * @param sol a struct type SODEsol
 */
void delSODEsol(SODEsol sol)
{
	free(sol->t);
	for (int i = 0; i < sol->step; i++)
	{
		free(sol->y[i]);
	}
	free(sol->y);
	free(sol);
}


