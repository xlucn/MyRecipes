/**
 * @file ODEsol.c
 * @brief ODEsol.c
 */

#include <stdlib.h>
#include "NR.h"

/**
 * @brief a struct type to contain the solution of a system of ODEs.
 */
struct _SODEsol{
	int step; /**< the steps used in the integration. It is useful in variable step size methods.*/
	double *t; /**< the pointer to the independent virable list */
	double **y; /**< the pointer to the dependent variables list */
};

/* the typedef following is in NR.h, but moving the comment here out of
 simplicity of the header file */
/**
 * @typedef typedef struct _SODEsol *SODEsol
 * @brief pointer to struct type containing solution of an SODE
 */

/**
 * @brief a struct type to contain the solution of an ODE.
 */
struct _ODEsol{
	int step; /**< the steps used in the integration. It is useful in variable step size methods. */
	double *t; /**< the pointer to the independent virable list */
	double *y; /**< the pointer to the dependent variable list */
};

/* the typedef following is in NR.h, but moving the comment here out of
 simplicity of the header file */
/**
 * @typedef typedef struct _ODEsol *ODEsol
 * @brief pointer to struct type containing solution of an ODE
 */

/**
 * @brief create a variable of type SODEsol
 * @param step step
 * @param t variable values
 * @param y integrals for every variable
 * @returns an SODEsol variable
 */
SODEsol newSODEsol(int step, double *t, double **y)
{
	SODEsol s = (SODEsol)malloc(sizeof(struct _SODEsol));
	s->step = step;
	s->t = t;
	s->y = y;
	return s;
}

/**
 * @brief create a variable of type ODEsol
 * @param step step
 * @param t variable values
 * @param y integrals for every variable
 * @returns an ODEsol variable
 */
ODEsol newODEsol(int step, double *t, double *y)
{
	ODEsol s = (ODEsol)malloc(sizeof(struct _ODEsol));
	s->step = step;
	s->t = t;
	s->y = y;
	return s;
}

/**
 * @brief get the step member of an ODEsol variable
 * @param sol an ODEsol variable
 * @returns the step member of sol
 */
int ODEsolGetStep(ODEsol sol)
{
	return sol->step;
}

/**
 * @brief get the t member of an ODEsol variable
 * @param sol an ODEsol variable
 * @returns the t member of sol
 */
double *ODEsolGetT(ODEsol sol)
{
	return sol->t;
}

/**
 * @brief get the y member of an ODEsol variable
 * @param sol an ODEsol variable
 * @returns the y member of sol
 */
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

/**
 * @brief get the step member of an SODEsol variable
 * @param sol an SODEsol variable
 * @returns the step member of sol
 */
int SODEsolGetStep(SODEsol sol)
{
	return sol->step;
}

/**
 * @brief get the t member of an SODEsol variable
 * @param sol an SODEsol variable
 * @returns the t member of sol
 */
double *SODEsolGetT(SODEsol sol)
{
	return sol->t;
}

/**
 * @brief get the y member of an SODEsol variable
 * @param sol an SODEsol variable
 * @returns the y member of sol
 */
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


